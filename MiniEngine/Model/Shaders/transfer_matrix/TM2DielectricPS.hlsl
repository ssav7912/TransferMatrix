#include "HenyeyGreenstein.hlsli"
#include "tm2.hlsli"
#include "MatrixOperators.hlsli"
#include "../Common.hlsli"

#define LAYERS_MAX 5
#define NUM_LAYERS 2
#define NUM_LOBES (NUM_LAYERS + 1)

#define NO_SECOND_UV 1

//Lookup table for Total Internal Reflection
Texture3D<float3> TIR_LUT : register(t14);

//LUT for Karis Split-sum FGD approximation.
Texture2D<float2> FGD_LUT : register(t15);

SamplerState TIR_Sampler : register(s0);


cbuffer GlobalConstants : register(b1)
{
    float4x4 ViewProj;
    float4x4 SunShadowMatrix;
    float3 ViewerPos;
    float3 SunDirection;
    float3 SunIntensity;
    float _pad;
    float IBLRange;
    float IBLBias;
}


struct VSOutput
{
    float4 position : SV_POSITION;
    float3 normal : NORMAL;
#ifndef NO_TANGENT_FRAME
    float4 tangent : TANGENT;
#endif
    float2 uv0 : TEXCOORD0;
#ifndef NO_SECOND_UV
    float2 uv1 : TEXCOORD1;
#endif
    float3 worldPos : TEXCOORD2;
    float3 sunShadowCoord : TEXCOORD3;
};

//PRNG hack for sampling
float nrand(float2 uv)
{
    return frac(sin(dot(uv, float2(12.9898, 78.233))) * 43758.5453);
}


//TODO: Probably drop TIR for being too expensive.
float3 TIR_lookup(float3 coords)
{
    return TIR_LUT.Sample(TIR_Sampler, coords);
}

//Lagarde 2011. Compute fresnel reflectance at 0 degrees from IOR
float f0(float ior)
{
    return ((ior - 1) * (ior - 1)) / ((ior + 1) * (ior + 1));

}

void albedos(float cti, float alpha, float ior_ij, out float3 r_ij, out float3 t_ij, out float3 r_ji, out float3 t_ji)
{
    if (abs(ior_ij - 1.0f) < 1e-3f)
    {
        r_ij = 0.0f.xxx;
        r_ji = 0.0f.xxx;
        
        t_ij = 1.0f.xxx;
        t_ji = 1.0f.xxx;
        return;
    }
    
    //FGD 4D LUT parameterised by elevation, roughness and complex IOR
    //(Karis 2013) approximation parameterised by Roughness and angle
    //outputs scale & bias to F0.
    //TODO: work out 2D split sum approx.
    //https://cdn2.unrealengine.com/Resources/files/2013SiggraphPresentationsNotes-26915738.pdf
    //Ok split sum approx contains FGD_1 & FGD_2 (whatever that means)
    //tristimulus energy output?
    float2 splitsum = FGD_LUT.Sample(TIR_Sampler, float2(alpha, cti));
    
    r_ij = splitsum.xxx;
    r_ji = r_ij;
    
    t_ij = 1.0f.xxx - r_ij;
    t_ji = t_ij;
    
}

void dielectric_transfer_factors(float3 incident, float ior, float alpha, out layer_components_tm2 ops)
{
    float ior_ji = 0.0f;
    float s_t_ij = 0.0f;
    float s_t_ji = 0.0f;
    
    if (abs(ior - 1.0f) < 1e-5f)
    {
        ops.component_type = TM2_TYPE_NOCOMPONENT;
        ops.transmission_down.mean = -incident;
        
        return;
    }
    
    ops.component_type = TM2_TYPE_DIELECTRICINTERFACE;
    
    ior_ji = 1.0f / ior;
    
    ops.reflection_down.asymmetry = ggx_to_hg(alpha);
    ops.reflection_down.mean = reflectZ(incident);
    
    ops.reflection_up.asymmetry = ops.reflection_down.asymmetry;
    ops.reflection_up.mean = refractZ(incident, ior);

    s_t_ij = abs(ior_ji * ops.reflection_down.mean.z + ops.reflection_up.mean.z) / ops.reflection_up.mean.z;
    s_t_ji = abs(ior * ops.reflection_up.mean.z + ops.reflection_down.mean.z) / ops.reflection_down.mean.z;

    ops.transmission_down.asymmetry = ggx_to_hg(0.5f * s_t_ij * alpha);
    ops.transmission_down.mean = ops.reflection_up.mean;
    
    ops.transmission_up.asymmetry = ggx_to_hg(0.5f * s_t_ji * alpha);
    ops.transmission_up.mean = ops.reflection_down.mean;
    
    albedos(abs(incident.z), alpha, ior, ops.reflection_down.norm, ops.transmission_down.norm, ops.reflection_up.norm, ops.transmission_up.norm);

}


void components_transfer_factors(float3 incident, inout float3 iors[LAYERS_MAX], inout float alphas[LAYERS_MAX], out layer_components_tm2 ops[LAYERS_MAX])
{
    float3 ior_ij;
    
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        //TODO: FGD lookup
        dielectric_transfer_factors(incident, float3_average(iors[i + 1] / iors[i]), alphas[i + 1], ops[i]);
        incident -= ops[i].transmission_down.mean;
    }

}

void outgoing_lobes(float3 incident, float3 iors[LAYERS_MAX], float alphas[LAYERS_MAX], out henyey_greenstein lobes[LAYERS_MAX], out float3 lobe_incident[LAYERS_MAX])
{
    
    //interface attributes
    float ior_ij = 0;
    float ior_ji = 0;
    float3 e_r_i = 0;
    float e_r_i_avg = 0;
    
    //cumulated attributes
    float3 e_r_0h = 0.f.xxx;
    float3 e_r_0i = 0.f.xxx;
    float3 e_t_0i = 0.f.xxx;
    
    float e_t_0i_avg = 0.0f;
    float g_r_0h = 0.0f;
    float g_r_0i = 0.0f;
    float g_t_0i = 0.0f;
    
    //first order terms
    float g_T_0i = 1.0f; 
    float g_T_0j = 0.0f;
    float g_T_0j_R = 0.0f;
    float g_T_0j_RT = 0.0f;
    
    float3 tir_norm = 0.0f.xxx;
    
    //identity transfer matrices
    tensor3d2x2 etm_0i = {     1.0f.xxx, 0.0f.xxx,
                               0.0f.xxx, 1.0f.xxx 
        };
    
    float2x2 gtm_0i = float2x2(1, 0, 
                               0, 1);
    
    //TODO: calc transfer factors
    layer_components_tm2 ops[LAYERS_MAX];
    components_transfer_factors(incident, iors, alphas, ops);
   
    //component iterations
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        ior_ij = float3_average(iors[i + 1] / iors[i]);
        ior_ji = 1.0f / ior_ij;
        
        //transmission correction
        
        g_T_0j = hg_refract(g_T_0i, ior_ji) * ops[i].transmission_down.asymmetry;
        g_T_0j_R = i < (NUM_LAYERS - 1) ? g_T_0j * ops[i + 1].reflection_down.asymmetry : 1.0f;
        g_T_0j_RT = hg_refract(g_T_0j_R, ior_ij) * ops[i].transmission_up.asymmetry;
        
        ops[i].transmission_down.asymmetry = g_T_0i != 0.f ? g_T_0j / g_T_0i : 0.0f;
        
        ops[i].transmission_up.asymmetry = g_T_0j_R != 0.0f ? g_T_0j_RT / g_T_0j_R : 0.0f;
        
        if (ior_ij < 1.0f)
        {
            //TODO: TIR
            tir_norm = TIR_lookup(float3(abs(ops[i].reflection_down.mean.z), hg_to_ggx(g_T_0i), ior_ij)) * ops[i].transmission_down.norm;

            ops[i].reflection_down.norm += tir_norm;
            ops[i].transmission_down.norm -= tir_norm;
    
        }
        else
        {
            tir_norm = TIR_lookup(float3(abs(ops[i].transmission_down.mean.z), hg_to_ggx(g_T_0j_R), ior_ji)) * ops[i].transmission_up.norm;
            
            ops[i].reflection_up.norm += tir_norm;
            ops[i].transmission_up.norm -= tir_norm;

        }
        
        //transfer matrix eval
        etm_0i = mul(etm_0i, energy_matrix(ops[i]));
        gtm_0i = mul(gtm_0i, asymmetry_matrix(ops[i]));
        
        //reflected lobe
        
        //Norm
        e_r_0i = reflection_energy_tm2(etm_0i);
        e_r_i = e_r_0i - e_r_0h;
        e_r_i_avg = float3_average(e_r_i);
        
        //asymmetry
        //g_r_0i = reflection_energy_tm2(gtm_0i);
        g_r_0i = gtm_0i._21 / gtm_0i._11;

        //reflected lobe
        lobes[i].norm = e_r_i;
        lobes[i].asymmetry = e_r_i_avg > 0.f ? min((g_r_0i - g_r_0h) / e_r_i_avg, 1.0f) : 0.0f;
        lobe_incident[i] = incident;
        
        //top layer attributes
        e_r_0h = e_r_0i;
        g_r_0h = g_r_0i;
        
        //first order terms
        g_T_0i = g_T_0j;
    }
    
    //transmitted lobe
    
    //norm
    e_t_0i = transmission_energy_tm2(etm_0i);
    e_t_0i_avg = float3_average(e_t_0i);
    
    //asymmetry
    //g_t_0i = transmission_energy_tm2(gtm_0i);
    g_t_0i = 1.0f / gtm_0i._11;
    
    lobes[NUM_LAYERS].norm = e_t_0i;
    lobes[NUM_LAYERS].asymmetry = e_t_0i_avg > 0.0f ? min(g_t_0i / e_t_0i_avg, 1.0f) : 0.0f;
    lobe_incident[NUM_LAYERS] = reflectZ(ops[NUM_LAYERS - 1].transmission_down.mean);

}

float3 eval(sample_record rec, int measure, float3 iors[LAYERS_MAX], float alphas[LAYERS_MAX])
{
    if (measure != TM2_MEASURE_SOLID_ANGLE || rec.incident.z > 0)
    {
        return 0.0f.xxx;
    }
    
    bool reflect = rec.incident.z * rec.outgoing.z > 0;
    
    /* OUTGOING LOBES */
    henyey_greenstein lobes[LAYERS_MAX];
    float3 lobes_incident[LAYERS_MAX];
    
    //TODO: outgoing_lobes
    outgoing_lobes(rec.incident, iors, alphas, lobes, lobes_incident); 
    
    
    /* THROUGHPUT */
    const float factor = reflect ? 1.f : float3_average((iors[0] / iors[NUM_LAYERS]));
    
    float3 throughput = 0.0f.xxx;
    
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        if (lobes[i].norm.x == 0.0f && lobes[i].norm.y == 0.0f && lobes[i].norm.z == 0.0f)
        {
            continue;
        }
        
        if (lobes_incident[i].z * rec.outgoing.z <= 0.0f)
        {
            continue;
        }
        
        float3 incoming = lobes_incident[i];
        float3 outgoing = rec.outgoing;
        
        incoming.z = abs(incoming.z);
        outgoing.z = abs(outgoing.z);
        
        const float3 H = normalize(incoming + outgoing);
        
        const float roughness = hg_to_ggx(lobes[i].asymmetry);
        //TODO: microfacet distribution
        const float G2 = smithG1(incoming, H, roughness) * smithG1(outgoing, H, roughness);
        const float D = D_GGX(H, roughness);
        
        const float f = G2 * D / (4.0f / incoming.z);
        
        
        throughput += lobes[i].norm * f * factor * factor;
    }
    
    return throughput;

}

float3 sample(inout sample_record rec, float2 samplePos, float3 iors[LAYERS_MAX], float alphas[LAYERS_MAX])
{
    /* LOBES */
	
    henyey_greenstein lobes[LAYERS_MAX];
    float3 lobes_incoming_direction[LAYERS_MAX];
	
	//TODO: Outgoing lobes
	
    float w[LAYERS_MAX];
    float w_sum = 0.0f;
	
    for (int i = 0; i < NUM_LOBES; i++)
    {
        w[i] = float3_average(lobes[0].norm);
        w_sum += w[i];

    }
	
    /* LOBE SELECTION */
    
	//TODO: random sampling?
    float sel_w = nrand(samplePos)* w_sum - w[0];
    int sel_i = 0;
    for (sel_i = 0; sel_w > 0.f && sel_i < NUM_LOBES; ++sel_i)
    {
        sel_w -= w[sel_i + 1];    
    }
    
    
    /* SAMPLING */
    
    const bool reflection = rec.incident.z * lobes_incoming_direction[sel_i].z > 0.f;
    
    //TODO: Microfacet distribution.
    const float sel_rough = hg_to_ggx(lobes[sel_i].asymmetry);
    
    
	
    //upper hemisphere reflection evaluation direction
    float3 sel_wi = lobes_incoming_direction[sel_i];
    sel_wi.z = abs(sel_wi.z);
   
   
    //normal sampling
    //TODO: microfacet NDF
    //TODO: use a sampleVisible rather than SampleAll distribution?
    float pdf;
    float3 H = sample_GGX(samplePos, sel_rough, pdf);
    
    if (pdf <= 0.0f)
    {
        return 0.0f.xxx;
    }
    
    //sampled direction
    float3 sel_wo = reflect(H, sel_wi);
    sel_wo.z *= sign(rec.incident.z) * (reflection ? 1.f : -1.f);
    
    //sampling record
    rec.outgoing = sel_wo;
    rec.ior = float3_average(iors[reflection ? 0 : NUM_LAYERS]);
    rec.is_reflection_sample = reflection ? false : true;
    rec.sample_type = reflection ? TM2_SAMPLE_TYPE_GLOSSY_REFLECTION : TM2_SAMPLE_TYPE_GLOSSY_TRANSMISSION;
    
    pdf = 0.0f;
    if (pdf <= 0.f)
    {
        return float3(0.0f, 0.0f, 0.0f);
    }
    
    /* PDF */
    for (int i = 0; i < NUM_LOBES; i++)
    {
        //likely to cause divergence
        //TODO: profile and refactor.
        if (w[i] == 0.f)
        {
            continue;
        }
        
        if (lobes_incoming_direction[i].z * rec.outgoing.z <= 0.f)
        {
            continue;
        }
        
        float3 incoming = lobes_incoming_direction[i];
        float3 outgoing = rec.outgoing;
        
        incoming.z = abs(incoming.z);
        outgoing.z = abs(outgoing.z);
        
        const float3 H = normalize(incoming + outgoing);
        
        const float rough = hg_to_ggx(lobes[i].asymmetry);
        //TODO: Microfacet NDF
        
        const float G1 = smithG1(incoming, H, rough);
        
        const float D = D_GGX(H, rough);
        
        pdf += w[i] * G1 * D / (4.0f * incoming.z);
        
    }
    
    pdf /= w_sum;
    
    /* THROUGHPUT */
    
    float3 throughput = eval(rec, TM2_MEASURE_SOLID_ANGLE, iors, alphas);
    
    
    return pdf > 0.f ? throughput / pdf : float3(0.0f, 0.0f, 0.0f);
    
}

//projects tangent space vector onto normal plane
float2 tangent_to_projective(float3 tangent_in)
{
    return float2(tangent_in.x, tangent_in.y);

}

float4 main(VSOutput vsOutput) : SV_TARGET
{
    //interface between air and arbitrarysurface.
    float3 iors[LAYERS_MAX] = { float3(1.3, 1.3, 1.3f), float3(1.4f, 1.4, 1.4f), 1.0.xxx, 1.0.xxx, 1.0.xxx };
    float alphas[LAYERS_MAX] = { 0.0f, 0.8f, 1.0f, 1.0f, 1.0f };
	
    
    float3 normal = normalize(vsOutput.normal);
    
    float3 tangent = normalize(vsOutput.tangent.xyz);
    float3 bitangent = normalize(cross(normal, tangent)) * vsOutput.tangent.w;
    float3x3 tangentFrame = float3x3(tangent, bitangent, normal);
    
    float3 SunDirectionTangent = mul(SunDirection, tangentFrame);
    
    //do 4 samples?
    sample_record rec = 
    {
        SunDirectionTangent,
        reflect(SunDirectionTangent, float3(0, 0, 1)),
        1.0f,
        true,
        TM2_SAMPLE_TYPE_GLOSSY_REFLECTION
    };
    
    float3 accum = sample(rec, float2(nrand(vsOutput.uv0), nrand(vsOutput.uv0)), iors, alphas);
    rec.sample_type = TM2_SAMPLE_TYPE_GLOSSY_TRANSMISSION;
    
    accum += sample(rec, float2(nrand(vsOutput.uv0), nrand(vsOutput.uv0)), iors, alphas);
    rec.sample_type = TM2_SAMPLE_TYPE_GLOSSY_REFLECTION;
    
    accum += sample(rec, float2(nrand(vsOutput.uv0), nrand(vsOutput.uv0)), iors, alphas);
    rec.sample_type = TM2_SAMPLE_TYPE_GLOSSY_TRANSMISSION;
    
	
	return float4(accum, 1.0f);
}
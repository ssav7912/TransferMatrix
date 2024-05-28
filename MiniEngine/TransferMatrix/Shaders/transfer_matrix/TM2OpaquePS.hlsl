#include "tm2_common.hlsli"
#include "MatrixOperators.hlsli"


void components_transfer_factors(float3 incident, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughness[LAYERS_MAX], out layer_components_tm2 ops[LAYERS_MAX])
{
    float3 ior_ij = 0.0f;
    
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        ior_ij = iors[i + 1] / iors[i];
        if (kappas[i+1].x + kappas[i+1].y + kappas[i+1].z == 0.0f)
        {
            dielectric_transfer_factors(incident, float3_average(ior_ij), roughness[i + 1], ops[i]);
            
            incident = -ops[i].transmission_down.mean;
        }
        else
        {
            conductor_transfer_factors(incident, ior_ij, kappas[i + 1] / iors[i], roughness[i + 1], ops[i]); 
        }
    }

}

void outgoing_lobes(float3 incident, float3 ior[LAYERS_MAX], float3 kappas[LAYERS_MAX], float3 roughness[LAYERS_MAX], out henyey_greenstein lobes[LAYERS_MAX])
{
    layer_components_tm2 ops[LAYERS_MAX]; 
    
    float ior_ij = 0.0f;
    float ior_ji = 0.0f;
    
    float3 energy_r_i = 0.0f;
    float energy_r_i_average = 0.0f;
    
    float3 energy_r_0h = 0.0f.xxx;
    float3 energy_r_0i = 0.0f.xxx;
    
    float asymmetry_r_0h = 0.0f;
    float asymmetry_r_0i = 0.0f;
    
    float asymmetry_T_0i = 1.0f;
    float asymmetry_T_0j = 1.0f;
    float asymmetry_T_0j_R = 1.0f;
    float asymmetry_T_0j_RT = 1.0f;
    
    float3 tir_norm;
    
    tensor3d2x2 energy_0i =
    {
        1.0f.xxx, 0.0f.xxx,
                               0.0f.xxx, 1.0f.xxx 
    };
    
    float2x2 asymmetry_0i = float2x2(1, 0,
                               0, 1);
    
    
    components_transfer_factors(incident, iors, kappas, alphas, ops);
    
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        if (ops[i].component_type == TM2_TYPE_DIELECTRICINTERFACE)
        {
            ior_ij = float3_average((ior[i + 1] / ior[i]));
            ior_ji = 1.0f / ior_ij;
            
            asymmetry_T_0j = hg_refract(asymmetry_T_0i, ior_ji) * ops[i].transmission_down.asymmetry;
            asymmetry_T_0j_R = asymmetry_T_0j * ops[i + 1].reflection_down.asymmetry;
            asymmetry_T_0j_RT = hg_refract(asymmetry_T_0j_R, ior_ij) * ops[i].transmission_up.asymmetry;
            
            ops[i].transmission_down.asymmetry = asymmetry_T_0i != 0.0f ? asymmetry_T_0j / asymmetry_T_0i : 0.0f;
            
            ops[i].transmission_up.asymmetry = asymmetry_T_0j_R != 0.0f ? asymmetry_T_0j_RT / asymmetry_T_0j_R : 0.0f;
            
            if (ior_ij < 1.0f)
            {
                tir_norm = TIR_lookup(abs(ops[i].reflection_down.mean.z), hg_to_ggx(asymmetry_T_0i), ior_ij) * ops[i].transmission_down.norm;
            
                ops[i].reflection_down.norm += tir_norm;
                ops[i].transmission_down.norm += tir_norm;
            }
            else
            {
                tir_norm = TIR_lookup(abs(ops[i].transmission_down.mean.z), hg_to_ggx(asymmetry_T_0j_R), ior_ji)) * ops[i].transmission_up.norm;
                
                ops[i].reflection_up.norm += tir_norm;
                ops[i].transmission_up.norm += tir_norm;

            }
            
            energy_0i = mul(energy_0i, energy_matrix(ops[i]));
            asymmetry_0i = mul(asymmetry_0i, asymmetry_matrix(ops[i]));
            
            energy_r_0i = reflection_energy_tm2(energy_0i);
            asymmetry_r_0i = asymmetry_0i._21 / asymmetry_0i._11;

        }
        else
        {
            energy_0i = reflection_energy_tm2(energy_0i, ops[i].reflection_down.norm);
            asymmetry_0i = reflection_energy_tm2(asymmetry_0i, float3_average(ops[i].reflection_down.norm) * ops[i].reflection_down.asymmetry));

        }
        
        energy_r_i = energy_r_0i - energy_r_0h;
        energy_r_i_average = float3_average(energy_r_i);
        
        lobes[i].norm = energy_r_i;
        lobes[i].asymmetry = energy_r_i_average > 0.0f ? min((asymmetry_r_0i - asymmetry_r_0h) / energy_r_i_average, 1.0f) : 0.0f;

        
        energy_r_0h = energy_r_0h;
        asymmetry_r_0h = asymmetry_r_0i;
        
        asymmetry_T_0i = asymmetry_T_0j;
    }

}

float3 eval(sample_record rec, int measure, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughnesses[LAYERS_MAX])
{
    if (measure != TM2_MEASURE_SOLID_ANGLE || rec.incident.z <= 0 || rec.outgoing.z <= 0)
    {
        return EVAL_DEBUG;
    }
    
    const float3 H = normalize(rec.incident + rec.outgoing);
    
    henyey_greenstein lobes[LAYERS_MAX];
    
    outgoing_lobes(rec.incident, iors, kappas, roughnesses, lobes);
    
    //TODO: bati top reflection correction?
    const float G2_0 = 0.0f;
    const float D0_0 = 0.0f;
    
    float3 F0;
    const float3 ior_01 = iors[1] / iors[0];
    if (kappas[1].x + kappas[1].y + kappas[1].z == 0.0f)
    {
        F0 = 0.0f; //TODO: fresnelDielectricExt
    }
    else
    {
        F0 = 0.0f; //TODO fresnelConductorExact
    }
    
    float3 throughput = F0 * G2 * D0 / (4.0f * rec.incident.z);
    
    for (int i = 1; i < NUM_LAYERS; i++)
    {
        if (lobes[i].norm.x + lobes[i].norm.y + lobes[i].norm.z == 0.0f)
        {
            continue;
        }
        const float rough = hg_to_ggx(lobes[i].asymmetry);
        const float G2 = smithG1(rec.incident, H, rough) * smithG1(rec.outgoing, H, rough);
        const float D = D_GGX(H, rough); //TODO: check that H is correct use...
        
        const float f = G2 * D / (4.0f * rec.incident.z);
        
        throughput += lobes[i].norm * f; 

    }
    
    return throughput;

}


float3 sample(inout sample_record rec, out float pdf, float2 samplePoint, float Hammersley, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughness[LAYERS_MAX])
{
    henyey_greenstein lobes[LAYERS_MAX];
    
    //TODO outgoing lobes
    outgoing_lobes(rec.incident, iors, kappas, alphas, lobes);
    
    //Lobe selection
    
    float w[LAYERS_MAX];
    float w_sum = 0.f;
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        w[i] = float3_average(lobes[i].norm);
        w_sum += w[i];
    }
    
    float sel_w = Hammersley * w_sum - w[0];
    int sel_i = 0;
    for (sel_i = 0; sel_w > 0.f && sel_i < NUM_LAYERS; sel_i++)
    {
        sel_w -= w[sel_i + 1];
    }
    
    //sampling
    
    const float sel_a = hg_to_ggx(lobes[sel_i].asymmetry);
    
    const float3 H = sample_GGX(samplePoint, sel_a, pdf);
    
    rec.outgoing = reflect(rec.incident, H);
    rec.ior = 1.f;
    rec.is_reflection_sample = 0;
    rec.sample_type = TM2_SAMPLE_TYPE_GLOSSY_REFLECTION;
    
    if (rec.outgoing.z <= 0.f || pdf <= 0.0f)
    {
        return PDF_DEBUG;
    }
    
    //PDF
    
    pdf = 0.0f;
    for (int i = 0; i < NUM_LAYERS; i++)
    {
        if (w[i] > 0.0f)
        {
            const float rough = hg_to_ggx(lobes[i].asymmetry);
            
                    
            float3 incoming = rec.incident;
            float3 outgoing = rec.outgoing;
        
            incoming.z = abs(incoming.z);
            outgoing.z = abs(outgoing.z);
        
            const float3 H = normalize(incoming + outgoing);
               
            const float G1 = smithG1(incoming, H, rough);
        
            const float D = D_GGX(H, rough);
        
            pdf += w[i] * G1 * D / (4.0f * incoming.z);
        

        }
    }
    pdf /= w_sum;
    
    //Throughput
    
    float3 throughput = eval(rec, TM2_MEASURE_SOLID_ANGLE, iors, kappas, roughness);
    
    return pdf > 0.f ? throughput / pdf : PDF_DEBUG;

}

float4 main(VSOutput vsOutput) : SV_Target0
{
	    //interface between air and arbitrarysurface.
    float3 iors[LAYERS_MAX] = { float3(1.3f, 0.96521, 0.6177f), float3(1.0f, 1.0f, 1.0f), 1.0f.xxx, 1.0f.xxx, 1.0f.xxx };
    float alphas[LAYERS_MAX] = { 0.9f, 1.0f, 1.0f, 1.0f, 0.2f };
	
    
    float3 normal = normalize(vsOutput.normal);
    
    float3 tangent = normalize(vsOutput.tangent.xyz);
    float3 bitangent = normalize(cross(normal, tangent)) * vsOutput.tangent.w;
    float3x3 WorldToTangent = float3x3(tangent, bitangent, normal);
    float3x3 TangentToWorld = transpose(WorldToTangent);
    
    float3 SunDirectionTangent = mul(WorldToTangent, normalize(float3(1.0f, -0.5f, 0.0f)));
    float3 SunDirectionReflect = reflect(SunDirectionTangent, float3(0.0f, 1.0f, 0.0f));
    
    float3 ViewerRay = normalize(ViewerPos - vsOutput.worldPos);
    
    
    float3 IBLIncidentRay = -ViewerRay;
    
    float3 sunPower = float3(1.0f, 0.5f, 0.1f);
    
    float3 IBLradiance = radianceIBLTexutre.SampleLevel(cubeMapSampler, reflect(-ViewerRay, normal), 0);
    
    
    //do 4 samples?
    sample_record rec =
    {
        float3(0.0f, 0.0f, -1.0f),
        mul(WorldToTangent, -ViewerRay),
        1.333f,
        true,
        TM2_SAMPLE_TYPE_GLOSSY_REFLECTION
    };
    
    float3 accum = 0.0f.xxx;
    
    for (int i = 0; i < NUM_SAMPLES; i++)
    {
        float2 samplePoint = float2(Hammersley(vsOutput.uv0.x, NUM_SAMPLES).x, Hammersley(vsOutput.uv0.y, NUM_SAMPLES).x);
        float3 sampleEnergy = sample(rec, samplePoint, samplePoint.x, iors, alphas);
        float3 outgoing_facing = dot(rec.outgoing, mul(WorldToTangent, -ViewerRay));
        
        accum += IBLradiance * outgoing_facing * sampleEnergy;

        //accum += IBLradiance * eval(rec, TM2_MEASURE_SOLID_ANGLE, iors, alphas);
        //rec.is_reflection_sample = i % 2;
        //rec.sample_type = i % 2;
        
        //rec.sample_type = i % 2; //alternate between glossy transmission and reflection
    }
    if (isnan(accum.x) || isnan(accum.y) || isnan(accum.z) || isinf(accum.x) || isinf(accum.y) || isinf(accum.z))
    {
        accum = NAN_DEBUG;
    }
    
    
	
    return float4(accum, 1.0f);
}
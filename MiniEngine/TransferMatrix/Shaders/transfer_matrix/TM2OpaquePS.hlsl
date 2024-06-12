#include "tm2_common.hlsli"
#include "MatrixOperators.hlsli"


void components_transfer_factors(float3 incident, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughness[LAYERS_MAX], out layer_components_tm2 ops[LAYERS_MAX])
{
    float3 ior_ij = 0.0f;
    
    for (int i = 0; i < NumLayers; i++)
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

void outgoing_lobes(float3 incident, float3 ior[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughness[LAYERS_MAX], out henyey_greenstein lobes[LAYERS_MAX])
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
    
    float3 tir_norm = 0.0f.xxx;
    
    tensor3d2x2 energy_0i =
    {
                               1.0f.xxx, 0.0f.xxx,
                               0.0f.xxx, 1.0f.xxx 
    };
    
    float2x2 asymmetry_0i = float2x2(
                               1, 0,
                               0, 1);
    
    
    components_transfer_factors(incident, ior, kappas, roughness, ops);
    
    for (int i = 0; i < NumLayers; i++)
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
                tir_norm = TIR_lookup(float3(abs(ops[i].reflection_down.mean.z), hg_to_ggx(asymmetry_T_0i), ior_ij)) * ops[i].transmission_down.norm;
            
                ops[i].reflection_down.norm += tir_norm;
                ops[i].transmission_down.norm -= tir_norm;
            }
            else
            {
                tir_norm = TIR_lookup(float3(abs(ops[i].transmission_down.mean.z), hg_to_ggx(asymmetry_T_0j_R), ior_ji)) * ops[i].transmission_up.norm;
                
                ops[i].reflection_up.norm += tir_norm;
                ops[i].transmission_up.norm -= tir_norm; //could go negative if transfer_factors produces negative norm?

            }
            
            energy_0i = mul(energy_0i, energy_matrix(ops[i]));
            asymmetry_0i = mul(asymmetry_0i, asymmetry_matrix(ops[i]));
            
            energy_r_0i = reflection_energy_tm2(energy_0i);
            asymmetry_r_0i = reflection_energy_tm2(asymmetry_0i);

        }
        else
        {
            energy_r_0i = reflection_energy_tm2(energy_0i, ops[i].reflection_down.norm);
            asymmetry_r_0i = reflection_energy_tm2(asymmetry_0i, float3_average(ops[i].reflection_down.norm) * ops[i].reflection_down.asymmetry);

        }
        
        energy_r_i = energy_r_0i - energy_r_0h;
        energy_r_i_average = float3_average(energy_r_i);
        
        lobes[i].norm = energy_r_i;
        lobes[i].asymmetry = energy_r_i_average > 0.0f ? min((asymmetry_r_0i - asymmetry_r_0h) / energy_r_i_average, 1.0f) : 0.0f;

        
        energy_r_0h = energy_r_0i;
        asymmetry_r_0h = asymmetry_r_0i;
        
        asymmetry_T_0i = asymmetry_T_0j;
    }

}

float3 eval(sample_record rec, int measure, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughnesses[LAYERS_MAX])
{
    if (rec.incident.z <= 0)
    {
        return
        float3(0.0f, 0.0f, 0.0f); 
    }
    
    if (measure != TM2_MEASURE_SOLID_ANGLE || rec.outgoing.z <= 0)
    {
        return EVAL_DEBUG;
    }
    
    
    
    const float3 H = normalize(rec.incident + rec.outgoing);
       
    henyey_greenstein lobes[LAYERS_MAX];
    
    outgoing_lobes(rec.incident, iors, kappas, roughnesses, lobes);
    
    //TODO: bati top reflection correction?
    const float G2_0 = smithG(rec.incident, rec.outgoing, H, roughnesses[1]);
    const float D0_0 = D_GGX(H, roughnesses[1]);
    
    float3 F0;
    const float3 ior_01 = iors[1] / iors[0];
    if (kappas[1].x + kappas[1].y + kappas[1].z == 0.0f)
    {
        F0 = fresnelDielectric(dot(rec.incident, H), float3_average(ior_01));

    }
    else
    {
        //F0 = Fresnel_Shlick(f0(ior_01), 1.0f.xxx, dot(rec.incident, H));
        F0 = fresnelConductorExact(dot(rec.incident, H), ior_01, kappas[1] / iors[0]);
    }
    
    float3 throughput = F0 * G2_0 * D0_0 / (4.0f * rec.incident.z);
    
    for (int i = 1; i < NumLayers; i++)
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



float3 sample(inout sample_record rec, out float pdf, out float lobe_rough, float2 samplePoint, float Hammersley, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float roughness[LAYERS_MAX])
{
    if (rec.incident.z < 0)
    {
        return PDF_DEBUG;
    }
    
    henyey_greenstein lobes[LAYERS_MAX];
    
    //TODO outgoing lobes
    outgoing_lobes(rec.incident, iors, kappas, roughness, lobes);
    
    //Lobe selection
    
    float weights[LAYERS_MAX];
    float weight_sum = 0.f;
    for (int i = 0; i < NumLayers; i++)
    {
        weights[i] = float3_average(lobes[i].norm);
        weight_sum += weights[i];
    }
    
    //scale and bias?
    float selected_weight = Hammersley * weight_sum - weights[0];
    int selection_index = 0;
    for (selection_index = 0; selected_weight > 0.f && selection_index < NumLayers; selection_index++)
    {
        selected_weight -= weights[selection_index + 1];
    }
    
    
    //sampling - take random lobe for outgoing direction.
    //should outgoing roughness be an average of all the sampled lobes?
    //or should it just be the top lobe?? 
    
    //lobe_rough = hg_to_ggx(lobes[selection_index].asymmetry);
    float rough = hg_to_ggx(lobes[selection_index].asymmetry);
    
    const float3 H = sample_GGX(samplePoint, rough, pdf);

    rec.outgoing = reflectSpherical(rec.incident, H);
    rec.ior = 1.f;
    rec.is_reflection_sample = false;
    rec.sample_type = TM2_SAMPLE_TYPE_GLOSSY_REFLECTION;
    
    if (rec.outgoing.z <= 0.f || pdf <= 0.0f)
    {
        return PDF_DEBUG;
    }
    
    //PDF
    lobe_rough = 0.0f;
    pdf = 0.0f;
    for (int i = 0; i < NumLayers; i++)
    {
        if (weights[i] > 0.0f)
        {
            const float rough = hg_to_ggx(lobes[i].asymmetry);
            //compute outgoing lobes roughness as weighted average for IBL evaluation.
            lobe_rough += rough * weights[i];
            
            
                    
            float3 incoming = rec.incident;
            float3 outgoing = rec.outgoing;
        
            incoming.z = abs(incoming.z);
            outgoing.z = abs(outgoing.z);
        
            const float3 H = normalize(incoming + outgoing);
               
            const float G1 = smithG1(incoming, H, rough);
        
            const float D = D_GGX(H, rough);
        
            pdf += weights[i] * G1 * D / (4.0f * incoming.z);
        

        }
    }
    
    pdf /= weight_sum;
    
    //Throughput
    
    float3 throughput = eval(rec, TM2_MEASURE_SOLID_ANGLE, iors, kappas, roughness);
    
    return pdf > 0.f ? throughput / pdf : PDF_DEBUG;

}

float compute_pdf(sample_record rec, int measure, float3 iors[LAYERS_MAX], float3 kappas[LAYERS_MAX], float alphas[LAYERS_MAX])
{
    if (measure != TM2_MEASURE_SOLID_ANGLE || rec.incident.z <= 0 || rec.outgoing.z <= 0)
    {
        return 0.0f;
    }
    
    const float3 H = normalize(rec.incident + rec.outgoing);
    
    henyey_greenstein lobes[LAYERS_MAX];
    outgoing_lobes(rec.incident, iors, kappas, alphas, lobes);
    
    float w_sum = 0.0f;
    float wpdf_sum = 0.0f;
    
    for (int i = 0; i < NumLayers; i++)
    {
        if (lobes[i].norm.x + lobes[i].norm.y + lobes[i].norm.z == 0)
        {
            continue;
        }
        
        float rough = hg_to_ggx(lobes[i].asymmetry);
        
        float G1 = smithG1(rec.incident, H, rough);
        float D = D_GGX(H, rough);
        
        float norm = float3_average(lobes[i].norm);
        
        w_sum += norm;
        wpdf_sum += norm * G1 * D / (4.0f * rec.incident.z);

    }
    
    return w_sum > 0.0f ? wpdf_sum / w_sum : 0.0f;
    
    
}


float4 main(VSOutput vsOutput) : SV_Target0
{
    float3 normal = normalize(vsOutput.normal);
    
    float3 tangent = normalize(vsOutput.tangent.xyz);
    float3 bitangent = normalize(cross(normal, tangent)) * vsOutput.tangent.w;
    float3x3 WorldToTangent = float3x3(tangent, bitangent, normal);
    float3x3 TangentToWorld = transpose(WorldToTangent);
    
    float3 ViewerRay = normalize(ViewerPos - vsOutput.worldPos);
    
    sample_record rec =
    {
        //feel like there's something wrong with the CRS i'm providing,
        //but I don't know what.
        cartesianTSToMitsubaLS(mul(WorldToTangent, ViewerRay)),
        normalize(cartesianTSToMitsubaLS(mul(WorldToTangent, reflect(-ViewerRay, normal)))),

        1.0f,
        true,
        TM2_SAMPLE_TYPE_GLOSSY_REFLECTION
    };
    
    float3 accumulated_energy = 0.0f.xxx;
    float accumulated_rough = 0.0f;
    float3 accumulated_sample_direction = 0.0f.xxx;
    
   // for (int i = 0; i < NumSamples; i++)
   {
        float2 samplePoint = 0.f.xx; //Hammersley(i, NumSamples);
        float pdf;
        float3 sampleEnergy = sample(rec, pdf, accumulated_rough, samplePoint, Hammersley(0, 5).y, IORs, Kappas, Roughs);
      
        float3 outgoingWS = mul(TangentToWorld, MitsubaLSToCartesianTS(rec.outgoing));
        
        accumulated_energy += sampleEnergy;
        accumulated_sample_direction = outgoingWS;
        
       
    }
    
    float rough = accumulated_rough;
    float3 energy = accumulated_energy;
    float3 direction = accumulated_sample_direction;
    
    float lod = rough * IBLRange + IBLBias;
    float3 IBLSample = radianceIBLTexutre.SampleLevel(cubeMapSampler, direction, lod);
    
    
    float3 output = accumulated_energy * IBLSample;
    
	
    return float4(output, 1.0f);
}
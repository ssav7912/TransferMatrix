#include "tm2_common.hlsli"
#include "MatrixOperators.hlsli"


void components_transfer_factors(float3 incident, real3 iors[LAYERS_MAX], real3 kappas[LAYERS_MAX], real roughness[LAYERS_MAX], out layer_components_tm2 ops[LAYERS_MAX])
{
    real3 ior_ij = 0.0;
    
    [loop]
    for (int i = 0; i < NumLayers; i++)
    {
        ops[i] = zero_init_tm2_components();
        ior_ij = iors[i + 1] / iors[i];
        if (kappas[i+1].x + kappas[i+1].y + kappas[i+1].z == 0.0)
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

void outgoing_lobes(float3 incident, real3 ior[LAYERS_MAX], real3 kappas[LAYERS_MAX], real roughness[LAYERS_MAX], out hg_nomean lobes[LAYERS_MAX])
{
    layer_components_tm2 ops[LAYERS_MAX]; 
    
    real ior_ij = 0.0;
    real ior_ji = 0.0;
    
    real3 energy_r_i = 0.0;
    real energy_r_i_average = 0.0;
    
    real3 energy_r_0h = 0.0.xxx;
    real3 energy_r_0i = 0.0.xxx;
    
    real asymmetry_r_0h = 0.0;
    real asymmetry_r_0i = 0.0;
    
    real asymmetry_T_0i = 1.0;
    real asymmetry_T_0j = 1.0;
    real asymmetry_T_0j_R = 1.0;
    real asymmetry_T_0j_RT = 1.0;
    
    real3 tir_norm = 0.0.xxx;
    
    tensor3d2x2 energy_0i =
    {
                               1.0.xxx, 0.0.xxx,
                               0.0.xxx, 1.0.xxx 
    };
    
    real2x2 asymmetry_0i = real2x2(
                               1.0, 0.0,
                               0.0, 1.0);
    
    
    components_transfer_factors(incident, ior, kappas, roughness, ops);
    
    [loop]
    for (int i = 0; i < NumLayers; i++)
    {
        
        
        //goddamn numerical robustness
        //wonder what the cost of throwing safe_divs everywhere is...
        lobes[i] = zero_hg_nomean();
        if (ops[i].component_type == TM_TYPE_DIELECTRICINTERFACE)
        {
            ior_ij = float3_average((ior[i + 1] / ior[i]));
            ior_ji = 1.0 / ior_ij;
            
            asymmetry_T_0j = hg_refract(asymmetry_T_0i, ior_ji) * ops[i].transmission_down.asymmetry;
            asymmetry_T_0j_R = asymmetry_T_0j * ops[i + 1].reflection_down.asymmetry; //uninitialised??
            asymmetry_T_0j_RT = hg_refract(asymmetry_T_0j_R, ior_ij) * ops[i].transmission_up.asymmetry;
            
            ops[i].transmission_down.asymmetry = safe_div(asymmetry_T_0j, asymmetry_T_0i);
            
            ops[i].transmission_up.asymmetry = safe_div(asymmetry_T_0j_RT, asymmetry_T_0j_R);
        
            //toggle to disable the TIR correction. Breaks conservation but seems to massively improve performance.
#if !defined(DISABLE_TIR) || DISABLE_TIR == 0
            if (ior_ij < 1.0)
            {
#if !defined ANALYTIC_TIR || ANALYTIC_TIR == 0
                tir_norm = TIR_lookup(float3(abs(ops[i].reflection_down.mean.z), hg_to_ggx(asymmetry_T_0i), ior_ij)) * ops[i].transmission_down.norm;
      
#else
                tir_norm = TIR_analytical(abs(ops[i].reflection_down.mean.z), hg_to_ggx(asymmetry_T_0i), ior_ij, 1.0/float3_average(ior[i])) * ops[i].transmission_down.norm;
#endif
                ops[i].reflection_down.norm += tir_norm;
                ops[i].transmission_down.norm -= tir_norm;
            }
            else
            {
#if !defined ANALYTIC_TIR || ANALYTIC_TIR == 0

                tir_norm = TIR_lookup(float3(abs(ops[i].transmission_down.mean.z), hg_to_ggx(asymmetry_T_0j_R), ior_ji)) * ops[i].transmission_up.norm;

#else            
                tir_norm = TIR_analytical(abs(ops[i].transmission_down.mean.z), hg_to_ggx(asymmetry_T_0j_R), ior_ji, 1.0 / float3_average(ior[i + 1])) * ops[i].transmission_up.norm;
#endif

                
                ops[i].reflection_up.norm += tir_norm;
                ops[i].transmission_up.norm -= tir_norm; //could go negative if transfer_factors produces negative norm?
            
            }
#endif
            
            
            energy_0i = mul(energy_0i, energy_matrix(ops[i])); //transmission_down.norm == 0?
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

        
        lobes[i].asymmetry = min(safe_div((asymmetry_r_0i - asymmetry_r_0h), energy_r_i_average), 1.0);
        
        energy_r_0h = energy_r_0i;
        asymmetry_r_0h = asymmetry_r_0i;
        
        asymmetry_T_0i = asymmetry_T_0j;
    }

}

real3 eval_lobe(const sample_record rec, const float3 H, const hg_nomean lobe)
{
    

    if (lobe.norm.x + lobe.norm.y + lobe.norm.z == 0.0)
    {
        return 0.0;
    }
    const real rough = hg_to_ggx(lobe.asymmetry);
#if !defined(USE_EARL_G2) || USE_EARL_G2 == 0
    const real G2 = smithG(rec.incident, rec.outgoing, H, rough);
#else
    const real G2 = smithG2(rec.outgoing, H, rec.incident, rough); 
#endif
    
#if !defined(USE_D_KARIS) || USE_D_KARIS == 0
    const real D = D_GGX(H, rough); //TODO: check that H is correct use...
#else
    const real3 N = float3(0,0,1);
    const real D = D_GGX_Karis(saturate(dot(N, H)), rough);
#endif
    const real3 f = (G2 * D) / (4.0 * rec.incident.z);
        
    float3 throughput = lobe.norm * f;

    
    return throughput;
}

real3 sample_preintegrated(sample_record rec, float3x3 TangentToWorld, LayerProperties props)
{
    if (rec.incident.z <= 0 || rec.outgoing.z <= 0)
    {
        return EVAL_DEBUG;
    }
    
    const float3 N = float3(0., 0., 1.); //H == N when not using specular_dominat. 

        
    //get the outgoing lobes
    hg_nomean lobes[LAYERS_MAX];
    outgoing_lobes(rec.incident, props.iors, props.kappas, props.rough, lobes);
    
    //shift the outgoing lobe direction to correct for lobe mean.
       
    
    real3 throughput = 0.0;
    float3 IBLSamples = 0.0;
    real weight_sum = 0.0;
    
    //compute PDF
        real pdf = 0.0;
    {
        for (int j = 0; j < NumLayers; j++)
        {

            const real lobe_weight = float3_average(lobes[j].norm);
            weight_sum += lobe_weight;
            if (lobe_weight > 0.0)
            {
                
            
                const real rough = hg_to_ggx(lobes[j].asymmetry);
            
                float3 incoming = rec.incident;
                //model lobe shift due to high roughness
                float3 outgoing = specular_dominant(N, rec.outgoing, saturate(dot(N, rec.incident)), rough);
                const float3 H = normalize(incoming + outgoing);
            
#if SCHLICK_G == 1
                const real G1 = SchlickG1(incoming, rough);
#else
                const real G1 = smithG1(incoming, H, rough);
                
#endif
                const real D = D_GGX(H, rough);
            
                pdf += ((lobe_weight) * (G1 * D / (4.0 * incoming.z)));
            }
        }
        pdf = safe_div(pdf, weight_sum);
    }
    
    
    
    [loop]
    for (int i = 0; i < NumLayers; i++)
    {
        if (!isZero(lobes[i].norm))
        {
            
        
            real lobe_weight = float3_average(lobes[i].norm);
            const real rough = hg_to_ggx(lobes[i].asymmetry);
        
            if (i == 0)
            {
                
                //top reflection correction. [Bati 2019]
                const real TopRough = props.rough[1];
                const float3 TopLobeOutgoing = specular_dominant(N, rec.outgoing, saturate(dot(N, rec.incident)), TopRough);
                const float3 H = normalize(rec.incident + TopLobeOutgoing);
#if !defined(USE_EARL_G2) || USE_EARL_G2 == 0
                const real G2_0 = smithG(rec.incident, TopLobeOutgoing, H, TopRough);
#else
                const real G2_0 = smithG2(TopLobeOutgoing, H, rec.incident, TopRough);
#endif
                
                const float3 topOutgoingWS = mul(TangentToWorld, MitsubaLSToCartesianTS(TopLobeOutgoing));
                const float BottomLOD = TopRough * IBLRange + IBLBias;
                const float3 TopIBLSample = radianceIBLTexutre.SampleLevel(cubeMapSampler, topOutgoingWS, BottomLOD);

#if !defined(USE_D_KARIS) || USE_D_KARIS == 0
                const real D0_0 = D_GGX(H, TopRough);
#else
                const real D0_0 = D_GGX_Karis(dot(N, H), rough);
#endif
                real3 F0 = 0.0.xxx;
                const real3 ior_01 = props.iors[1] / props.iors[0];
                if (isZero(props.kappas[1]))
                {
                    F0 = fresnelDielectric(saturate(dot(rec.incident, H)), float3_average(ior_01));
                    //F0 = fresnelDielectric(saturate(dot(rec.incident, H)), ior_01);
                }
                else
                {
                    F0 = fresnelConductorExact(saturate(dot(rec.incident, H)), ior_01, props.kappas[1] / props.iors[0]);
                }

                IBLSamples += (TopIBLSample / (real) NumLayers);
                throughput += ((F0 * G2_0 * D0_0/ (4.0 * rec.incident.z)));
            }
              
             //only evaluate layers when not in air-substrate interface.
            if (i >= 1)
            {
                
                const float3 lobeOutgoing = specular_dominant(N, rec.outgoing, saturate(dot(N, rec.incident)), rough);
                const float3 H = normalize(rec.incident + lobeOutgoing);
#if !defined(USE_EARL_G2) || USE_EARL_G2 == 0
                const real G = smithG(rec.incident, normalize(lobeOutgoing), H, rough);
#else
                const real G = smithG2(lobeOutgoing, H, rec.incident, rough);
#endif
                
                const float3 lobeOutgoingWS = mul(TangentToWorld, MitsubaLSToCartesianTS(lobeOutgoing));
                
                const float LOD = rough * IBLRange + IBLBias;
                const float3 IBLSample = radianceIBLTexutre.SampleLevel(cubeMapSampler, lobeOutgoingWS, LOD);
                //throughput
                
                const sample_record lobe_rec = { rec.incident, normalize(lobeOutgoing)  };
                
                const real3 individual_lobe = eval_lobe(lobe_rec, H, lobes[i]);
                throughput += individual_lobe;
            
                IBLSamples += (IBLSample / (real) NumLayers);
                
            }
 
        }
    }
    
   
    return !all(isnan(throughput)) && !isnan(pdf) ?
    safe_div(throughput, pdf) * IBLSamples : NAN_DEBUG;
}



float4 main(VSOutput vsOutput) : SV_Target0
{
    float3 normal = normalize(vsOutput.normal);
    
    float3 tangent = normalize(vsOutput.tangent.xyz);
    float3 bitangent = normalize(cross(normal, tangent)) * vsOutput.tangent.w;
    float3x3 WorldToTangent = float3x3(tangent, bitangent, normal);
    float3x3 TangentToWorld = transpose(WorldToTangent);
    
    const float3 ViewerRay = normalize(ViewerPos - vsOutput.worldPos);
    
    const float3 incident = cartesianTSToMitsubaLS(mul(WorldToTangent, ViewerRay));
    const float3 outgoing = reflectSpherical(incident, float3(0., 0., 1.0)); //reflect about N
    
    
    
    
    sample_record rec =
    {
        incident,
        outgoing,
    };
    
    
    LayerProperties props = truncate_layer_parameters();
   
    float3 output = sample_preintegrated(rec, TangentToWorld, props);
    
    //layer_components_tm2 ops[LAYERS_MAX];
    //components_transfer_factors(rec.incident, props.iors, props.kappas, props.rough, ops);
    //
    //hg_nomean lobes[LAYERS_MAX];
    //outgoing_lobes(rec.incident, props.iors, props.kappas, props.rough, lobes);
    //
    //output = lobes[NumLayers-1].norm.xxx;
    
    //output = sample_FGD(1.0, props.rough[NumLayers], props.iors[NumLayers], props.kappas[NumLayers]);
    
    
    return float4(output, 1.0);
}
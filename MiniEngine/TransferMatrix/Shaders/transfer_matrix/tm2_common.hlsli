#pragma once

#include "TransferMatrixCommon.hlsli"
#include "HenyeyGreenstein.hlsli"
#include "tensor3d.hlsli"
#include "../Common.hlsli"




struct layer_components_tm2
{
    henyey_greenstein reflection_down;
    henyey_greenstein transmission_down;
    henyey_greenstein reflection_up;
    henyey_greenstein transmission_up;
    int component_type;
    
};

layer_components_tm2 zero_init_tm2_components()
{
    layer_components_tm2 x;
    x.reflection_down = zero_hg();
    x.transmission_down = zero_hg();
    x.reflection_up = zero_hg();
    x.transmission_up = zero_hg();
    x.component_type = TM_TYPE_NOCOMPONENT;
    return x;
}


tensor3d2x2 energy_matrix(layer_components_tm2 ops)
{
    //if the transmission energy is 0, should we make this the identity matrix?
    tensor3d2x2 e = { 0.0.xxx, 0.0.xxx, 0.0.xxx, 0.0.xxx };
    //can i saturate this.. 
    const real3 t_inv = safe_div(1.0.xxx, ops.transmission_down.norm);
    
    e._11 = t_inv;
    e._12 = -ops.reflection_up.norm * t_inv;
    e._21 = ops.reflection_down.norm * t_inv;
    e._22 = (-ops.reflection_down.norm * ops.reflection_up.norm + ops.transmission_down.norm * ops.transmission_up.norm) * t_inv;
    
    return e;
}

real2x2 asymmetry_matrix(layer_components_tm2 ops)
{
    real2x2 asymmetry;
    
    const real r_asymmetry = ops.reflection_down.asymmetry * float3_average(ops.reflection_down.norm);
    const real t_asymmetry = ops.transmission_down.asymmetry * float3_average(ops.transmission_down.norm);
    const real rp_asymmetry = ops.reflection_up.asymmetry * float3_average(ops.reflection_up.norm);
    const real tp_asymmetry = ops.transmission_up.asymmetry * float3_average(ops.transmission_up.norm);
    
    const real t_inv = safe_div(1.0, t_asymmetry);
    
    asymmetry._11 = t_inv;
    asymmetry._12 = -rp_asymmetry * t_inv;
    asymmetry._21 = r_asymmetry * t_inv;
    asymmetry._22 = (-r_asymmetry * rp_asymmetry + t_asymmetry * tp_asymmetry) * t_inv;
    
    return asymmetry; 
}



void dielectric_transfer_factors(float3 incident, real ior, real alpha, out layer_components_tm2 ops)
{
    real ior_ji = 0.0;
    real s_t_ij = 0.0;
    real s_t_ji = 0.0;
    ops = zero_init_tm2_components();

    if (abs(ior - 1.0) < 1e-5)
    {
        ops.component_type = TM_TYPE_NOCOMPONENT;
        ops.transmission_down.mean = -incident;
        
        return;
    }
    
    ops.component_type = TM_TYPE_DIELECTRICINTERFACE;
    
    ior_ji = 1.0 / ior;
    
    ops.reflection_down.asymmetry = ggx_to_hg(alpha);
    ops.reflection_down.mean = reflectZ(incident);
    
    ops.reflection_up.asymmetry = ops.reflection_down.asymmetry;
    ops.reflection_up.mean = refractZ(incident, ior);

    s_t_ij = abs(ior_ji * ops.reflection_down.mean.z + ops.reflection_up.mean.z) / ops.reflection_up.mean.z;
    s_t_ji = abs(ior * ops.reflection_up.mean.z + ops.reflection_down.mean.z) / ops.reflection_down.mean.z;

    ops.transmission_down.asymmetry = ggx_to_hg(0.5 * s_t_ij * alpha);
    ops.transmission_down.mean = ops.reflection_up.mean;
    
    ops.transmission_up.asymmetry = ggx_to_hg(0.5 * s_t_ji * alpha);
    ops.transmission_up.mean = ops.reflection_down.mean;
    
    albedos(abs(incident.z), alpha, ior, ops.reflection_down.norm, ops.transmission_down.norm, ops.reflection_up.norm, ops.transmission_up.norm);

}

void conductor_transfer_factors(float3 incident, real3 ior, real3 kappa, real rough, out layer_components_tm2 ops)
{
    ops = zero_init_tm2_components();
    ops.component_type = TM_TYPE_CONDUCTORINTERFACE;
    
    ops.reflection_down.asymmetry = ggx_to_hg(rough);
    ops.reflection_down.mean = reflectZ(incident);
    
    albedo(abs(incident.z), rough, ior, kappa, ops.reflection_down.norm);
}


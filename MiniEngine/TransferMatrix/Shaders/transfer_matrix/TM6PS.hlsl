#include "TransferMatrixCommon.hlsli"
#include "tensor3d.hlsli"

#define TM6_LOBE_MAX (LAYERS_MAX * 3)
#define TM_EXP_ARG_MAX 45.0f

static const uint TM_TRANSFER_FACTOR_NONE = 0;
static const uint TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE = 0x01;
static const uint TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE = 0x02;
static const uint TM_TRANSFER_FACTOR_SECONDARY_BACKWARD_REFLECTANCE = 0x04;

static const int NumComponents = NumLayers * 2 - 1;

struct LayerProperties
{
    float3 iors[LAYERS_MAX];
    float3 kappas[LAYERS_MAX];
    float rough[LAYERS_MAX];
    float depths[LAYERS_MAX];
    float3 sigma_s[LAYERS_MAX];
    float3 sigma_k[LAYERS_MAX];
    float gs[LAYERS_MAX];
};

struct layer_components_tm6
{
    uint type;
    struct {
        henyey_greenstein reflection_down;
        henyey_greenstein transmission_down;
        henyey_greenstein reflection_up;
        henyey_greenstein transmission_up;
        
        henyey_greenstein reflection_down_secondary;
        henyey_greenstein transmission_down_secondary;
        henyey_greenstein reflection_up_secondary;
        henyey_greenstein transmission_up_secondary;
        
    } interface_type;
    
    struct
    {
        henyey_greenstein primary_flux_transmission;
        henyey_greenstein secondary_flux_backward_reflection;
        henyey_greenstein secondary_flux_forward_transmission;
        
    } media_type;
};


tensor3d6x6 energy_matrix(layer_components_tm6 ops)
{
    tensor3d6x6 energy;
    energy.type = ops.type;
    
    if (ops.type == TM_TYPE_DIELECTRICINTERFACE)
    {
        energy._11 = 1.0f.xxx / ops.interface_type.transmission_down.norm;
        energy._12 = -ops.interface_type.reflection_down.norm * energy._11;
        energy._21 = ops.interface_type.reflection_down.norm * energy._11;
        energy._22 = -ops.interface_type.reflection_down.norm * ops.interface_type.reflection_up.norm * energy._11 + ops.interface_type.transmission_up.norm;
        
        energy._33 = 1.0f.xxx / ops.interface_type.transmission_down_secondary.norm;
        energy._34 = -ops.interface_type.reflection_up_secondary.norm * energy._33;
        energy._43 = ops.interface_type.reflection_down_secondary.norm * energy._33;
        energy._44 = -ops.interface_type.reflection_down_secondary.norm * ops.interface_type.reflection_up.norm * energy._33 + ops.interface_type.transmission_up_secondary.norm;
    }
    else if (ops.type == TM_TYPE_HOMOGENOUSMEDIUM)
    {
        energy._11 = 1.0f.xxx / ops.media_type.primary_flux_transmission.norm;
        energy._22 = ops.media_type.primary_flux_transmission.norm;
        energy._33 = 1.0f.xxx / ops.media_type.secondary_flux_forward_transmission.norm;
        energy._44 = -(ops.media_type.secondary_flux_backward_reflection.norm * ops.media_type.secondary_flux_backward_reflection.norm) * energy._33 + ops.media_type.secondary_flux_forward_transmission.norm;
        energy._52 = energy._36 - ops.media_type.secondary_flux_backward_reflection.norm * energy._33;
        energy._31 = energy._33 - energy._11;
        energy._42 = energy._44 - ops.media_type.primary_flux_transmission.norm;
        
        energy._61 = -energy._36;
        energy._45 = -energy._36;
        
            
    }
    
    return energy;
}

tensor3d6x6 asymmetry_matrix(layer_components_tm6 ops)
{
    tensor3d6x6 asymmetry = tensor3d6x6_identity();
    asymmetry.type = ops.type;
    
    if (ops.type == TM_TYPE_DIELECTRICINTERFACE)
    {
        const float3 reflection_down_asymmetry = ops.interface_type.reflection_down.norm * ops.interface_type.reflection_down.asymmetry;
        const float3 transmission_down_asymmetry = ops.interface_type.transmission_down.norm * ops.interface_type.transmission_down.asymmetry;
        const float3 reflection_up_asymmetry = ops.interface_type.reflection_up.norm * ops.interface_type.reflection_up.asymmetry;
        const float3 transmission_up_asymmetry = ops.interface_type.transmission_up.norm * ops.interface_type.transmission_up.asymmetry;
        
        const float3 reflection_down_secondary_asymmetry = ops.interface_type.reflection_down_secondary.norm * ops.interface_type.reflection_down_secondary.asymmetry;
        const float3 transmission_down_secondary_asymmetry = ops.interface_type.transmission_down_secondary.norm * ops.interface_type.transmission_down_secondary.asymmetry;
        const float3 reflection_up_secondary_asymmetry = ops.interface_type.reflection_up_secondary.norm * ops.interface_type.reflection_up_secondary.asymmetry;
        const float3 transmission_up_secondary_asymmetry = ops.interface_type.transmission_up_secondary.norm * ops.interface_type.transmission_up_secondary.asymmetry;
        
        asymmetry._11 = 1.0f.xxx / transmission_down_asymmetry;
        asymmetry._12 = -reflection_up_asymmetry * asymmetry._11;
        asymmetry._21 = reflection_down_asymmetry * asymmetry._11;
        asymmetry._22 = -reflection_down_asymmetry * reflection_up_asymmetry * asymmetry._11 + transmission_up_asymmetry;
        
        asymmetry._33 = 1.0f.xxx / transmission_up_secondary_asymmetry;
        asymmetry._34 = -reflection_up_secondary_asymmetry * asymmetry._33;
        asymmetry._43 = reflection_down_secondary_asymmetry * asymmetry._33;
        asymmetry._44 = -reflection_down_secondary_asymmetry * reflection_up_secondary_asymmetry * asymmetry._33 + transmission_up_secondary_asymmetry;
        
    }
    else if (ops.type == TM_TYPE_HOMOGENOUSMEDIUM)
    {
        const float3 transmission_asymmetry = ops.media_type.primary_flux_transmission.norm * ops.media_type.primary_flux_transmission.asymmetry;
        const float3 backwards_reflection_asymmetry = ops.media_type.secondary_flux_backward_reflection.norm * ops.media_type.secondary_flux_backward_reflection.asymmetry;
        const float3 transmission_secondary_asymmetry = ops.media_type.secondary_flux_forward_transmission.norm * ops.media_type.secondary_flux_forward_transmission.asymmetry;
        
        asymmetry._11 = 1.0f.xxx / transmission_asymmetry;
        asymmetry._22 = transmission_asymmetry;
        asymmetry._33 = 1.0f.xxx / transmission_secondary_asymmetry;
        asymmetry._44 = -(backwards_reflection_asymmetry * backwards_reflection_asymmetry) * asymmetry._33 + transmission_secondary_asymmetry;
        
        asymmetry._36 = -backwards_reflection_asymmetry * asymmetry._33;
        asymmetry._52 = asymmetry._36;
        
        asymmetry._31 = asymmetry._33 - asymmetry._11;
        asymmetry._42 = asymmetry._44 - transmission_asymmetry;
        
        asymmetry._45 = -asymmetry._36;
        asymmetry._61 = asymmetry._45;
    }
    
    return asymmetry;
}

void matrix_factors(tensor3d6x6 mat, inout float3 r, inout float3 r_f, inout float3 r_b, uint req)
{
    const float3 x0 = 1.0f.xxx / mat._11;
    if (req & TM_TRANSFER_FACTOR_SECONDARY_BACKWARD_REFLECTANCE)
    {
        const float3 x1 = -mat._31 * mat._33 + mat._35 * mat._51;
        const float3 x2 = mat._31 * mat._35 - mat._33 * mat._51;
        const float3 x3 = x0 / (mat._33 * mat._33 - mat._35 * mat._35);
        const float3 x4 = mat._43 * x3;
        const float3 x5 = mat._45 * x3;
        
        r_b = mat._61 * x0 + x1 * x5 + x4 * x2;
        if (req & TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE)
        {
            r_f = mat._41 * x0 + x1 * x4 + x2 * x5;
        }

    }
    if (req & TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE)
    {
        r = mat._21 * x0;
    }
}

void matrix_factors(tensor3d6x6 mat, inout float3 r, inout float3 r_f, inout float3 r_b, float3 r_cond, uint req)
{
    const float3 x0 = 1.0f.xxx / (mat._11 + mat._12 * r_cond);
    if (req & TM_TRANSFER_FACTOR_SECONDARY_BACKWARD_REFLECTANCE)
    {
        const float3 x1 = mat._45 + mat._46 * r_cond;
        const float3 x2 = mat._31 + mat._32 * r_cond;
        const float3 x3 = mat._35 + mat._36 * r_cond;
        const float3 x4 = mat._33 + mat._34 * r_cond;
        const float3 x5 = mat._51 + mat._52 * r_cond;
        const float3 x6 = x0 / (-x3 * x3 + x4 * x4);
        const float3 x7 = x6 * (x2 * x3 - x4 * x5);
        const float3 x8 = mat._43 + mat._44 * r_cond;
        const float3 x9 = x6 * (x2 * x4 - x3 * x5);
        
        r_b = x0 * (mat._61 + mat._62 * r_cond) - x1 * x9 + x7 * x8;
        if (req & TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE)
        {
            r_f = x0 * (mat._41 + mat._42 * r_cond) + x1 * x7 - x8 * x9;
        }
    }
    
    if (req & TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE)
    {
        r = x0 * (mat._21 + mat._22 * r_cond);
    }
}


void dielectric_transfer_factors(float3 incident, float ior_ij, float rough, out layer_components_tm6 ops)
{
    float ior_ji = 0.0f;
    float s_t_ij = 0.0f;
    float s_t_ji = 0.0f;
   
    if (abs(ior_ij - 1.0f) < 1e-5f)
    {
        ops.type = TM_TYPE_NOCOMPONENT;
        ops.interface_type.transmission_down.mean = -incident;
        return;
    }
    
    ops.type = TM_TYPE_DIELECTRICINTERFACE;
    ior_ji = 1.0f / ior_ij;
    
    ops.interface_type.reflection_down.asymmetry = ggx_to_hg(rough);
    ops.interface_type.reflection_down.mean = reflectZ(incident);
    
    ops.interface_type.reflection_up.asymmetry = ops.interface_type.reflection_down.asymmetry;
    ops.interface_type.reflection_up.mean = refractZ(incident, ior_ij);
    
    s_t_ij = abs((ior_ji * ops.interface_type.reflection_down.mean.z + ops.interface_type.reflection_up.mean.z) / ops.interface_type.reflection_up.mean.z);
    s_t_ji = abs((ior_ij * ops.interface_type.reflection_up.mean.z + ops.interface_type.reflection_down.mean.z) / ops.interface_type.reflection_down.mean.z);
    
    ops.interface_type.transmission_down.asymmetry = ggx_to_hg(0.5f * s_t_ij * rough);
    ops.interface_type.transmission_down.mean = ops.interface_type.reflection_up.mean;
    
    ops.interface_type.transmission_up.asymmetry = ggx_to_hg(0.5f * s_t_ji * rough);
    ops.interface_type.transmission_up.mean = ops.interface_type.reflection_down.mean;
    
    albedos(abs(incident.z), rough, ior_ij, ops.interface_type.reflection_down.norm, ops.interface_type.transmission_down.norm, ops.interface_type.reflection_up.norm, ops.interface_type.transmission_up.norm);
    
    ops.interface_type.reflection_down_secondary = ops.interface_type.reflection_down;
    ops.interface_type.transmission_down_secondary = ops.interface_type.transmission_down;
    ops.interface_type.reflection_up_secondary = ops.interface_type.reflection_up;
    ops.interface_type.transmission_down_secondary = ops.interface_type.transmission_down;
}

void conductor_transfer_factors(float3 incident, float3 ior, float3 kappa, float rough, out layer_components_tm6 ops)
{
    ops.type = TM_TYPE_CONDUCTORINTERFACE;
    ops.interface_type.reflection_down.asymmetry = ggx_to_hg(rough);
    ops.interface_type.reflection_down.mean = reflectZ(incident);
    
    albedo(abs(incident.z), rough, ior, kappa, ops.interface_type.reflection_down.norm);
    
    ops.interface_type.reflection_down_secondary = ops.interface_type.reflection_down;
}

void medium_transfer_factors(float3 incident, float depth, float3 sigma_s, float3 sigma_k, float g, layer_components_tm6 ops)
{
    ops.type = TM_TYPE_HOMOGENOUSMEDIUM;
    if (depth == 0.0f)
    {
        ops.type = TM_TYPE_NOCOMPONENT;
        ops.interface_type.transmission_down.mean = -incident;
    }
    
    float tau = depth / incident.z;
    
    const float3 sigma_sb = sigma_s * hg_lh_norm(g);
    const float3 sigma_sf = sigma_s - sigma_sb;
    const float3 sigma_ext = sigma_s + sigma_k;
    {
        const float3 alpha = sigma_ext - sigma_sf;
        const float3 beta = sigma_sb;
        const float3 gamma = sqrt((alpha * alpha - beta * beta));
        
        tau = min(tau, TM_EXP_ARG_MAX / float3_max(gamma));
        tau = min(tau, TM_EXP_ARG_MAX / float3_max(sigma_ext));
        
        const float3 S = sinh(gamma * tau);
        const float3 C = sqrt(1.0f.xxx + S * S);
        
        const float3 x0 = exp(-sigma_ext * tau);
        const float3 x1 = 1.0f.xxx / (C * gamma + S * alpha);
        
        ops.media_type.primary_flux_transmission.norm = x0;
        ops.media_type.secondary_flux_backward_reflection.norm = S * beta * x1;
        ops.media_type.secondary_flux_forward_transmission.norm = gamma * x1; 

    }
    
    {
        const float3 alpha = sigma_ext - sigma_sf * g;
        const float3 beta = sigma_sb * -g;
        const float3 gamma = sqrt(alpha * alpha - beta * beta);
        
        const float3 S = sinh(gamma * tau);
        const float3 C = sqrt(1.0f.xxx + S * S);
        
        const float3 x0 = 1.0f.xxx / (C * gamma + S * alpha);
        
        ops.media_type.primary_flux_transmission.asymmetry = 1.0f;
        //TODO safe_div
        ops.media_type.secondary_flux_backward_reflection.asymmetry = float3_average((S * beta * x0) / (ops.media_type.secondary_flux_backward_reflection.norm));

        ops.media_type.secondary_flux_forward_transmission.asymmetry = float3_average((gamma * x0) / (ops.media_type.secondary_flux_forward_transmission.norm));
        
    }
}

void components_transfer_factors_tm6(float3 incident, LayerProperties props, out layer_components_tm6 ops[TM6_LOBE_MAX])
{
    float3 ior_ij;
    
    for (int i = 0; i < NumLayers; i++)
    {
        ior_ij = props.iors[i + 1] / props.iors[i];
        
        if (isZero(props.kappas[i + 1]))
        {
            dielectric_transfer_factors(incident, float3_average(ior_ij), props.rough[i + 1], ops[i * 2]);
            incident = -ops[i * 2].interface_type.transmission_down.mean;
            
            medium_transfer_factors(incident, props.depths[i + 1], props.sigma_s[i + 1], props.sigma_k[i + 1], props.gs[i + 1], ops[i * 2 + 1]);
        }
        else
        {
            conductor_transfer_factors(incident, ior_ij, props.kappas[i + 1] / props.iors[i], props.rough[i + 1], ops[i * 2]);
        }

    }


}


int outgoing_lobes(float3 incident, LayerProperties props, out henyey_greenstein lobes[TM6_LOBE_MAX], out float3 lobe_incident[TM6_LOBE_MAX])
{
    const float3 incident_reflect = reflectZ(incident);
    
    int count = 0;
    
    layer_components_tm6 ops[TM6_LOBE_MAX];
    
    float ior_ij = 0.0f;
    float ior_ji = 0.0f;
    float3 energy_reflect_i = 0.0f.xxx;
    
    float3 energy_reflect_0h = 0.0f.xxx;
    float3 energy_reflect_0i = 0.0f.xxx;
    
    float3 energy_reflect_f_0h = 0.0f.xxx;
    float3 energy_reflect_f_0i = 0.0f.xxx;
    
    float3 energy_reflect_b_0h = 0.0f.xxx;
    float3 energy_reflect_b_0i = 0.0f.xxx;
    
    
    float3 asymmetry_reflect_0h = 0.0f.xxx;
    float3 asymmetry_reflect_0i = 0.0f.xxx;
    float3 asymmetry_reflect_f_0h = 0.0f.xxx;
    float3 asymmetry_reflect_f_0i = 0.0f.xxx;
    float3 asymmetry_reflect_b_0h = 0.0f.xxx;
    float3 asymmetry_reflect_b_0i = 0.0f.xxx;
    
    float asymmetry_T_0i = 1.0f;
    float asymmetry_T_0j = 1.0f;
    float asymmetry_T_0j_R = 1.0f;
    float asymmetry_T_0j_RT = 1.0f;
    
    float asymmetry_s_T_0i = 1.0f;
    float asymmetry_s_T_0j = 1.0f;
    float asymmetry_s_T_0j_R = 1.0f;
    
    float asymmetry_s_T_0j_RT = 1.0f; 
    
    
    float3 TIR_norm;
    
    uint req = TM_TYPE_NOCOMPONENT; //what is this?
    
    tensor3d6x6 energy_matrix_0i = tensor3d6x6_identity();
    tensor3d6x6 asymmetry_matrix_0i = tensor3d6x6_identity();
    
    
    components_transfer_factors_tm6(incident, props, ops);
    
    
    for (int c = 0; c < NumComponents; c++)
    {
        if (ops[c].type == TM_TYPE_NOCOMPONENT)
        {
            lobes[count].norm = 0.0f.xxx;
            lobes[count].asymmetry = 0.0f;
            lobe_incident[count++] = incident;
            continue;
        }
        else if (ops[c].type == TM_TYPE_DIELECTRICINTERFACE)
        {
            const int i = count >> 1; //hmmm
            ior_ij = float3_average((props.iors[i + 1] / props.iors[i]));
            ior_ji = 1.0f / ior_ij;
            
            asymmetry_T_0j = hg_refract(asymmetry_T_0i, ior_ji) * ops[c].interface_type.transmission_down.asymmetry;
            asymmetry_T_0j_R = asymmetry_T_0j + ops[c+ 2].interface_type.reflection_down.asymmetry;
            asymmetry_T_0j_RT = hg_refract(asymmetry_T_0j_R, ior_ij) * ops[c].interface_type.transmission_up.asymmetry;
            
            asymmetry_s_T_0j = hg_refract(asymmetry_s_T_0i, ior_ji) * ops[c].interface_type.transmission_down_secondary.asymmetry;
            asymmetry_s_T_0j_R = asymmetry_s_T_0j * ops[c + 1].media_type.secondary_flux_forward_transmission.asymmetry * ops[c + 2].interface_type.reflection_down_secondary.asymmetry * ops[c + 1].media_type.secondary_flux_forward_transmission.asymmetry;
            asymmetry_s_T_0j_RT = hg_refract(asymmetry_s_T_0j_R, ior_ij) * ops[c].interface_type.transmission_up_secondary.asymmetry;
            
            //downward transmission;
            ops[c].interface_type.transmission_down.asymmetry = asymmetry_T_0i != 0.f ? asymmetry_T_0j/ asymmetry_T_0i : 0.0f;
            ops[c].interface_type.transmission_down_secondary.asymmetry = asymmetry_s_T_0i != 0.0f ? asymmetry_s_T_0j / asymmetry_s_T_0i : 0.0f;
            
            //upward transmission
            ops[c].interface_type.transmission_up.asymmetry = asymmetry_T_0j_R != 0.0f ? asymmetry_T_0j_RT / asymmetry_T_0j_R : 0.0f;
            ops[c].interface_type.transmission_up_secondary.asymmetry = asymmetry_s_T_0j_R != 0.0f ? asymmetry_s_T_0j_RT / asymmetry_s_T_0j_R : 0.0f;
            
            
            //simplification? drop TIR correction?
            if (ior_ij < 1.0f)
            {
                TIR_norm = TIR_lookup(float3(abs(ops[c].interface_type.reflection_down.mean.z), hg_to_ggx(asymmetry_T_0i), ior_ij)) * ops[c].interface_type.transmission_down.norm;
                
                ops[c].interface_type.reflection_down.norm += TIR_norm;
                ops[c].interface_type.transmission_down.norm -= TIR_norm;

                TIR_norm = TIR_lookup(float3(abs(ops[c].interface_type.reflection_down_secondary.mean.z), hg_to_ggx(asymmetry_s_T_0i), ior_ij)) * ops[c].interface_type.transmission_down_secondary.norm;
                
                ops[c].interface_type.reflection_down_secondary.norm += TIR_norm;
                ops[c].interface_type.transmission_down_secondary.norm += TIR_norm;
            }
            else
            {
                TIR_norm = TIR_lookup(float3(abs(ops[c].interface_type.transmission_down.mean.z), hg_to_ggx(asymmetry_T_0j_R), ior_ji)) * ops[c].interface_type.transmission_up.norm;
                
                ops[c].interface_type.reflection_up.norm += TIR_norm;
                ops[c].interface_type.transmission_up.norm -= TIR_norm;
                
                TIR_norm = TIR_lookup(float3(abs(ops[c].interface_type.transmission_down_secondary.mean.z), hg_to_ggx(asymmetry_s_T_0j_R), ior_ji)) * ops[c].interface_type.transmission_up_secondary.norm;
                
                ops[c].interface_type.reflection_up_secondary.norm += TIR_norm;
                ops[c].interface_type.transmission_up_secondary.norm -= TIR_norm;
                
            }
            
            req |= TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE;
            req |= req & TM_TRANSFER_FACTOR_SECONDARY_BACKWARD_REFLECTANCE ? TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE : TM_TRANSFER_FACTOR_NONE;
            
            energy_matrix_0i = mul(energy_matrix_0i, energy_matrix(ops[c]));
            asymmetry_matrix_0i = mul(asymmetry_matrix_0i, asymmetry_matrix(ops[c]));
            
            matrix_factors(energy_matrix_0i, energy_reflect_0i, energy_reflect_f_0i, energy_reflect_b_0i, req);
            matrix_factors(asymmetry_matrix_0i, asymmetry_reflect_0i, asymmetry_reflect_f_0i, asymmetry_reflect_b_0i, req);
            
            
        }
        else if (ops[c].type == TM_TYPE_HOMOGENOUSMEDIUM)
        {
            req |= TM_TRANSFER_FACTOR_SECONDARY_BACKWARD_REFLECTANCE;
            req |= req & TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE ? TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE : TM_TRANSFER_FACTOR_NONE;
            
            energy_matrix_0i = mul(energy_matrix_0i, energy_matrix(ops[c]));
            asymmetry_matrix_0i = mul(asymmetry_matrix_0i, asymmetry_matrix(ops[c]));
            
            //TODO
            matrix_factors(energy_matrix_0i, energy_reflect_0i, energy_reflect_f_0i, energy_reflect_b_0i, req);
            
            matrix_factors(asymmetry_matrix_0i, asymmetry_reflect_0i, asymmetry_reflect_f_0i, asymmetry_reflect_b_0i, req);
            

        }
        else
        {
            req |= TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE;
            req |= req & TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE ? TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE : TM_TRANSFER_FACTOR_NONE;
            
            //TODO
            matrix_factors(energy_matrix_0i, energy_reflect_0i, energy_reflect_f_0i, energy_reflect_b_0i, ops[c].interface_type.reflection_down.norm, req);
            matrix_factors(asymmetry_matrix_0i, asymmetry_reflect_0i, asymmetry_reflect_f_0i, asymmetry_reflect_b_0i, ops[c].interface_type.reflection_down.norm * ops[c].interface_type.reflection_down.asymmetry, req);
        }
        
        if (req & TM_TRANSFER_FACTOR_PRIMARY_REFLECTANCE)
        {
            energy_reflect_i = energy_reflect_0i - energy_reflect_0h;
            
            lobes[count].norm = max(energy_reflect_i, 0.0f.xxx);
            //TODO: save division
            lobes[count].asymmetry = float3_average((asymmetry_reflect_0i - asymmetry_reflect_0h) / energy_reflect_i);

            lobe_incident[count++] = incident;
        }
        
        if (req & TM_TRANSFER_FACTOR_SECONDARY_FORWARD_REFLECTANCE)
        {
            energy_reflect_i = energy_reflect_f_0i - energy_reflect_f_0h;
            
            lobes[count].norm = max(energy_reflect_i, 0.0f);
            lobes[count].asymmetry = float3_average((asymmetry_reflect_f_0i - asymmetry_reflect_f_0h )/ energy_reflect_i);
            lobe_incident[count++] = incident;
        }
        
        if (req & TM_TRANSFER_FACTOR_SECONDARY_BACKWARD_REFLECTANCE)
        {
            energy_reflect_i = energy_reflect_b_0i - energy_reflect_b_0h;
            
            lobes[count].norm = max(energy_reflect_i, 0.0f.xxx);
            lobes[count].asymmetry = float3_average((asymmetry_reflect_b_0i - asymmetry_reflect_b_0h) / energy_reflect_i);
            lobe_incident[count++] = incident_reflect;

        }
        
        energy_reflect_0h = energy_reflect_0i;
        asymmetry_reflect_0h = asymmetry_reflect_0i;
        
        energy_reflect_f_0h = energy_reflect_f_0i;
        asymmetry_reflect_f_0h = asymmetry_reflect_f_0i;
        
        energy_reflect_b_0h = energy_reflect_b_0i;
        asymmetry_reflect_b_0h = asymmetry_reflect_b_0i;
        
        asymmetry_T_0i = asymmetry_T_0j;
        asymmetry_s_T_0i = asymmetry_s_T_0j;
    }

    return count;

}

float3 eval(sample_record rec, uint measure, LayerProperties props)
{
    if (measure != TM_MEASURE_SOLID_ANGLE || rec.incident.z <= 0 || rec.outgoing.z <= 0)
    {
        return EVAL_DEBUG;
    }
    
    float3 H = normalize(rec.incident + rec.outgoing);
    
    henyey_greenstein lobes[TM6_LOBE_MAX];
    float3 lobes_incident[TM6_LOBE_MAX];
    
    const int lobe_count = outgoing_lobes(rec.incident, props, lobes, lobes_incident);
    
    float3 throughput = 0.f.xxx;
    
    if (!isZero(lobes[0].norm))
    {
        const float G2 = smithG(rec.incident, rec.outgoing, H, props.rough[1]);
        const float D = D_GGX(H, props.rough[1]);
        
        float3 F;
        float3 ior_01 = props.iors[1] / props.iors[0];
        if (isZero(props.kappas[1]))
        {
            F = fresnelDielectric(dot(rec.incident, H), float3_average(ior_01)).xxx;
        }
        else
        {
            F = fresnelConductorExact(dot(rec.incident, H), ior_01, props.kappas[1] / props.iors[0]);
        }
        throughput += F * G2 * D / (4.0f * rec.incident.z);
    }
    
    for (int i = 1; i < lobe_count; i++)
    {
        if (isZero(lobes[i].norm))
        {
            continue;
        }
        
        const float3 incident = lobes_incident[i];
        
        H = normalize(incident + rec.outgoing);
        
        const float rough = hg_to_ggx(lobes[i].asymmetry);
        
        const float G2 = smithG(incident, rec.outgoing, H, rough);
        const float D = D_GGX(H, rough);
        
        //TODO albedo normalisation.
        const float essi = 1.0f;
        
        const float f = G2 * D / (4.0f * incident.z * essi);
        
        throughput += lobes[i].norm * f;
        
    }

    return throughput;
}


float3 sample(sample_record rec, out float pdf, out float out_rough, float2 samplePoint, float Hammersley, LayerProperties props)
{
    henyey_greenstein lobes[TM6_LOBE_MAX];
    float3 lobe_incident[TM6_LOBE_MAX];
    
    const int lobe_count = outgoing_lobes(rec.incident, props, lobes, lobe_incident);
    
    float weight[TM6_LOBE_MAX];
    float weight_sum = 0.0f;
    
    for (int i = 0; i < lobe_count; i++)
    {
        weight[i] = float3_average(lobes[i].norm);
        weight_sum += weight[i];
    }
    
    float selection_weight = Hammersley * weight_sum - weight[0];
    int selection_index = 0;
    for (selection_index = 0; selection_weight > 0.0f && selection_index < lobe_count; selection_index++)
    {
        selection_weight -= weight[selection_index + 1];
    }
    
    const float selection_rough = hg_to_ggx(lobes[selection_index].asymmetry);
    
    float3 selection_incident = lobe_incident[selection_index];
    
    float3 H = sample_GGX(samplePoint, selection_rough, pdf);
    
    rec.outgoing = reflectSpherical(selection_incident, H);
    rec.ior = 1.0f;
    rec.is_reflection_sample = 0;
    rec.sample_type = TM_SAMPLE_TYPE_GLOSSY_REFLECTION;
    
    if (rec.outgoing.z <= 0.0f || pdf <= 0.0f)
    {
        return PDF_DEBUG;
    }
    
    pdf = 0.0f;
    for (int i = 0; i < lobe_count; i++)
    {
        if (weight[i] == 0.0f)
        {
            continue;
        }
        
        float3 incident = lobe_incident[i];
        H = normalize(incident + rec.outgoing);
        
        const float rough = hg_to_ggx(lobes[i].asymmetry);
        out_rough += rough * weight[i];
        
        const float G1 = smithG1(incident, H, rough); 
        const float D = D_GGX(H, rough);
        
        pdf += weight[i] * G1 * D / (4.0f * incident.z);
    }
    
    pdf /= weight_sum;
    
    float3 throughput = eval(rec, TM_MEASURE_SOLID_ANGLE, props);
    
    return pdf > 0.0f ? throughput / pdf : PDF_DEBUG;
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
        TM_SAMPLE_TYPE_GLOSSY_REFLECTION
    };
    
    float3 accumulated_energy = 0.0f.xxx;
    float accumulated_rough = 0.0f;
    float3 accumulated_sample_direction = 0.0f.xxx;
    
   // for (int i = 0; i < NumSamples; i++)
   {
        float2 samplePoint = 0.f.xx; //Hammersley(i, NumSamples);
        float pdf;
        
        LayerProperties props;
        props.iors = IORs;
        props.kappas = Kappas;
        props.depths = Depths;
        props.gs = G;
        props.rough = Roughs;
        props.sigma_k = Sigma_K;
        props.sigma_s = Sigma_S;
        
        
        
        //float3 sampleEnergy = sample(rec, pdf, accumulated_rough, samplePoint, Hammersley(0, 5).y, props);
        float3 sampleEnergy = 0.0f;
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
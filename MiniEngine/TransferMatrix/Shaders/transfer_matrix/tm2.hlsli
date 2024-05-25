#pragma once

#include "HenyeyGreenstein.hlsli"
#include "tensor3d.hlsli"

#define TM2_SAMPLE_TYPE_GLOSSY_REFLECTION 0
#define TM2_SAMPLE_TYPE_GLOSSY_TRANSMISSION 1

#define TM2_MEASURE_SOLID_ANGLE 0


#define TM2_TYPE_NOCOMPONENT 0
#define TM2_TYPE_DIELECTRICINTERFACE 1

struct layer_components_tm2
{
    henyey_greenstein reflection_down;
    henyey_greenstein transmission_down;
    henyey_greenstein reflection_up;
    henyey_greenstein transmission_up;
    int component_type;
    
};

struct sample_record
{
    float3 incident;
    float3 outgoing;
    float ior;
    bool is_reflection_sample;
    bool sample_type;
};

float float3_average(float3 f)
{
    return (f.x + f.y + f.z) / 3.0f;

}

float copy_sign(float x, float s)
{
    return (s >= 0) ? abs(x) : -abs(x);

}


float3 reflectZ(float3 f)
{
    return float3(-f.xy, f.z);
}

float3 refractZ(float3 f, float ior)
{
    return refract(f, float3(0.0f, 0.0f, copy_sign(1.0f, f.z)), ior);

}

float tanTheta(float3 v)
{
    return max(0, sqrt(1 - v.z * v.z) / v.z);
}

float sinTheta(float3 v)
{
    return max(0, sqrt(1.0f - v.z * v.z));
}

//TODO: shortcut
float sinTheta2(float3 v)
{
    return 1.0f - v.z * v.z;
}

float cosTheta2(float3 v)
{
    return v.z * v.z;
}


//float projectRoughness_GGX(float3 v, float3 normal, float alpha)
//{
//    float invSinTheta2 = 1.0f / sinTheta2(v, normal);
//}

float3 sample_GGX(float2 sample, float alpha, out float pdf)
{
    float cosThetaM = 0.0f;
    float sinPhiM = 0.0f;
    float cosPhiM = 0.0f;
    float alphasquare = 0.0f;
    
    //assume isotropic
    sincos(2.0f * PI * sample.y, sinPhiM, cosPhiM);
    alphasquare = alpha * alpha;
    
    float tanthetaMSqr = alphasquare * sample.x / (1.0f - sample.x);
    cosThetaM = 1.0f / sqrt(1.0f + tanthetaMSqr);
    
    float temp = 1.0f + tanthetaMSqr / alphasquare;
    
    //TODO: clamp?
    pdf = (1 / PI) / (alphasquare * cosThetaM * cosThetaM * cosThetaM * temp * temp);

    float sinThetaM = sqrt(max(0.0f, 1 - cosThetaM * cosThetaM));
    
    return float3(sinThetaM * cosPhiM, sinThetaM * sinPhiM, cosThetaM);
}

float D_GGX(float3 m, float alpha)
{
    float costheta2 = cosTheta2(m);
    float beckmann = ((m.x * m.x) / (alpha * alpha) + (m.y * m.y) / (alpha * alpha)) / costheta2;
    
    float root = (1.0f + beckmann) * costheta2;
    return 1.0f / (PI * alpha * alpha * root * root);

}

//assume isotropic
float smithG1(float3 v, float3 m, float alpha)
{
    //v.z == Frame::cosTheta
    if (dot(v, m) * v.z <= 0)
    {
        return 0.0f;
    }
    
    float tantheta = abs(tanTheta(v));
    if (tantheta == 0.0f) //TODO: almost equal op...
    {
        return 1.0f;
    }
    
    
    float root = alpha * tantheta;
                    //hypot2
    return 2.0f / (1.0f + 1.0f * 1.0f + root * root);

}

tensor3d2x2 energy_matrix(layer_components_tm2 ops)
{
    tensor3d2x2 e = { 0.0f.xxx, 0.0f.xxx, 0.0f.xxx, 0.0f.xxx };
    const float3 t_inv = 1.0f / ops.transmission_down.norm;
    
    e._11 = t_inv;
    e._12 = -ops.reflection_up.norm * t_inv;
    e._21 = ops.reflection_down.norm * t_inv;
    e._22 = (-ops.reflection_down.norm * ops.reflection_up.norm + ops.transmission_down.norm * ops.transmission_up.norm) * t_inv;
    
    return e;
}

float2x2 asymmetry_matrix(layer_components_tm2 ops)
{
    float2x2 asymmetry;
    
    const float r_asymmetry = ops.reflection_down.asymmetry * float3_average(ops.reflection_down.norm);
    const float t_asymmetry = ops.transmission_down.asymmetry * float3_average(ops.transmission_down.norm);
    const float rp_asymmetry = ops.reflection_up.asymmetry * float3_average(ops.reflection_up.norm);
    const float tp_asymmetry = ops.transmission_up.asymmetry * float3_average(ops.transmission_up.norm);
    
    const float t_inv = 1.0f / t_asymmetry;
    
    asymmetry._11 = t_inv;
    asymmetry._12 = -rp_asymmetry * t_inv;
    asymmetry._21 = r_asymmetry * t_inv;
    asymmetry._22 = (-r_asymmetry * rp_asymmetry + t_asymmetry * tp_asymmetry) * t_inv;
    
    return asymmetry; 
}

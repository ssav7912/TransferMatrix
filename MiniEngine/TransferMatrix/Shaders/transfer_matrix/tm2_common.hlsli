#pragma once

#include "HenyeyGreenstein.hlsli"
#include "tensor3d.hlsli"
#include "../Common.hlsli"

static const uint TM2_SAMPLE_TYPE_GLOSSY_REFLECTION = 0;
static const uint TM2_SAMPLE_TYPE_GLOSSY_TRANSMISSION = 1;


static const uint TM2_MEASURE_SOLID_ANGLE = 0;


static const uint TM2_TYPE_NOCOMPONENT = 0;
static const uint TM2_TYPE_DIELECTRICINTERFACE = 1;


#define PDF_DEBUG float3(1.0f, 0.0f, 0.0f)
#define EVAL_DEBUG float3(0.0f, 1.0f, 0.0f)
#define NAN_DEBUG float3(1.0f, 0.0f, 1.0f)

#define NUM_SAMPLES 5

#define LAYERS_MAX 5
#define NUM_LAYERS 2
#define NUM_LOBES (NUM_LAYERS + 1)

#define NO_SECOND_UV 1


//Lookup table for Total Internal Reflection
Texture3D<float3> TIR_LUT : register(t18);

//LUT for Karis Split-sum FGD approximation.
Texture2D<float2> FGD_LUT : register(t19);

//IBL
TextureCube<float3> radianceIBLTexutre : register(t10);
TextureCube<float3> irradianceIBLTexture : register(t10);


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
//https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2017/Presentations/Hammon_Earl_PBR_Diffuse_Lighting.pdf
//Earl H. (2017)
float smithG2(float3 V, float3 N, float3 L, float rough)
{
    float numerator = 2.0f * (abs(dot(N, L)) * abs(dot(N, V)));
    float denominator = lerp(2.0f * abs(dot(N, L)) * abs(dot(N, V)), abs(dot(N, L)) + abs(dot(N, V)), rough);
    return numerator / denominator;
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


//PRNG hack for sampling
float nrand(float2 uv)
{
    return frac(sin(dot(uv, float2(12.9898, 78.233))) * 43758.5453);
}

///Hammersley Generator
///https://www.shadertoy.com/view/4lscWj
///Tien-Tsin Wong et al. (1997) Sampling with Hammersley and Halton Points
///Pharr et al. (2021) PBRT.  
float2 Hammersley(float i, float numSamples)
{
    uint b = uint(i);
    b = (b << 16u) | (b >> 16u);
    b = ((b & 0x55555555u) << 1u) | ((b & 0xAAAAAAAAu) >> 1u);
    b = ((b & 0x33333333u) << 2u) | ((b & 0xCCCCCCCCu) >> 2u);
    b = ((b & 0x0F0F0F0Fu) << 4u) | ((b & 0xF0F0F0F0u) >> 4u);
    b = ((b & 0x00FF00FFu) << 8u) | ((b & 0xFF00FF00u) >> 8u);
    
    float radicalInverse = float(b) * 2.3283064365386963e-10;
    
    return float2(i / numSamples, radicalInverse);
    
}



float3 cartesian_to_spherical(float3 cartesian)
{
    float3 cartesian2 = cartesian * cartesian;
    float phi = acos(cartesian.z / (sqrt(cartesian2.x + cartesian2.y + cartesian2.z)));
    return float3(sqrt(cartesian2.x + cartesian2.y + cartesian2.z), atan(cartesian.y / cartesian.x), phi);

}


//TODO: Probably drop TIR for being too expensive.
float3 TIR_lookup(float3 coords)
{
    return TIR_LUT.Sample(defaultSampler, coords);
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
    float2 splitsum = FGD_LUT.Sample(defaultSampler, float2(alpha, cti));
    
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
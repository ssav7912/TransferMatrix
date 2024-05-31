#pragma once

#include "HenyeyGreenstein.hlsli"
#include "tensor3d.hlsli"
#include "../Common.hlsli"

static const uint TM2_SAMPLE_TYPE_GLOSSY_REFLECTION = 0;
static const uint TM2_SAMPLE_TYPE_GLOSSY_TRANSMISSION = 1;


static const uint TM2_MEASURE_SOLID_ANGLE = 0;


static const uint TM2_TYPE_NOCOMPONENT = 0;
static const uint TM2_TYPE_DIELECTRICINTERFACE = 1;
static const uint TM2_TYPE_CONDUCTORINTERFACE = 2;


#define PDF_DEBUG float3(0.0f, 0.0f, 0.0f)
#define EVAL_DEBUG float3(0.0f, 0.0f, 0.0f)
#define NAN_DEBUG float3(0.0f, 0.0f, 0.0f)

#define NUM_SAMPLES 100

#define LAYERS_MAX 5
#define NUM_LAYERS 2
#define NUM_LOBES (NUM_LAYERS + 1)

#define NO_SECOND_UV 1

#define EPSILON 1e-6f


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

struct IBLLight
{
    float3 radiance;
    float3 irradiance;
};

float float3_average(float3 f)
{
    return (f.x + f.y + f.z) / 3.0f;

}

float copy_sign(float x, float s)
{
    return (s >= 0) ? abs(x) : -abs(x);

}

float3 reflectSpherical(float3 incident, float3 normal)
{
    return 2 * dot(incident, normal) * normal - incident;
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

float sinTheta2(float3 v)
{
    return 1.0f - v.z * v.z;
}

float cosTheta2(float3 v)
{
    return v.z * v.z;
}

//based off mitsuba 3 implementation
float fresnelDielectric(float incidentCosTheta, float ior)
{
    float transmittedCosTheta;
    if (ior == 1.0f)
    {
        return 0.0f;
    }
    
    float scale = (incidentCosTheta > 0) ? 1.0f / ior : ior;
    float transmittedcosTheta2 = 1.0f - (1.0f - incidentCosTheta * incidentCosTheta) * (scale * scale);
    
    if (transmittedcosTheta2 <= 0.0f)
    {
        return 1.0f;
    }
    
    incidentCosTheta = abs(incidentCosTheta);
    transmittedCosTheta = sqrt(transmittedcosTheta2);
    
    float Rs = (incidentCosTheta - ior * transmittedCosTheta) / (incidentCosTheta + ior * transmittedCosTheta);
    float Rp = (ior * incidentCosTheta - transmittedCosTheta) / (ior * incidentCosTheta + transmittedCosTheta);
    
    return 0.5f * (Rs * Rs + Rp * Rp);
}

float3 fresnelConductorExact(float incidentCosTheta, float3 ior, float3 kappa)
{
    float incidentCosTheta2 = incidentCosTheta * incidentCosTheta;
    float sinTheta2 = 1 - incidentCosTheta2;
    float sinTheta4 = sinTheta2 * sinTheta2;
    
    float3 temp1 = ior * ior - kappa * kappa - sinTheta2.xxx;
    float3 a2pb2 = sqrt((temp1 * temp1 + kappa * kappa * ior * ior * 4.0f.xxx));
    float3 a = sqrt(((a2pb2 + temp1) * 0.5f));
    
    float3 term1 = a2pb2 + incidentCosTheta2.xxx;
    float3 term2 = a * (2 * incidentCosTheta);
    
    float3 Rs2 = (term1 - term2) / (term1 + term2);
    
    float3 term3 = a2pb2 * incidentCosTheta2 + sinTheta4.xxx;
    float3 term4 = term2 * sinTheta2;
    
    float3 Rp2 = Rs2 * (term3 - term4) / (term3 + term4);
    
    return 0.5f * (Rp2 + Rs2);
}

float3 fresnelConductorApprox(float cosThetaI, float3 ior, float3 kappa)
{
    float cosThetaI2 = cosThetaI * cosThetaI;
    
    float3 temp = (ior * ior + kappa * kappa) * cosThetaI2;
    
    float3 Rp2 = (temp - (ior * (2.0f * cosThetaI) + 1.0f.xxx)) / (temp + (ior * (2 * cosThetaI) + 1.0f.xxx));
    
    float3 tmpF = ior * ior + kappa * kappa;
    
    float3 Rs2 = (tmpF - (ior * (2.0f * cosThetaI) + cosThetaI2.xxx)) / (tmpF + (ior * (2.0f * cosThetaI2) + cosThetaI2.xxx));
    
    return 0.5f * (Rp2 + Rs2);
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


float smithG(float3 incident, float3 outgoing, float3 m, float rough)
{
    return smithG1(incident, m, rough) * smithG1(outgoing, m, rough);
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



///Translate cartesian to spherical coordinates. Note the cartesian basis is Y-up.
float3 cartesian_to_spherical(float3 cartesian)
{
    float InUp = cartesian.y;
    float InRight = cartesian.x;
    float InForward = cartesian.z;
    
    float3 cartesian2 = cartesian * cartesian;
    //float theta = atan2(cartesian.y, cartesian.x);
    //float phi = atan2(sqrt(cartesian2.x + cartesian2.y + cartesian2.z), cartesian.z);
    float elevation = atan2(InRight, InForward);
    float azimuth = atan2(sqrt(InForward * InForward + InRight * InRight), InUp);
    
    return float3(sqrt(cartesian2.x + cartesian2.y + cartesian2.z), azimuth, elevation);

}

///Takes coordinates of form (radius, Azimuth, Elevation), returns (Right, Up, Forward).
float3 spherical_to_cartesian(float3 spherical)
{
    float InRadius = spherical.x;
    float InElevation = spherical.z;
    float InAzimuth = spherical.y;
    
    float SinAzimuth = sin(InAzimuth);
    float CosAzimuth = cos(InAzimuth);
    
    float CosElevation = cos(InElevation);
    float SinElevation = sin(InElevation);
   
    
    float Forward = InRadius * CosElevation * SinAzimuth;
    float Right = InRadius * SinElevation * SinAzimuth;
    float Up = InRadius * CosAzimuth;
    
    return float3(Right, Up, Forward);
}



//takes tangent space cartesian coords and transforms it to mitsuba local space (spherical?)
float3 cartesianTSToMitsubaLS(float3 cartesian)
{
    const float3 s = float3(1, 0, 0);
    const float3 t = float3(0, 0, 1);
    const float3 n = float3(0, 1, 0);
    return float3(dot(cartesian, s), dot(cartesian, n), dot(cartesian, t));

}

float3 MitsubaLSToCartesianTS(float3 mitsuba)
{
    const float3 s = float3(1, 0, 0);
    const float3 t = float3(0, 0, 1);
    const float3 n = float3(0, 1, 0);
    
    return s * mitsuba.x + n * mitsuba.y + t * mitsuba.z;
}


//TODO: Probably drop TIR for being too expensive.
float3 TIR_lookup(float3 coords)
{
    return TIR_LUT.Sample(defaultSampler, coords);
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


//Lagarde 2011. Compute fresnel reflectance at 0 degrees from IOR
float3 f0(float3 ior)
{
    return ((ior - 1) * (ior - 1)) / ((ior + 1) * (ior + 1));
}

void albedo(float cti, float alpha, float3 ior_ij, float3 kappa_ij, out float3 r_ij)
{
    //TODO: investigates
    float2 FGD = FGD_LUT.Sample(defaultSampler, float2(alpha, cti)).x;
    
    r_ij = FGD.xxx + f0(ior_ij) * FGD.y;

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

void conductor_transfer_factors(float3 incident, float3 ior, float3 kappa, float rough, out layer_components_tm2 ops)
{
    ops.component_type = TM2_TYPE_CONDUCTORINTERFACE;
    
    ops.reflection_down.asymmetry = ggx_to_hg(rough);
    ops.reflection_down.mean = reflectZ(incident);
    
    albedo(abs(incident.z), rough, ior, kappa, ops.reflection_down.norm);
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
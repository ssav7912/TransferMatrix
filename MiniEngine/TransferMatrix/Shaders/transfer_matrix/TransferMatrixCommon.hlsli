#pragma once
#include "../Common.hlsli"
#include "HenyeyGreenstein.hlsli"
#include "Constants.hlsli"



#define PDF_DEBUG float3(0.0f, 0.0f, 0.0f)
#define EVAL_DEBUG float3(0.0f, 0.0f, 0.0f)
#define NAN_DEBUG float3(1.0f, 0.0f, 1.0f)

#define NUM_SAMPLES 5

#define LAYERS_MAX 5

#define NO_SECOND_UV 1

#define EPSILON 1e-6f
#define HALFEPSILON 2E-10

#define IOR_AIR 1.003.xxx
#define KAPPA_AIR 0.0.xxx
#define ALPHA_AIR 0.0


//Preprocessor switch validations,
#if USE_KARIS_FGD == 1 && USE_BELCOUR_FGD == 1
    #error "Only 1 of USE_KARIS_FGD or USE_BELCOUR_FGD can be enabled, not both!"
#endif

#if ANALYTIC_TIR == 1 && DISABLE_TIR == 1
    #error "Only 1 of USE_ANALYTIC_TIR or DISABLE_TIR can be enabled, not both!" 
#endif

#if USE_EARL_G2 == 1 && SCHLICK_G == 1
    #error "Only 1 of USE_EARL_G2 or SCHLICK_G can be enabled, not both!"

#endif


struct sample_record
{
    float3 incident;
    float3 outgoing;
};




//Texture Arrays for textured impl.
//Texture2DArray<float3> TexIORs : register(t0);
//Texture2DArray<float3> TexKappas : register(t5);
//#ifdef TM6
//Texture2DArray<float3> TexSigma_S : register(t10);
//Texture2DArray<float3> TexSigma_K : register(t15);
//Texture2DArray<float> TexDepths : register(t20);
//Texture2DArray<float> TexPhase : register (t25);
//#endif
//Texture2DArray<float> TexRoughs : register(t30);


//Lookup table for Total Internal Reflection
Texture3D<real> TIR_LUT : register(t43);

//compile time switch for (Karis 2013). split sum FGD approx, the Belcour (2020) splitsum approx, or the (Belcour 2018) FGD LUT.
#if USE_KARIS_FGD == 1
//LUT for Karis Split-sum FGD approximation.
Texture2D<float2> FGD_LUT : register(t44);

#elif USE_BELCOUR_FGD == 1
Texture2D<float4> FGD_Belcour_LUT : register(t47);

#else
Texture3D<float> FGD_LUT : register(t45);
#endif

Texture2D<float> GD_LUT : register(t46);


//IBL
TextureCube<float3> radianceIBLTexutre : register(t35);
TextureCube<float3> irradianceIBLTexture : register(t35);

sampler LUTSampler : register(s13);

//Layer Parameters
cbuffer LayerParameters : register(b2)
{
    float3 IORs[LAYERS_MAX];
    float3 Kappas[LAYERS_MAX];
    //used for TM6
    float3 Sigma_S[LAYERS_MAX];
    float3 Sigma_K[LAYERS_MAX];
    float Depths[LAYERS_MAX];
    float G[LAYERS_MAX];
    float Roughs[LAYERS_MAX];
    //common
    double PADDING;
    int NumLayers;
    int NumSamples;
}

struct LayerProperties
{
    real3 iors[LAYERS_MAX];
    real3 kappas[LAYERS_MAX];
    real rough[LAYERS_MAX];
#ifdef TM6
    real depths[LAYERS_MAX];
    real3 sigma_s[LAYERS_MAX];
    real3 sigma_k[LAYERS_MAX];
    real gs[LAYERS_MAX];
#endif
    int NumLayers;
};


LayerProperties zero_initialise_layers()
{
    LayerProperties x;
    for (int i = 0; i < LAYERS_MAX; i++)
    {
        x.iors[i] = 0.0.xxx;
        x.kappas[i] = 0.0.xxx;
#ifdef TM6
         x.sigma_s[i] = 0.0.xxx;
        x.sigma_k[i] = 0.0.xxx;
        x.depths[i] = 0.0;
        x.gs[i] = 0.0;
#endif
        x.rough[i] = 0.0;
    }
    return x;

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


//populates the layer properties structure for the sampled UV coordinate using either 
//the texture arrays or the uniform buffer, depending on the provided layermask.
//if the bit in the i'th digit of layermask is set, this means sample from the relevant texturearray.
//LayerProperties sample_layer_textures(float2 uv, uint layermask)
//{
//    LayerProperties x = zero_initialise_layers();
//    for (int i = 0; i <= NumLayers + 1; i++) //indexing starts at 1,
//    {
//        x.iors[i] = (1 << i) & layermask ? TexIORs.Sample(defaultSampler, float3(uv, i)) : IORs[i];
//        x.kappas[i] = (1 << i) & layermask ? TexKappas.Sample(defaultSampler, float3(uv, i)) : Kappas[i];
//#ifdef TM6
//        x.sigma_s[i] = (1 << i) & layermask ? TexSigma_S.Sample(defaultSampler, float3(uv, i)) : Sigma_S[i];
//        x.sigma_k[i] = (1 << i) & layermask ? TexSigma_K.Sample(defaultSampler, float3(uv, i)) : Sigma_K[i];
//        x.depths[i] = (1 << i) & layermask ? TexDepths.Sample(defaultSampler, float3(uv, i)) : Depths[i];
//        x.gs[i] = (1 << i) & layermask ? TexPhase.Sample(defaultSampler, float3(uv, i)) : G[i];
//#endif        
//        x.rough[i] = (1 << i) & layermask ? TexRoughs.Sample(defaultSampler, float3(uv, i)) : Roughs[i];
//
//    }
//    
//    x.NumLayers = NumLayers;
//    return x;
//
//}


//helper function to truncate full precision constant buffer to
//half precision layerprops struct, workaround for type mismatch.
LayerProperties truncate_layer_parameters()
{
    LayerProperties props;
    
    for (int i = 0; i < LAYERS_MAX; i++)
    {
        props.iors[i] = (real3) IORs[i];
        props.kappas[i] = (real3) Kappas[i];
#ifdef TM6
        props.sigma_s[i] = (real3)Sigma_S[i];
        props.sigma_k[i] = (real3)Sigma_K[i];
        props.gs[i] = (real)G[i];
        props.depths[i] = (real)Depths[i];
#endif
        
        props.rough[i] = (real) Roughs[i];
    }

    return props;
}


real hg_lh_norm(real g)
{
    const bool g_neg = g < 0.;
    g = abs(g);
    const float n = clamp(0.5039 - 0.8254 * g + 0.3226 * g * g, 0., 0.5);
    return g_neg ? 1.0 - n : n;
}

real safe_div(real n, real d)
{
    return abs(d) > EPSILON ? n / d : 0.;
}

real3 safe_div(real3 n, real3 d)
{
    float3 r;
    r.x = abs(d.x) > EPSILON ? n.x / d.x : 0.0;
    r.y = abs(d.y) > EPSILON ? n.y / d.y : 0.0;
    r.z = abs(d.z) > EPSILON ? n.z / d.z : 0.0;
    return r;
}

bool isZero(real3 vec)
{
    return vec.x + vec.y + vec.z == 0.0;
}


real float3_max(real3 v)
{
    return max(v.x, max(v.y, v.z));
}

real float3_average(real3 f)
{
    
    return (f.x + f.y + f.z) / 3.0;

}


real copy_sign(real x, real s)
{
    return (s >= 0) ? abs(x) : -abs(x);

}


float3 reflectSpherical(float3 incident, float3 normal)
{
    return 2 * dot(incident, normal) * normal - incident;
}

//min16float3 reflectZ(min16float3 f)
//{
//    return min16float3(-f.xy, f.z);
//}

float3 reflectZ(float3 f)
{
    return float3(-f.xy, f.z);
}

//min16float3 refractZ(min16float3 f, min16float ior)
//{
//    return refract(f, min16float3(0.0, 0.0, copy_sign(1.0, f.z)), ior);
//}

float3 refractMitsuba(float3 wi, float3 n, real ior)
{
    if (ior == 1.0)
    {
        return -wi;
    }
    
    float cosThetaI = dot(wi, n);
    if (cosThetaI > 0)
    {
        ior = 1 / ior;
    }
    
    float cosThetaTsqr = 1 - (1 - cosThetaI * cosThetaI) * (ior * ior);
    
    if (cosThetaTsqr <= 0.0)
    {
        return 0.0.xxx;
    }
    
    return n * (cosThetaI * ior - sign(cosThetaI) * sqrt(cosThetaTsqr)) - wi * ior;
}

float3 refractZ(float3 f, real ior)
{
    return refractMitsuba(f, float3(0.0, 0.0, copy_sign(1.0, f.z)), ior);

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


float Pow5(float x)
{
    float xSq = x * x;
    return xSq * xSq * x;
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

//mitsuba standard basis
static const float3 S = float3(1, 0, 0);
static const float3 T = float3(0, 1, 0);
static const float3 N = float3(0, 0, 1);

//takes tangent space cartesian coords and transforms it to mitsuba local space (spherical?)
float3 cartesianTSToMitsubaLS(float3 cartesian)
{
    cartesian = cartesian.xzy;
    
    
    return float3(dot(cartesian, S), dot(cartesian, N), dot(cartesian, T));

}

float3 MitsubaLSToCartesianTS(float3 mitsuba)
{

    
    return (S * mitsuba.x + N * mitsuba.y + T * mitsuba.z).xzy;
}



//Lagarde 2011. Compute fresnel reflectance at 0 degrees from IOR
float3 f0(float3 ior)
{
    return ((ior - 1.0f) * (ior - 1.0f)) / ((ior + 1.0f) * (ior + 1.0f));
}

//based off mitsuba 3 implementation
real fresnelDielectric(float incidentCosTheta, real ior)
{
    if (ior == 1.0)
    {
        return 0.0;
    }
    
    real scale = (incidentCosTheta > 0) ? 1.0 / ior : ior;
    real transmittedcosTheta2 = 1.0 - (1.0 - incidentCosTheta * incidentCosTheta) * (scale * scale);
    
    if (transmittedcosTheta2 <= 0.0)
    {
        return 1.0;
    }
    
    incidentCosTheta = abs(incidentCosTheta);
    real transmittedCosTheta = sqrt(transmittedcosTheta2);
    
    real Rs = (incidentCosTheta - ior * transmittedCosTheta) / (incidentCosTheta + ior * transmittedCosTheta);
    real Rp = (ior * incidentCosTheta - transmittedCosTheta) / (ior * incidentCosTheta + transmittedCosTheta);
    
    return 0.5 * (Rs * Rs + Rp * Rp);
}

float3 fresnelDielectric(float incidentCosTheta, float3 ior)
{
    //if ior == 1 early exit
   
    
    float3 scale = (incidentCosTheta > 0) ? 1.0.xxx / ior : ior;
    float3 transmittedcosTheta2 = 1.0.xxx - (1.0 - incidentCosTheta * incidentCosTheta) * (scale * scale);

    
    incidentCosTheta = abs(incidentCosTheta);
    float3 transmittedCosTheta = sqrt(transmittedcosTheta2);

    float3 Rs = (incidentCosTheta - ior * transmittedCosTheta) / (incidentCosTheta + ior * transmittedCosTheta);
    float3 Rp = (ior * incidentCosTheta - transmittedCosTheta) / (ior * incidentCosTheta + transmittedCosTheta);
   
    float3 output = 0.5 * (Rs * Rs + Rp * Rp);
   
    bool3 iormask = ior == 1.0;
    bool3 tcti2mask = transmittedcosTheta2 <= 0.0;
    
    output = lerp(output, 1.0.xxx, tcti2mask); //sets values where tcti2 <= 0 to 1.
    output = lerp(output, 0.0.xxx, iormask); //sets values where ior == 1 to 0.
    
    
    return output;
    
}

// Shlick's approximation of Fresnel
float3 Fresnel_Shlick(float3 F0, float3 F90, float cosine)
{
    return lerp(F0, F90, Pow5(1.0 - cosine));
}

float Fresnel_Shlick(float F0, float F90, float cosine)
{
    return lerp(F0, F90, Pow5(1.0 - cosine));
}

real3 fresnelConductorExact(real incidentCosTheta, real3 ior, real3 kappa)
{
    real incidentCosTheta2 = incidentCosTheta * incidentCosTheta;
    real sinTheta2 = 1.0 - incidentCosTheta2;
    real sinTheta4 = sinTheta2 * sinTheta2;
    
    real3 temp1 = (ior * ior) - (kappa * kappa) - sinTheta2.xxx;
    real3 a2pb2 = sqrt(((temp1 * temp1) + (kappa * kappa) * (ior * ior) * 4.0.xxx));
    real3 a = sqrt(((a2pb2 + temp1) * 0.5));
    
    real3 term1 = a2pb2 + incidentCosTheta2.xxx;
    real3 term2 = a * (2.0 * incidentCosTheta);
    
    real3 Rs2 = (term1 - term2) / (term1 + term2);
    
    real3 term3 = a2pb2 * incidentCosTheta2.xxx + sinTheta4.xxx;
    real3 term4 = term2 * sinTheta2;
    
    real3 Rp2 = Rs2 * (term3 - term4) / (term3 + term4);
    
    return 0.5 * (Rp2 + Rs2);
}

float3 fresnelConductorApprox(float cosThetaI, float3 ior, float3 kappa)
{
    float cosThetaI2 = cosThetaI * cosThetaI;
    
    float3 temp = (ior * ior + kappa * kappa) * cosThetaI2;
    
    float3 Rp2 = (temp - (ior * (2.0 * cosThetaI) + 1.0.xxx)) / (temp + (ior * (2 * cosThetaI) + 1.0.xxx));
    
    float3 tmpF = ior * ior + kappa * kappa;
    
    float3 Rs2 = (tmpF - (ior * (2.0 * cosThetaI) + cosThetaI2.xxx)) / (tmpF + (ior * (2.0 * cosThetaI2) + cosThetaI2.xxx));
    
    return 0.5 * (Rp2 + Rs2);
}


//Karis [2013]
real SchlickG1(float3 v, real alpha)
{
    real k = ((alpha + 1.0) * (alpha + 1.0) / 8.0);
    float3 N = float3(0, 0, 1);
    float NdotV = saturate(dot(N, v));
    
    return NdotV / (NdotV * (1.0 - k) + k);
    
}


real SchlickG2(float3 l, float3 v, float3 h, real alpha)
{
    return SchlickG1(l, alpha) * SchlickG1(v, alpha);

}

real smithG1(float NdotV, real alpha)
{
    real a2 = alpha * alpha;
    return 2.0 / (1.0 + sqrt(1.0 + a2 * (1.0 - NdotV * NdotV)) / (NdotV * NdotV));

}

//[Lagarde & De Rousiers, 2014]
//masking-shadowing for G = G1*G2 (incoming and outgoing)
real smithGCorrelated(real NdotL, real NdotV, real alpha)
{
    NdotL = saturate(NdotL);
    NdotV = saturate(NdotV);
    
    real alpha2 = alpha * alpha;
    real LambdaV = NdotL * sqrt((-NdotV * alpha2 + NdotV) * NdotV + alpha2);
    real LambdaL = NdotV * sqrt((-NdotL * alpha2 + NdotL) * NdotL + alpha2);
    
    return 0.5 / (LambdaV + LambdaL);
}

//[Pharr et al 2023]
real G2Correlated(real3 incident, real3 outgoing, real alpha)
{
    real3 H = float3(0, 0, 1);
    
    if (dot(incident, H) * incident.z <= 1e-05 || dot(outgoing, H) * outgoing.z <= 1e-05)
    {
        return 0.0;
    }
    
    if (tanTheta(incident) <= 1e-05 || tanTheta(outgoing) <= 1e-05)
    {
        return 1.0;
    }
    
    
    real lambdaI = (sqrt(1 + alpha * alpha * tanTheta(incident) * tanTheta(incident)) - 1.0)/2.0;
    real lambdaO = (sqrt(1 + alpha * alpha * tanTheta(outgoing) * tanTheta(outgoing)) - 1.0) / 2.0;
    
    return 1.0 / (1.0 + lambdaO + lambdaI);
   
}

real G1Correlated(real3 incident,real alpha)
{   
    real3 H = float3(0, 0, 1);
    if (dot(incident, H) * incident.z <= 1e-05)
    {
        return 0.0;
    }
    
    if (tanTheta(incident) <= 1e-05)
    {
        return 1.0;
    }
    real lambdaI = (sqrt(1 + alpha * alpha * tanTheta(incident) * tanTheta(incident)) - 1.0) / 2.0;
    
    return 1.0 / (1.0 + lambdaI);
   
}

//assume isotropic
real smithG1(float3 v, float3 m, real alpha)
{
    alpha = clamp(alpha, EPSILON, 1.0);
    
    //v.z == Frame::cosTheta
    if (dot(v, m) * v.z <= 0.0)
    {
         return 0.0;
    }
    
    float tantheta = abs(tanTheta(v));
    if (tantheta <= 0.0) //TODO: almost equal op...
    {
        return 1.0;
    }
    
    
    real root = alpha * tantheta;
                    //hypot2
    return (2.0 / (1.0 + (1.0 * 1.0 + root * root)));

}


float D_GGX(float3 m, real alpha)
{
   if (m.z <= 0.0)
    {
        return 0.0;
    }
    alpha = max(EPSILON, alpha);
    
    float costheta2 = cosTheta2(m);
    real beckmann = (((m.x * m.x) / (alpha * alpha) + (m.y * m.y) / (alpha * alpha)) / costheta2);
    
    real root = (1.0 + beckmann) * costheta2;
    return (1.0 / (PI * alpha * alpha * root * root));

}

//Karis 2013 distribution.
real D_GGX_Karis(float NdotH, real rough)
{
    rough = max(EPSILON, rough);
    real rough_square = rough * rough;
    real  f = (NdotH * NdotH) * (rough_square - 1.0) + 1.0;
    return rough_square / (PI * f * f);
}

real smithG(float3 incident, float3 outgoing, float3 m, real rough)
{
    return smithG1(incident, m, rough) * smithG1(outgoing, m, rough);
}


//https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2017/Presentations/Hammon_Earl_PBR_Diffuse_Lighting.pdf
//Earl H. (2017)


real smithG2(float3 V, float3 N, float3 L, real rough)
{
    real numerator = 2.0f * (abs(dot(N, L)) * abs(dot(N, V)));
    real denominator = lerp(2.0f * abs(dot(N, L)) * abs(dot(N, V)), abs(dot(N, L)) + abs(dot(N, V)), rough);
    return numerator / denominator;
}


// [Lagarde & De Rousiers, 2014] Moving Frostbite to Physically Based Rendering 3.0
//Lobe off-peak shift correction (Section 4.9.3)
float3 lobe_mean_shift(float3 incident, float masking_shadowing, float3 Normal)
{
    return masking_shadowing * Normal + (1.0 - masking_shadowing) * incident;
}

//[Lagarde & De Rousiers, 2014] correlated GSmith specular dominant.
float3 specular_dominant(float3 normal, float3 outgoing, float NdotV, real rough)
{
    //correlated version
    float lerpFactor = pow(1 - NdotV, 10.8649) * (1 - 0.298475 * log(39.4115 - 39.0029 * rough))
        + 0.298475 * log(39.4115 - 39.0029 * rough);
    
    return lerp(normal, outgoing, lerpFactor);
}

// [Lagarde & De Rousiers, 2014]
float3 specular_dominant_approx(float3 normal, float3 outgoing, float rough)
{
    float smooth = saturate(1.0 - rough);
    float lerpFactor = smooth * (sqrt(smooth) + rough);
    return lerp(normal, outgoing, lerpFactor);
}

//TODO: Probably drop TIR for being too expensive.
real TIR_lookup(float3 coords)
{
    return TIR_LUT.Sample(LUTSampler, coords);
}

real TIR_analytical(float cti, real rough, real ior_12, real ior_10)
{
    float3 wo = 0.xxx;
    wo.x = sin(0.25 * PI * (1.0 + cti));
    wo.z = sqrt(1.0 - (wo.x * wo.x));
    
    real alpha = max(rough, 0.01);
    float3 wi = 0.xxx;
    wi.z = cti;
    
    float3 normal = float3(0, 0, 1);
    
    float3 off_specular = specular_dominant(normal, wo, dot(normal, wi), alpha);
    
    
    float3 ws = -reflect(wo, off_specular);
    
    
    real value = 0;
   
    real R = fresnelDielectric(dot(wo, off_specular), ior_12);
#if USE_SMITH_G2 == 1 
    value += R * smithG1(ws, normal, alpha);
#else
    value += R * G1Correlated(ws, alpha);
#endif
    value *= (1.0 - fresnelDielectric(wo.z, ior_10)) * step(ws.z, 0.0); //branch elimination, make value 0 if ws.z <= 0 
             
   

    return value;
}


//[Belcour 2020]
//calculate coefficients for split-sum LUT
real4 fresnel_to_coefficients(real ior, real kappa)
{
    //critical angles used to constrain the basis functions
    static const float ct1 = 0.12812812812812813;
    static const float ct2 = 0.43243243243243246;
    
    //basis functions sampled at ct1 and ct2.
    static const float i1_basis1 = 0.4268566850109629;
    static const float i2_basis1 = 0.07262985737860257;
    static const float i1_basis2 = 1.0;
    static const float i2_basis2 = 0.4435120658488671;
    static const float i1_basis3 = 0.6181475431224497;
    static const float i2_basis3 = -0.640692996763135;
    
    real4 coeffs = 0.xxxx;
    real2 dF = 0.xx;
    
    real fresnel1 = kappa == 0.0 ? fresnelDielectric(1.0, float3_average(ior)) : fresnelConductorExact(1.0, ior.xxx, kappa.xxx).x;
    real fresnelct1 = kappa == 0.0 ? fresnelDielectric(ct1, float3_average(ior)) : fresnelConductorExact(ct1, ior.xxx, kappa.xxx).x;
    real fresnelct2 = kappa == 0.0 ? fresnelDielectric(ct2, float3_average(ior)) : fresnelConductorExact(ct2, ior.xxx, kappa.xxx).x;
    
    coeffs.x = fresnel1;
    coeffs.y = 1.0 - coeffs.x;

    dF.x = fresnelct1 - (coeffs.x + coeffs.y * i1_basis1);
    dF.y = fresnelct2 - (coeffs.x + coeffs.y * i2_basis1);
    
    /*
    float2x2 F12 = float2x2(i1_basis2, i2_basis2,
                            i1_basis3, i2_basis3);
        
    A is precomputed inverse of F12 above ^ 
    */
    real2x2 A = real2x2(0.70032658, 0.4847927,
                          0.67568267, -1.09307669); 
    
    coeffs.zw = mul(coeffs.zw, A);
    
    return coeffs;

}

real3 sample_FGD(float cti, real alpha, real3 ior, real3 kappa)
{
    real3 output = NAN_DEBUG;
#if USE_KARIS_FGD == 1
    //FGD 4D LUT parameterised by elevation, roughness and complex IOR
    //(Karis 2013) approximation parameterised by Roughness and angle
    //outputs scale & bias to F0.
    //TODO: work out 2D split sum approx.
    //https://cdn2.unrealengine.com/Resources/files/2013SiggraphPresentationsNotes-26915738.pdf
    //Ok split sum approx contains FGD_1 & FGD_2 (whatever that means)
    //tristimulus energy output?
    real2 splitsum = FGD_LUT.Sample(LUTSampler, float2(cti, alpha));
    
    //this crushes up the energy for dielectrics...
    //why doesn't this warn?
     const real3 F0 = f0(ior);
    
     output.x = (F0.x * splitsum.x + (1-F0.x) * splitsum.y);
     output.y = (F0.y * splitsum.x + (1-F0.y) * splitsum.y);
     output.z = (F0.z * splitsum.x + (1-F0.z) * splitsum.y);
    
    return output;
     
#elif USE_BELCOUR_FGD == 1
    //Belcour (2020) split sum model for complex IOR.

    const real4 coeffsX = fresnel_to_coefficients(ior.x, kappa.x);
    const real4 coeffsY = fresnel_to_coefficients(ior.y, kappa.y);
    const real4 coeffsZ = fresnel_to_coefficients(ior.z, kappa.z);
    
    real4 splitsum = FGD_Belcour_LUT.Sample(LUTSampler, float2(cti, alpha));
    
    output.x = (coeffsX.x * splitsum.x) + (coeffsX.y * splitsum.y) + (coeffsX.z * splitsum.z) + (coeffsX.w * splitsum.w);
    output.y = (coeffsY.x * splitsum.x) + (coeffsY.y * splitsum.y) + (coeffsY.z * splitsum.z) + (coeffsY.w * splitsum.w);
    output.z = (coeffsZ.x * splitsum.x) + (coeffsZ.y * splitsum.y) + (coeffsZ.z * splitsum.z) + (coeffsZ.w * splitsum.w);
    
    return output;
    
#else
    static const uint NumDimensions = 4;
    static const uint DimensionXY = 64;
    static const uint DimensionZ = DimensionXY * 32;
    static const uint DimensionW = 32;
    
    static const float DimensionXMin = 0.0;
    static const float DimensionXMax = 1.0;
    static const float DimensionYMin = 0.0;
    static const float DimensionYMax = 1.0;
    static const float DimensionZMin = 0.0;
    static const float DimensionZMax = 4.0;
    static const float DimensionWMin = 0.0;
    static const float DimensionWMax = 4.0;
    
   
    //remap ranges to 0-1 normalised index.
    const float cti_index = saturate((cti - DimensionXMin) / (DimensionXMax - DimensionXMin));
    const float alpha_index = saturate((alpha - DimensionYMin) / (DimensionYMax - DimensionYMin));
    
    //remap Z & W index from being within [0,64] range to [0,2048]
    //Z index is normalised [0,64] range stretched to [0,2048]
    const float3 ior_norm = (ior - DimensionZMin) / (DimensionZMax - DimensionZMin);
    
    const float3 ior_index = (ior_norm * (DimensionZ - 1));
    
    //W index is normalised [0,32] range + Z index.
    const float3 kappa_norm = (kappa - DimensionWMin) / (DimensionWMax - DimensionWMin);
    
    const float3 kappa_index = (kappa_norm * (DimensionW - 1));
   
    const float3 ZIndex = (ior_index + kappa_index) / (DimensionZ - 1);

    
    //4D LUT is packed into 3D Texture, z & w share the same dimension. 
    output.x = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.x));
    output.y = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.y));
    output.z = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.z));
#endif
    
    return output;
}

//could i just reuse the FGD LUT for this?
float sample_GD(float cti, float rough)
{
    return GD_LUT.Sample(LUTSampler, float2(cti, rough));
}



void albedos(float cti, real alpha, real ior_ij, out real3 r_ij, out real3 t_ij, out real3 r_ji, out real3 t_ji)
{
    if (abs(ior_ij - 1.0) < 1e-3)
    {
        r_ij = 0.0.xxx;
        r_ji = 0.0.xxx;
        
        t_ij = 1.0.xxx;
        t_ji = 1.0.xxx;
        return;
    }
    

    r_ij = sample_FGD(cti , alpha, ior_ij.xxx, 0.0.xxx);
    
    
    r_ji = r_ij;
    
    t_ij = 1.0.xxx - r_ij; //can't be negative
    t_ji = t_ij;
    
}



void albedo(float cti, real alpha, real3 ior_ij, real3 kappa_ij, out real3 r_ij)
{
    //can't be negative! 
    
    r_ij = sample_FGD(cti, alpha, ior_ij, kappa_ij);

}

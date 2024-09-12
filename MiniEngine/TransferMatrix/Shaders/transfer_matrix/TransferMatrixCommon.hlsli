#pragma once
#include "../Common.hlsli"
#include "HenyeyGreenstein.hlsli"
#include "Constants.hlsli"



#define PDF_DEBUG float3(0.0f, 0.0f, 0.0f)
#define EVAL_DEBUG float3(0.0f, 0.0f, 0.0f)
#define NAN_DEBUG float3(1.0f, 0.0f, 1.0f)

#define NUM_SAMPLES 5

#define LAYERS_MAX 5
#define NUM_LAYERS 2
#define NUM_LOBES (NUM_LAYERS + 1)

#define NO_SECOND_UV 1

#define EPSILON 1e-6f
#define HALFEPSILON 2E-10



struct sample_record
{
    float3 incident;
    float3 outgoing;
    float ior;
    bool is_reflection_sample;
    bool sample_type;
};


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

//compile time switch for (Karis 2013). split sum FGD approx, or the (Belcour 2018). FGD LUT.
#if USE_FAST_FGD == 1
//LUT for Karis Split-sum FGD approximation.
Texture2D<float2> FGD_LUT : register(t44);
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
    
    [unroll]
    for (int i = 0; i < LAYERS_MAX; i++)
    {
        props.iors[i] = IORs[i];
        props.kappas[i] = Kappas[i];
#ifdef TM6
        props.sigma_s[i] = Sigma_S[i];
        props.sigma_k[i] = Sigma_K[i];
        props.gs[i] = G[i];
        props.depths[i] = Depths[i];
#endif
        
        props.rough[i] = Roughs[i];
    }

    return props;
}


min16float hg_lh_norm(min16float g)
{
    const bool g_neg = g < 0.f;
    g = abs(g);
    const min16float n = clamp(0.5039 - 0.8254 * g + 0.3226 * g * g, 0.0, 0.5);
    return g_neg ? 1.0 - n : n;
}

float hg_lh_norm(float g)
{
    const bool g_neg = g < 0.f;
    g = abs(g);
    const float n = clamp(0.5039f - 0.8254f * g + 0.3226f * g * g, 0.f, 0.5f);
    return g_neg ? 1.0f - n : n;
}

min16float safe_div(min16float n, min16float d)
{
    return abs(d) > EPSILON ? n / d : 0.0;
}

float safe_div(float n, float d)
{
    return abs(d) > EPSILON ? n / d : 0.f;
}

min16float3 safe_div(min16float3 n, min16float3 d)
{
    
    min16float3 r;
    //EPSILON should be 16bit float eps...
    r.x = abs(d.x) > EPSILON ? n.x / d.x : 0.0;
    r.y = abs(d.y) > EPSILON ? n.y / d.y : 0.0;
    r.z = abs(d.z) > EPSILON ? n.z / d.z : 0.0;
    return r;
    
}

float3 safe_div(float3 n, float3 d)
{
    float3 r;
    r.x = abs(d.x) > EPSILON ? n.x / d.x : 0.f;
    r.y = abs(d.y) > EPSILON ? n.y / d.y : 0.f;
    r.z = abs(d.z) > EPSILON ? n.z / d.z : 0.f;
    return r;
}

bool isZero(min16float3 vec)
{
    return vec.x + vec.y + vec.z == 0.0;
}

bool isZero(float3 vec)
{
    return vec.x + vec.y + vec.z == 0.0f;
}

min16float float3_max(min16float3 v)
{
    return max(v.x, max(v.y, v.z));
}

float float3_max(float3 v)
{
    return max(v.x, max(v.y, v.z));
}

float float3_average(real3 f)
{
    
    return (f.x + f.y + f.z) / 3.0;

}

min16float copy_sign(min16float x, min16float s)
{
    return (s >= 0) ? abs(x) : -abs(x);
}

float copy_sign(float x, float s)
{
    return (s >= 0) ? abs(x) : -abs(x);

}

min16float3 reflectSpherical(min16float3 incident, min16float3 normal)
{
    return 2.0 * dot(incident, normal) * normal - incident;
    
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

float3 refractMitsuba(float3 wi, float3 n, float ior)
{
    if (ior == 1)
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

float3 refractZ(float3 f, float ior)
{
    return refractMitsuba(f, float3(0.0f, 0.0f, copy_sign(1.0f, f.z)), ior);

}

min16float tanTheta(min16float3 v)
{
    return max(0, sqrt(1.0 - v.z * v.z) / v.z);
}

float tanTheta(float3 v)
{
    return max(0, sqrt(1 - v.z * v.z) / v.z);
}

min16float sinTheta(min16float3 v)
{
    return max(0.0, sqrt(1.0 - v.z * v.z));
}

float sinTheta(float3 v)
{
    return max(0, sqrt(1.0f - v.z * v.z));
}

min16float sinTheta2(min16float3 v)
{
    return 1.0 - v.z * v.z;
}

float sinTheta2(float3 v)
{
    return 1.0f - v.z * v.z;
}

min16float cosTheta2(min16float3 v)
{
    return v.z * v.z;
}

float cosTheta2(float3 v)
{
    return v.z * v.z;
}

min16float Pow5(min16float x)
{
    min16float xSq = x * x;
    return xSq * xSq * x;
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
min16float fresnelDielectric(min16float incidentCosTheta, min16float ior)
{
    min16float transmittedCosTheta = 0.0;
    if (ior == 1.0)
    {
        return 0.0;
    }
    
    min16float scale = (incidentCosTheta > 0) ? 1.0 / ior : ior;
    min16float transmittedcosTheta2 = 1.0 - (1.0 - incidentCosTheta * incidentCosTheta) * (scale * scale);
    
    if (transmittedcosTheta2 <= 0.0)
    {
        return 1.0;
    }
    
    incidentCosTheta = abs(incidentCosTheta);
    transmittedCosTheta = sqrt(transmittedcosTheta2);
    
    min16float Rs = (incidentCosTheta - ior * transmittedCosTheta) / (incidentCosTheta + ior * transmittedCosTheta);
    min16float Rp = (ior * incidentCosTheta - transmittedCosTheta) / (ior * incidentCosTheta + transmittedCosTheta);
    
    return 0.5 * (Rs * Rs + Rp * Rp);
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
    
    real3 temp1 = ior * ior - kappa * kappa - sinTheta2.xxx;
    real3 a2pb2 = sqrt((temp1 * temp1 + kappa * kappa * ior * ior * 4.0.xxx));
    real3 a = sqrt(((a2pb2 + temp1) * 0.5));
    
    real3 term1 = a2pb2 + incidentCosTheta2.xxx;
    real3 term2 = a * (2.0 * incidentCosTheta);
    
    real3 Rs2 = (term1 - term2) / (term1 + term2);
    
    real3 term3 = a2pb2 * incidentCosTheta2 + sinTheta4.xxx;
    real3 term4 = term2 * sinTheta2;
    
    real3 Rp2 = Rs2 * (term3 - term4) / (term3 + term4);
    
    return 0.5 * (Rp2 + Rs2);
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
    pdf = (1.0f / PI) / (alphasquare * cosThetaM * cosThetaM * cosThetaM * temp * temp);

    float sinThetaM = sqrt(max(0.0f, 1 - cosThetaM * cosThetaM));
    
    return float3(sinThetaM * cosPhiM, sinThetaM * sinPhiM, cosThetaM);
}

real smithG1(real NdotV, real alpha)
{
    real a2 = alpha * alpha;
    return 2.0 / (1.0 + sqrt(1.0 + a2 * (1.0 - NdotV * NdotV)) / (NdotV * NdotV));
}

//assume isotropic
min16float smithG1(min16float3 v, min16float3 m, min16float alpha)
{
    //v.z == Frame::cosTheta
    if (dot(v, m) * v.z <= 0.0f)
    {
        return 0.0;
    }
    
    min16float tantheta = abs(tanTheta(v));
    if (tantheta <= 0.0) //TODO: almost equal op...
    {
        return 1.0;
    }
    
    
    min16float root = alpha * tantheta;
                    //hypot2
    return saturate(2.0 / (1.0 + 1.0 * 1.0 + root * root));

}

//assume isotropic
float smithG1(float3 v, float3 m, float alpha)
{
    alpha = clamp(alpha, EPSILON, 1.0);
    
    //v.z == Frame::cosTheta
    if (dot(v, m) * v.z <= 0.0f)
    {
        return 0.0f;
    }
    
    float tantheta = abs(tanTheta(v));
    if (tantheta <= 0.0f) //TODO: almost equal op...
    {
        return 1.0f;
    }
    
    
    float root = alpha * tantheta;
                    //hypot2
    return saturate(2.0 / (1.0 + 1.0 * 1.0 + root * root));

}

min16float D_GGX(min16float3 m, min16float alpha)
{
    if (m.z <= 0.0)
    {
        return 0.0;
    }
    alpha = max(HALFEPSILON, alpha);
    
    min16float costheta2 = cosTheta2(m);
    min16float beckmann = (((m.x * m.x) / (alpha * alpha) + (m.y * m.y) / (alpha * alpha)) / costheta2);
    
    min16float root = (1.0 + beckmann) * costheta2;
    return (1.0 / (HALF_PI * alpha * alpha * root * root));

}


float D_GGX(float3 m, float alpha)
{
   if (m.z <= 0.0f)
    {
        return 0.0f;
    }
    alpha = max(EPSILON, alpha);
    
    float costheta2 = cosTheta2(m);
    float beckmann = (((m.x * m.x) / (alpha * alpha) + (m.y * m.y) / (alpha * alpha)) / costheta2);
    
    float root = (1.0f + beckmann) * costheta2;
    return (1.0f / (PI * alpha * alpha * root * root));

}

float3 sample_GGX_Visible(float3 incident, float2 samplePoint, float rough, out float pdf)
{
    float3 incident_weighted = normalize(float3(rough * incident.x, rough * incident.y, incident.z));
    
    float theta = 0.0f;
    float phi = 0.0f;
    
    if (incident.z < 0.99999f)
    {
        theta = acos(incident_weighted.z);
        phi = atan2(incident_weighted.y, incident_weighted.x);
    }
    
    float sinphi;
    float cosphi;
    
    sincos(phi, sinphi, cosphi);
    float2 slope;
    //sample visible
    {
        static const float SQRT_PI_INV = 1.0f / sqrt(PI);
        
        if (theta < 1e-4f)
        {
            float sinphi = 0.0f;
            float cosphi = 0.0f;
            
            float r = sqrt(samplePoint.x / (1 - samplePoint.x)); //TODO: safe sqrt
            sincos(2 * PI * samplePoint.y, sinphi, cosphi);
            slope = float2(r * cosphi, r * sinphi);
        }
        
        float tantheta = tan(theta);
        float a = 1 / tantheta;
        float G1 = 2.0f / (1.0f + sqrt(1.0f + 1.0f / (a * a)));
        
        float A = 2.0f * samplePoint.x / G1 - 1.0f;
        if( abs(A) == 1.0f)
        {
            A -= sign(A) * EPSILON;
        }
        float tmp = 1.0f / (A * A - 1.0f);
        float B = tantheta;
        float D = sqrt(B * B * tmp * tmp - (A * A - B * B) * tmp);
        float slope_x_1 = B * tmp - D;
        float slope_X_2 = B * tmp + D;
        
        slope.x = (A < 0.0f || slope_X_2 > 1.0f / tantheta) ? slope_x_1 : slope_X_2;
        
        float S = 0.0f;
        if (samplePoint.y > 0.5f)
        {
            S = 1.0f;
            samplePoint.y = 2.0f * (samplePoint.y - 0.5f);
        }
        else
        {
            S = -1.0f;
            samplePoint.y = 2.0f * (0.5f - samplePoint.y);
        }
        
        float z = (samplePoint.y * (samplePoint.y * (samplePoint.y * (-0.365728915865723f) + 0.790235037209296f) -
                            0.424965825137544f) + 0.000152998850436920f) /
                        (samplePoint.y * (samplePoint.y * (samplePoint.y * (samplePoint.y * 0.169507819808272f - 0.397203533833404f) -
                            0.232500544458471f) + 1.0f) - 0.539825872510702f);

        slope.y = S * z * sqrt(1.0f + slope.x * slope.x);
    }
    
    //pdf
    {
        if (incident.z == 0.0f)
        {
            pdf = 0.0f;
        }
        else
        {
            static const float3 M = float3(0.0f, 0.0f, 1.0f);
            pdf = smithG1(incident, M, rough) * abs(dot(incident, M) * D_GGX(M, rough)) / abs(incident.z);
        }
    }
    
    slope = float2(cosphi * slope.x - sinphi * slope.y,
                    sinphi * slope.x + cosphi * slope.y);
    
    slope.x *= rough;
    slope.y *= rough;
    
    float norm = 1.0f / sqrt(slope.x * slope.x + slope.y * slope.y + 1.0f);

    return float3(-slope.x * norm, -slope.y * norm, norm);
}

//Karis 2013 distribution.
min16float D_GGX_Karis(min16float NdotH, min16float rough)
{
    min16float rough_square = rough * rough;
    min16float f = (NdotH * NdotH) * (rough_square - 1.0) + 1.0;
    return rough_square / (HALF_PI * f * f);
}

float D_GGX_Karis(float NdotH, float rough)
{
    float rough_square = rough * rough;
    float f = (NdotH * NdotH) * (rough_square - 1.0) + 1.0;
    return rough_square / (PI * f * f);
}

float smithG(float3 incident, float3 outgoing, float3 m, float rough)
{
    return smithG1(incident, m, rough) * smithG1(outgoing, m, rough);
}


//https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2017/Presentations/Hammon_Earl_PBR_Diffuse_Lighting.pdf
//Earl H. (2017)
min16float smithG2(min16float3 V, min16float3 N, min16float3 L, min16float rough)
{
    min16float numerator = 2.0 * (abs(dot(N, L)) * abs(dot(N, V)));
    min16float denominator = lerp(2.0 * abs(dot(N, L)) * abs(dot(N, V)), abs(dot(N, L)) + abs(dot(N, V)), rough);
    return numerator / denominator;
}

float smithG2(float3 V, float3 N, float3 L, float rough)
{
    float numerator = 2.0f * (abs(dot(N, L)) * abs(dot(N, V)));
    float denominator = lerp(2.0f * abs(dot(N, L)) * abs(dot(N, V)), abs(dot(N, L)) + abs(dot(N, V)), rough);
    return numerator / denominator;
}


// [Lagarde & De Rousiers, 2014] Moving Frostbite to Physically Based Rendering 3.0
//Lobe off-peak shift correction (Section 4.9.3)
float3 lobe_mean_shift(float3 incident, float masking_shadowing, float3 Normal)
{
    return masking_shadowing * Normal + (1.0 - masking_shadowing) * incident;
}

float3 specular_dominant(float3 normal, float3 outgoing, float NdotV, float rough)
{
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

real TIR_analytical(real cti, real rough, real ior_12, real ior_10)
{
    float3 wo = 0.xxx;
    wo.x = sin(0.25 * PI * (1.0 + cti));
    wo.z = sqrt(1.0 - (wo.x * wo.x));
    
    float alpha = max(rough, 0.01);
    float3 wi = 0.xxx;
    wi.z = cti;
    
    float3 normal = float3(0, 0, 1);
    
    float3 off_specular = specular_dominant(normal, wo, dot(normal, wi), alpha);
    
    
    float3 ws = -reflect(wo, off_specular);
    
    
    real value = 0;
   
    real R = fresnelDielectric(dot(wo, off_specular), ior_12);
        
    value += R * smithG1(ws.z, alpha);
        
    value *= (1.0 - fresnelDielectric(wo.z, ior_10)) * step(ws.z, 0); //branch elimination, make value 0 if ws.z <= 0 
             
   

    return value;
}


//Remaps 4D index to a 3D index, assuming 4th dimension is packed in Z.
float3 Dim4ToDim3(float4 coord, uint dimension)
{
    
    return float3(coord.x, coord.y, coord.z + (coord.w * dimension));
}

real3 sample_FGD(float cti, real alpha, real3 ior, real3 kappa)
{
    float3 output = 0.0f.xxx;
#if USE_FAST_FGD == 1
    //FGD 4D LUT parameterised by elevation, roughness and complex IOR
    //(Karis 2013) approximation parameterised by Roughness and angle
    //outputs scale & bias to F0.
    //TODO: work out 2D split sum approx.
    //https://cdn2.unrealengine.com/Resources/files/2013SiggraphPresentationsNotes-26915738.pdf
    //Ok split sum approx contains FGD_1 & FGD_2 (whatever that means)
    //tristimulus energy output?
    float2 splitsum = FGD_LUT.Sample(LUTSampler, float2(cti, alpha));
    
    //this crushes up the energy for dielectrics...
    //why doesn't this warn?
     output = (splitsum.x + float3_average(f0(ior)) * splitsum.y).xxx;
#else
    static const uint NumDimensions = 4;
    static const uint DimensionXY = 64;
    static const uint DimensionZ = DimensionXY * 32;
    static const uint DimensionW = 32;
    
    static const float DimensionXMin = 0.0f;
    static const float DimensionXMax = 1.0f;
    static const float DimensionYMin = 0.0f;
    static const float DimensionYMax = 1.0f;
    static const float DimensionZMin = 0.0f;
    static const float DimensionZMax = 4.0f;
    static const float DimensionWMin = 0.0f;
    static const float DimensionWMax = 4.0f;
    
   
    //remap ranges to 0-1 normalised index.
    const float cti_index = saturate((cti - DimensionXMin) / (DimensionXMax - DimensionXMin));
    const float alpha_index =  saturate((alpha - DimensionYMin) / (DimensionYMax - DimensionYMin));
    
    //remap Z & W index from being within [0,64] range to [0,2048]
    //Z index is normalised [0,64] range stretched to [0,2048]
    const float3 ior_norm = (ior - DimensionZMin) / (DimensionZMax - DimensionZMin);
    const float3 ior_index = saturate((trunc(ior_norm * (DimensionXY - 1)) / (DimensionXY - 1)));
    
    //W index is normalised [0,32] range + Z index.
    const float3 kappa_norm = (kappa - DimensionWMin) / (DimensionWMax - DimensionWMin);
    const float3 kappa_index = saturate((trunc(kappa_norm * (DimensionW - 1)) / (DimensionW - 1)) / (DimensionZ - 1));
   
    const float3 ZIndex = saturate(ior_index + kappa_index);
    
    //4D LUT is packed into 3D Texture, z & w share the same dimension. 
    //since it's normalised, need to reproduce the power.
    output.x = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.x));
    output.y = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.y));
    output.z = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.z));
#endif
    
    return output; // float3(cti_index * (DimensionXY - 1), alpha_index * (DimensionXY - 1), ZIndex.x * (DimensionZ - 1))
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
    
    //somehow this is 0. FGD LUT does have 1.0 values... Possible to hit them??
    t_ij = 1.0.xxx - r_ij; //can't be negative
    t_ji = t_ij;
    
}



void albedo(float cti, real alpha, real3 ior_ij, real3 kappa_ij, out real3 r_ij)
{
    //can't be negative! 
    
    r_ij = sample_FGD(cti, alpha, ior_ij, kappa_ij);

}

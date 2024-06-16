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


struct sample_record
{
    float3 incident;
    float3 outgoing;
    float ior;
    bool is_reflection_sample;
    bool sample_type;
};


//Lookup table for Total Internal Reflection
Texture3D<float> TIR_LUT : register(t18);

//compile time switch for (Karis 2013). split sum FGD approx, or the (Belcour 2018). FGD LUT.
#if USE_FAST_FGD == 1
//LUT for Karis Split-sum FGD approximation.
Texture2D<float2> FGD_LUT : register(t19);
#else
Texture3D<float> FGD_LUT : register(t20);
#endif


//IBL
TextureCube<float3> radianceIBLTexutre : register(t10);
TextureCube<float3> irradianceIBLTexture : register(t10);

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

float hg_lh_norm(float g)
{
    const bool g_neg = g < 0.f;
    g = abs(g);
    const float n = clamp(0.5039f - 0.8254f * g + 0.3226f * g * g, 0.f, 0.5f);
    return g_neg ? 1.0f - n : n;
}


bool isZero(float3 vec)
{
    return vec.x + vec.y + vec.z == 0;
}

float float3_max(float3 v)
{
    return max(v.x, max(v.y, v.z));
}


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


//TODO: Probably drop TIR for being too expensive.
float3 TIR_lookup(float3 coords)
{
    return TIR_LUT.Sample(LUTSampler, coords);
}


//Lagarde 2011. Compute fresnel reflectance at 0 degrees from IOR
float3 f0(float3 ior)
{
    return ((ior - 1.0f) * (ior - 1.0f)) / ((ior + 1.0f) * (ior + 1.0f));
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

float3 fresnelConductorExact(float incidentCosTheta, float3 ior, float3 kappa)
{
    float incidentCosTheta2 = incidentCosTheta * incidentCosTheta;
    float sinTheta2 = 1.0f - incidentCosTheta2;
    float sinTheta4 = sinTheta2 * sinTheta2;
    
    float3 temp1 = ior * ior - kappa * kappa - sinTheta2.xxx;
    float3 a2pb2 = sqrt((temp1 * temp1 + kappa * kappa * ior * ior * 4.0f.xxx));
    float3 a = sqrt(((a2pb2 + temp1) * 0.5f));
    
    float3 term1 = a2pb2 + incidentCosTheta2.xxx;
    float3 term2 = a * (2.0f * incidentCosTheta);
    
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
    pdf = (1.0f / PI) / (alphasquare * cosThetaM * cosThetaM * cosThetaM * temp * temp);

    float sinThetaM = sqrt(max(0.0f, 1 - cosThetaM * cosThetaM));
    
    return float3(sinThetaM * cosPhiM, sinThetaM * sinPhiM, cosThetaM);
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
    return saturate(2.0f / (1.0f + 1.0f * 1.0f + root * root));

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
float smithG2(float3 V, float3 N, float3 L, float rough)
{
    float numerator = 2.0f * (abs(dot(N, L)) * abs(dot(N, V)));
    float denominator = lerp(2.0f * abs(dot(N, L)) * abs(dot(N, V)), abs(dot(N, L)) + abs(dot(N, V)), rough);
    return numerator / denominator;
}


//Remaps 4D index to a 3D index, assuming 4th dimension is packed in Z.
float3 Dim4ToDim3(float4 coord, uint dimension)
{
    
    return float3(coord.x, coord.y, coord.z + (coord.w * dimension));
}

float3 sample_FGD(float cti, float alpha, float3 ior, float3 kappa)
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
     output = (splitsum.x + f0(ior).x * splitsum.y).xxx;
#else
    static const uint DimensionXY = 64;
    static const uint DimensionZ = DimensionXY * 32;
    
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
    const float3 ior_index = saturate((ior - DimensionZMin) / (DimensionZMax - DimensionZMin));
    //W index is normalised [0,64] range + Z index.
    const float3 kappa_index = saturate((kappa - DimensionWMin) / (DimensionWMax - DimensionWMin));
   
    const float3 ZIndex = saturate(ior_index + kappa_index);
    
    //4D LUT is packed into 3D Texture, z & w share the same dimension. 
    //since it's normalised, need to reproduce the power.
    output.x = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.x));
    output.y = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.y));
    output.z = FGD_LUT.Sample(LUTSampler, float3(cti_index, alpha_index, ZIndex.z));
#endif
    
    return output;
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
    
    r_ij = sample_FGD(cti, alpha, ior_ij.xxx, 0.0f.xxx);
    
    
    r_ji = r_ij;
    
    t_ij = 1.0f.xxx - r_ij; //can't be negative
    t_ji = t_ij;
    
}



void albedo(float cti, float alpha, float3 ior_ij, float3 kappa_ij, out float3 r_ij)
{
    //can't be negative! 
    
    r_ij = sample_FGD(cti, alpha, ior_ij, kappa_ij);

}
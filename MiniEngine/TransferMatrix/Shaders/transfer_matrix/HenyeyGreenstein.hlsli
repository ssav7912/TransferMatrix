
#pragma once

static const float PI = 3.14159265f;


///RGB variant of henyey greenstein?
struct henyey_greenstein
{
    float asymmetry;
    float3 mean;
    float3 norm;
};

///Henyey-Greenstein phase function.
///Equation 13.
///@param theta deviation angle (radians) from
///forward direction
///@param g asymmetry paramter, must be [-1,1].
float henyey_greenstein_phase_function(float theta, float g)
{
    g = clamp(g, -1.0f, 1.0f);
    
    return (1.0f / (4.0f * PI)) * ((1.0f - g * g) / pow((1.0f + g * g - 2.0f * g * cos(theta)), 3.0f / 2.0f));
    
}

///Equation 15
///Translates ggx roughness param to henyey greenstein asymmetry param.
float ggx_to_hg(float rough)
{
    rough = saturate(rough);
    return -0.085f + ((1.085f) / (1.0f + pow((rough / 0.5f), 1.3f)));

}

float hg_to_ggx(float asymmetry)
{
    asymmetry = clamp(asymmetry,-1.0f, 1.0f);
    return clamp(0.5f * pow((1.085f) / (max(asymmetry, 0.2f) + 0.085f) - 1.0f, 1.0f / 1.3f), 1e-4f, 1.f);

}

float hg_refract(float asymmetry, float ior)
{
    return min(sqrt(max(1.0f - (1.0f - asymmetry * asymmetry) * pow(ior, 0.75f), 0.0f)), 1.0f);
}

//equation 16
//@param g_i incoming light distribution with asymmetry g_i
//@param rough roughness of dielectric interface.
float dielectric_reflected_asymmetry(float g_i, float rough)
{
    return clamp(g_i, -1.0f, 1.0f) * ggx_to_hg(rough);
}

///Equation 9. 
///Belcour scaling factor.
///n = refractive_index = refractive index of what? 
float fake_rough_scaling(float refractive_index, float theta_i, float theta_t)
{
    return (1.0f / 2.0f) * (1.0f + refractive_index * (cos(theta_i) / cos(theta_t)));
}

///Equation 18.
///@param g_i asymmetry parameter of incident phase function.
///@param rough roughness of inicident interface.
///@param n_i refractive index of incident media.
///@param n_t refractive index of transmitted media.
float dielectric_transmited_asymmetry(float g_i, float rough, float n_i, float n_t, float theta_i, float theta_t, float ior)
{
    //equation 17.
    const float t = (1.0f - g_i * g_i) * pow((n_i / n_t), 3.0f / 4.0f);
    const float h = sqrt(1.0f - max(0, min(t, 1.0f)));
    
    //belcour fake scaling factor equation 9.

    const float s = fake_rough_scaling(ior, theta_i, theta_t);
    
    return h * ggx_to_hg(s * rough);
}


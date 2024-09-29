
#pragma once
#include "../Common.hlsli"

static const float PI = 3.14159265;
static const min16float HALF_PI = 3.14159265;


///RGB variant of henyey greenstein?
struct henyey_greenstein
{
    real asymmetry;
    real3 mean; //not even used?
    real3 norm;
};

struct hg_nomean
{
    real asymmetry;
    real3 norm;
};


henyey_greenstein zero_hg()
{
    henyey_greenstein h;
    h.asymmetry = 0.0;
    h.mean = 0.0.xxx;
    h.norm = 0.0.xxx;
    return h;
}

hg_nomean zero_hg_nomean()
{
    hg_nomean h;
    h.asymmetry = 0.0;
    h.norm = 0.0.xxx;
    return h;
}

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
min10float ggx_to_hg(min16float rough)
{
    rough = saturate(rough);
    return saturate(-0.085 + ((1.085) / (1.0 + pow((rough / 0.5), 1.3))));

}

real ggx_to_hg(float rough)
{
    rough = saturate(rough);
    return saturate(-0.085 + ((1.085) / (1.0 + pow((rough / 0.5), 1.3))));

}

real hg_to_ggx(real asymmetry)
{
    asymmetry = clamp(asymmetry,-1.0, 1.0);
    return clamp(0.5 * pow(1.085 / (max(asymmetry, 0.2) + 0.085) - 1.0, 1.0 / 1.3), 1e-4, 1.);

}

real hg_refract(real asymmetry, real ior)
{
    return min(sqrt(max(1.0 - (1.0 - asymmetry * asymmetry) * pow(ior, 0.75), 0.0)), 1.0);
}
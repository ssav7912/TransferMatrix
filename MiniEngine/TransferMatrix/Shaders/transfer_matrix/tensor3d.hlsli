#pragma once
#include "Constants.hlsli"
#include "../Common.hlsli"

struct tensor3d2x2
{

    /*
     [ _11 _12 ]
     [ _21 _22 ]*/
    
    real3 _11;
    real3 _12;
    real3 _21;
    real3 _22;
    
};


tensor3d2x2 mul(tensor3d2x2 a, tensor3d2x2 b)
{
    tensor3d2x2 x =
    {
        a._11 * b._11 + a._12 * b._21,
        a._11 * b._12 + a._12 * b._22,
        
        a._21 * b._11 + a._22 * b._21,
        a._21 * b._12 + a._22 * b._22
    };

    return x;

}

//sparse representation.
//will absolutely contribute to register pressure...
struct tensor3d6x6
{
    uint type;
    
    
    //row 1
    real3 _11;
    real3 _21;
    real3 _31;
    real3 _41;
    real3 _51;
    real3 _61;
    
    //row 2
    real3 _12;
    real3 _22;
    real3 _32;
    real3 _33;
    real3 _34;
    real3 _35;
    real3 _36;

    real3 _42;
    real3 _43;
    real3 _44;
    real3 _45;
    real3 _46;
    
    real3 _52;
    real3 _62;
    
 
};

tensor3d6x6 zero_tensor3d6x6()
{
    tensor3d6x6 x;
    x.type = TM_TYPE_NOCOMPONENT;
    x._11 = 0.0f.xxx;
    x._21 = 0.0f.xxx;
    x._31 = 0.0f.xxx;
    x._41 = 0.0f.xxx;
    x._51 = 0.0f.xxx;
    x._61 = 0.0f.xxx;
    x._12 = 0.0f.xxx;
    x._22 = 0.0f.xxx;
    x._32 = 0.0f.xxx;
    x._33 = 0.0f.xxx;
    x._34 = 0.0f.xxx;
    x._35 = 0.0f.xxx;
    x._36 = 0.0f.xxx;
    x._42 = 0.0f.xxx;
    x._43 = 0.0f.xxx;
    x._44 = 0.0f.xxx;
    x._45 = 0.0f.xxx;
    x._46 = 0.0f.xxx;
    x._52 = 0.0f.xxx;
    x._62 = 0.0f.xxx;
    return x;
}

tensor3d6x6 tensor3d6x6_identity()
{
    tensor3d6x6 x = zero_tensor3d6x6();
    x.type = TM_TYPE_NOCOMPONENT;
    x._11 = 1.0f.xxx;
    x._22 = 1.0f.xxx;
    x._33 = 1.0f.xxx;
    x._44 = 1.0f.xxx;
    return x;
}


tensor3d6x6 mul(tensor3d6x6 a, tensor3d6x6 b)
{
    tensor3d6x6 result = tensor3d6x6_identity();
    
    if (a.type == TM_TYPE_NOCOMPONENT)
    {
        return b;
    }
    
    if (b.type == TM_TYPE_NOCOMPONENT)
    {
        return a;
    }
    
    else if (b.type == TM_TYPE_DIELECTRICINTERFACE)
    {
        result._11 = a._11 * b._11 + a._12 * b._21;
        result._12 = a._11 * b._12 + a._12 * b._22;
        
        result._21 = a._21 * b._11 + a._22 * b._21;
        result._22 = a._21 * b._12 + a._22 * b._22;
        
        result._31 = a._31 * b._11 + a._32 * b._21;
        result._32 = a._31 * b._12 + a._32 * b._22;
        result._33 = a._33 * b._33 + a._34 * b._43;
        result._34 = a._33 * b._34 + a._34 * b._44;
        result._35 = a._35 * b._33 + a._36 * b._43;
        result._36 = a._35 * b._34 + a._36 * b._44;
        
        result._41 = a._41 * b._11 + a._42 * b._21;
        result._42 = a._41 * b._12 + a._42 * b._22;
        result._43 = a._43 * b._33 + a._44 * b._43;
        result._44 = a._43 * b._34 + a._44 * b._44;
        result._45 = a._45 * b._33 + a._46 * b._43;
        result._46 = a._45 * b._34 + a._46 * b._44;

        result._51 = a._51 * b._11 + a._52 * b._21;
        result._52 = a._51 * b._12 + a._52 * b._22;

        result._61 = a._61 * b._11 + a._62 * b._21;
        result._62 = a._61 * b._12 + a._62 * b._22;
        
        
    }
    else if (b.type == TM_TYPE_HOMOGENOUSMEDIUM)
    {
        result._11 = a._11 * b._11;
        result._12 = a._12 * b._22;

        result._21 = a._21 * b._11;
        result._22 = a._22 * b._22;

        result._31 = a._31 * b._11 + a._33 * b._31 - a._36 * b._36;
        result._32 = a._32 * b._22 + a._34 * b._42 + a._35 * b._36;
        result._33 = a._33 * b._33 - a._36 * b._36;
        result._34 = a._34 * b._44 + a._35 * b._36;
        result._35 = -a._34 * b._36 + a._35 * b._33;
        result._36 = a._33 * b._36 + a._36 * b._44;

        result._41 = a._41 * b._11 + a._43 * b._31 - a._46 * b._36;
        result._42 = a._42 * b._22 + a._44 * b._42 + a._45 * b._36;
        result._43 = a._43 * b._33 - a._46 * b._36;
        result._44 = a._44 * b._44 + a._45 * b._36;
        result._45 = -a._44 * b._36 + a._45 * b._33;
        result._46 = a._43 * b._36 + a._46 * b._44;

        result._51 = -a._34 * b._36 + a._35 * b._31 + a._51 * b._11;
        result._52 = a._33 * b._36 + a._36 * b._42 + a._52 * b._22;

        result._61 = -a._44 * b._36 + a._45 * b._31 + a._61 * b._11;
        result._62 = a._43 * b._36 + a._46 * b._42 + a._62 * b._22;
    }
    
    return result;
}
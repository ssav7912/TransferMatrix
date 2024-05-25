#pragma once


struct tensor3d2x2
{

    /*
     [ _11 _12 ]
     [ _21 _22 ]*/
    
    float3 _11;
    float3 _12;
    float3 _21;
    float3 _22;
    
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
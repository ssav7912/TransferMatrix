#pragma once
#include "tensor3d.hlsli"
#include "TransferMatrixCommon.hlsli"

//computes downward reflection energy for transfer matrix of interface separating 2 media.
float reflection_energy_tm2(float2x2 m_ij)
{
    return safe_div(m_ij._21, m_ij._11);
}

float3 reflection_energy_tm2(tensor3d2x2 m_ij)
{
    return safe_div(m_ij._21, m_ij._11);
}

float reflection_energy_tm2(float2x2 m_ij, float cond)
{
    return safe_div((m_ij._21 + m_ij._22 * cond), (m_ij._11 + m_ij._12 * cond));
}

float3 reflection_energy_tm2(tensor3d2x2 m, float3 cond)
{
    return safe_div((m._21 + m._22 * cond), (m._11 + m._12 * cond));
}

//computes downward transmission energy for transfer matrix of interface separating 2 media.
float transmission_energy_tm2(float2x2 m_ij)
{
    return safe_div(1.0f, m_ij._11);
}

float3 transmission_energy_tm2(tensor3d2x2 m_ij)
{
    return safe_div(1.0f.xxx, m_ij._11);
}

///computes dielectric reflected energy for 2 energy matrices.
///Special case for 2 interfaces.
float dielectric_reflected_energy_tm2(float2x2 energy_matrix_01, float2x2 energy_matrix_12)
{
    return reflection_energy_tm2(mul(energy_matrix_01, energy_matrix_12));

}

float3 dielectric_reflected_energy_tm2(tensor3d2x2 energy_matrix_01, tensor3d2x2 energy_matrix_12)
{
    return reflection_energy_tm2(mul(energy_matrix_01, energy_matrix_12));
}

///Equation 24.
///Explicit reflectance operator for conducting base.
///@param top_transfer: transfer matrix of upper layers.
///@param base_reflectance_factor: reflectance transfer factors of conducting base. 
float explicit_reflectance_tm2(float2x2 top_transfer, float base_reflectance_transfer)
{
    return (top_transfer._21 + top_transfer._22 * base_reflectance_transfer) / (top_transfer._11 + top_transfer._12 * base_reflectance_transfer);
}

float3 explicit_reflectance_tm2(tensor3d2x2 top_transfer, float base_reflectance_transfer)
{
    return (top_transfer._21 + top_transfer._22 * base_reflectance_transfer) / (top_transfer._11 + top_transfer._12 * base_reflectance_transfer);
}


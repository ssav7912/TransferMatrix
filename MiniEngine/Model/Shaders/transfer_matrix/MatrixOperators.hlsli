
//computes downward reflection energy for transfer matrix of interface separating 2 media.
float reflection_energy_tm2(float2x2 m_ij)
{
    return m_ij._21 / m_ij._11;
}

//computes downward transmission energy for transfer matrix of interface separating 2 media.
float transmission_energy_tm2(float2x2 m_ij)
{
    return 1.0f / m_ij._11;
}

///computes dielectric reflected energy for 2 energy matrices.
///Special case for 2 interfaces.
float dielectric_reflected_energy_tm2(float2x2 energy_matrix_01, float2x2 energy_matrix_12)
{
    return reflection_energy_tm2(mul(energy_matrix_01, energy_matrix_12));

}
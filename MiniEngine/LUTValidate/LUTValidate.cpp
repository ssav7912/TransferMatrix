// LUTValidate.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <format>

#include "DirectXMath.h"
#include "LUT.h"

struct float3
{
    float x;
    float y;
    float z;
};

/*
* Recreation of the sampling routine used in HLSL.
* 
* */
static float3 sample_flattened_FGD(float cti, float alpha, float3 ior, float3 kappa, const lut3& FGD)
{

    //something wrong with this!! 
    const auto range = FGD.range();
    const auto size = FGD.size();

    constexpr int32_t DimensionW = 64;
    constexpr float DimensionWMax = 4.0;
    constexpr float DimensionWMin = 0.0;

    const float cti_index = (cti - range[0].min) / (range[0].max - range[0].min);
    const float alpha_index = (alpha - range[1].min) / (range[1].max - range[1].min);

    const float3 ior_norm = {
        (ior.x - range[2].min) / (range[2].max - range[2].min),
        (ior.y - range[2].min) / (range[2].max - range[2].min),
        (ior.z - range[2].min) / (range[2].max - range[2].min),
    };

    const float3 ior_index = {
        (ior_norm.x * static_cast<float>((size[2]))),
        (ior_norm.y * static_cast<float>((size[2]))),
        (ior_norm.z * static_cast<float>((size[2]))),
    };

    const float3 kappa_norm = {
        (kappa.x - DimensionW) / (DimensionWMax - DimensionWMin),
        (kappa.y - DimensionW) / (DimensionWMax - DimensionWMin),
        (kappa.z - DimensionW) / (DimensionWMax - DimensionWMin),

    };

    const float3 kappa_index = {
        (kappa_norm.x * static_cast<float>((DimensionW))),
        (kappa_norm.y * static_cast<float>((DimensionW))),
        (kappa_norm.z * static_cast<float>((DimensionW))),

    };

    const float3 ZIndex = {
        (ior_index.x + kappa_index.x) / static_cast<float>(size[2]),
        (ior_index.y + kappa_index.y) / static_cast<float>(size[2]),
        (ior_index.z + kappa_index.z) / static_cast<float>(size[2]),

    };

    return float3{
        FGD.range_get_interpolate(cti_index, alpha_index, ZIndex.x),
        FGD.range_get_interpolate(cti_index, alpha_index, ZIndex.y),
        FGD.range_get_interpolate(cti_index, alpha_index, ZIndex.z),

    };
}
static float sample_flattened_FGD(float cti, float alpha, float ior, float kappa, const lut3& FGD)
{
    const auto range = FGD.range();
    const auto size = FGD.size();

    constexpr int32_t DimensionW = 64;
    constexpr float DimensionWMax = 4.0;
    constexpr float DimensionWMin = 0.0;

    const float cti_index = (cti - range[0].min) / (range[0].max - range[0].min);
    const float alpha_index = (alpha - range[1].min) / (range[1].max - range[1].min);

    const float ior_norm = (ior - range[2].min) / (range[2].max - range[2].min);


    const float ior_index = (ior_norm * static_cast<float>((size[2])));

    const float kappa_norm = (kappa - DimensionW) / (DimensionWMax - DimensionWMin);
      

    const float kappa_index = (kappa_norm * static_cast<float>((DimensionW )));

    const float ZIndex = (ior_index + kappa_index) / static_cast<float>(size[2]);

    return FGD.range_get_interpolate(cti_index, alpha_index, ZIndex);
      
}



int main()
{ 

    lut3 FGD_resample{ "tm_FGD_resample.bin" };
    lut4 FGD{ "tm_FGD.bin" };
    lut4 FGD_validate{ "tm_FGD_validate.bin" };


    std::cout << "Validation Spot Tests" << std::endl;


    const float3 bottomlayerIOR = { 1.0, 1.0, 1.0 };
    const float3 bottomlayerKappa = { 1.0, 0.1, 0.1 };
    const float bottomlayerAlpha = 0.01;
    const float cti = 1;
    const auto fgd_ref_x = FGD.range_get_interpolate(cti, bottomlayerAlpha, bottomlayerIOR.x, bottomlayerKappa.x);
    const auto fgd_ref_y = FGD.range_get_interpolate(cti, bottomlayerAlpha, bottomlayerIOR.x, bottomlayerKappa.y);
    const auto fgd_ref_z = FGD.range_get_interpolate(cti, bottomlayerAlpha, bottomlayerIOR.x, bottomlayerKappa.z);

    const auto fgd_flat = sample_flattened_FGD(cti, bottomlayerAlpha, bottomlayerIOR, bottomlayerKappa, FGD_resample);

    std::cout << std::format("Reference FGD for cti = {0}, alpha = {1}, IOR = {2}, Kappa = [{3},{4},{5}] is [{6},{7},{8}]", cti, bottomlayerAlpha, bottomlayerIOR.x, 
        bottomlayerKappa.x, bottomlayerKappa.y, bottomlayerKappa.z,
        fgd_ref_x, fgd_ref_y, fgd_ref_z) << std::endl;

    std::cout << std::format("Flattened FGD for cti = {0}, alpha = {1}, IOR = {2}, Kappa = [{3},{4},{5}] is [{6},{7},{8}]", cti, bottomlayerAlpha, bottomlayerIOR.x,
        bottomlayerKappa.x, bottomlayerKappa.y, bottomlayerKappa.z,
        fgd_flat.x, fgd_flat.y, fgd_flat.z) << std::endl;


    const auto fgd_ref_zeroIOR_x = FGD.range_get_interpolate(cti, bottomlayerAlpha, 0, bottomlayerKappa.x);
    const auto fgd_ref_zeroIOR_y = FGD.range_get_interpolate(cti, bottomlayerAlpha, 0, bottomlayerKappa.y);
    const auto fgd_ref_zeroIOR_z = FGD.range_get_interpolate(cti, bottomlayerAlpha, 0, bottomlayerKappa.z);

    const auto fgd_flat_zeroIOR = sample_flattened_FGD(cti, bottomlayerAlpha, {0}, bottomlayerKappa, FGD_resample);

    std::cout << std::format("Reference FGD for cti = {0}, alpha = {1}, IOR = {2}, Kappa = [{3},{4},{5}] is [{6},{7},{8}]",
        cti, bottomlayerAlpha, 0, bottomlayerKappa.x, bottomlayerKappa.y, bottomlayerKappa.z, fgd_ref_zeroIOR_x, fgd_ref_zeroIOR_y, fgd_ref_zeroIOR_z) 
        << std::endl;


    std::cout << std::format("Flattened FGD for cti = {0}, alpha = {1}, IOR = {2}, Kappa = [{3},{4},{5}] is [{6},{7},{8}]",
        cti, bottomlayerAlpha, 0, bottomlayerKappa.x, bottomlayerKappa.y, bottomlayerKappa.z, fgd_flat_zeroIOR.x,fgd_flat_zeroIOR.y, fgd_flat_zeroIOR.z)
        << std::endl;

    const float3 goldIOR = { 0.1, 0.4, 1.4 };
    const float3 goldKappa = { 4.0, 2.4, 1.6 };
    const float goldAlpha = 0.03;
    
    const auto fgd_ref_gold_x = FGD.range_get_interpolate(cti, goldAlpha, goldIOR.x, goldKappa.x);
    const auto fgd_ref_gold_y = FGD.range_get_interpolate(cti, goldAlpha, goldIOR.y, goldKappa.y);
    const auto fgd_ref_gold_z = FGD.range_get_interpolate(cti, goldAlpha, goldIOR.z, goldKappa.z);

    std::cout << std::format("Reference FGD for cti = {0}, alpha = {1}, IOR = [{2}, {3}, {4}], Kappa = [{5},{6},{7}] is [{8},{9},{10}]",
        cti, goldAlpha,goldIOR.x, goldIOR.y, goldIOR.z,
        goldKappa.x, goldKappa.y, goldKappa.z, 
        fgd_ref_gold_x, fgd_ref_gold_y, fgd_ref_gold_z)
        << std::endl;


    uint64_t num_validate_errors = 0;
    uint64_t num_errors = 0;
    uint64_t num_samples = 0;

    float validate_max_diff = -std::numeric_limits<float>::infinity();
    float validate_min_diff = std::numeric_limits<float>::infinity();
    float max_diff = -std::numeric_limits<float>::infinity();
    float min_diff = std::numeric_limits<float>::infinity();
    float avg_error = 0.0;
    float avg_validate_error = 0.0;
    uint64_t border_error = 0;

    std::cout << "Validating entire LUT...." << std::endl;
    for (int i = 0; i < 64; i++)
    {
        for (int j = 0; j < 64; j++)
        {
            for (int k = 0; k < 64; k++)
            {
                for (int l = 0; l < 64; l++)
                {
                    const float cti = static_cast<float>(i) / 64.0;
                    const float alpha = static_cast<float>(j) / 64.0;
                    const float ior = (static_cast<float>(k) / 64.0) * 4.0;
                    const float kappa = (static_cast<float>(l) / 64.0) * 4.0;

                    //const float kappa_flat = static_cast<float>(l) / 32.0;


                    const auto fgd_sample = FGD.range_get_interpolate(cti, alpha, ior, kappa);
                    const auto fgd_validate_sample = FGD_validate.range_get_interpolate(cti, alpha, ior, kappa);
                    const auto fgd_flat_sample = sample_flattened_FGD(cti, alpha, ior, kappa, FGD_resample);

                    const float validate_epsilon = 0.001;
                    const float epsilon = 0.2;
                    const float diff = std::abs(fgd_sample - fgd_flat_sample);
                    const float validate_diff = std::abs(fgd_sample - fgd_validate_sample);

                    validate_max_diff = std::max(validate_max_diff, validate_diff);
                    validate_min_diff = std::min(validate_min_diff, validate_diff);
                    avg_validate_error += validate_diff;

                    max_diff = std::max(diff, max_diff);
                    min_diff = std::min(diff, min_diff);
                    avg_error += diff;

                    if (std::abs(fgd_sample - fgd_flat_sample) > epsilon)
                    {
                        if (i % 64 < 3 || j % 64 < 3 || k % 64 < 3 || l % 64 < 3)
                        {
                            border_error++;
                        }
                        // std::cout << std::format("Reference FGD and Flattened FGD differ by more than {0} (delta = {1}), \n" 
                            //"   for cti = {2}, alpha = {3}, IOR = {4} and Kappa = {5}",
                           // epsilon, diff, cti, alpha, ior, kappa) << std::endl;

                        num_errors++;
                    }

                    if (validate_diff > validate_epsilon)
                    {
                        num_validate_errors++;
                    }

                    num_samples++;


                }

            }

        }
    }
    std::cout << std::format("FGD vs Resampled: {0} Errors out of {1} Samples. {2} At axis borders. {3} Rate. Min Delta: {4}, Max: {5}, Avg: {6}", 
        num_errors, num_samples, border_error, static_cast<float>(num_errors) / static_cast<float>(num_samples), 
        min_diff, max_diff, avg_error / static_cast<float>(num_samples)) << std::endl;
    std::cout << "Complete." << std::endl;

    const float validate_error_rate = static_cast<float>(num_validate_errors) / static_cast<float>(num_samples);
    avg_validate_error = avg_validate_error / static_cast<float>(num_samples);

    std::cout << std::format("FGD vs FGD Validation: {0} Errors out of {1} Samples. \n" \
        "Error Rate: {2}\n" \
        "Min Delta: {3}, Max: {4}, Avg: {5}.", \
        num_validate_errors, num_samples, validate_error_rate, validate_min_diff, validate_max_diff, avg_validate_error) << std::endl;

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

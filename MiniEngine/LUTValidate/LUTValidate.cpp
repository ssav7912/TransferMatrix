// LUTValidate.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <format>

#include "DirectXMath.h"
#include "LUT.h"




int main()
{ 


    lut4 FGD{ "tm_FGD.bin" };

    std::cout << "Validation Spot Tests" << std::endl;


    const float bottomlayerIOR = 1.0;
    const float bottomlayerKappa[3] = { 1.0, 0.1, 0.1 };
    const float bottomlayerAlpha = 0.1;
    const float cti = 0.77;
    const auto fgd_ref_x = FGD.range_get_interpolate(cti, bottomlayerAlpha, bottomlayerIOR, bottomlayerKappa[0]);
    const auto fgd_ref_y = FGD.range_get_interpolate(cti, bottomlayerAlpha, bottomlayerIOR, bottomlayerKappa[1]);
    const auto fgd_ref_z = FGD.range_get_interpolate(cti, bottomlayerAlpha, bottomlayerIOR, bottomlayerKappa[2]);

    std::cout << std::format("Reference FGD for cti = {0}, alpha = {1}, IOR = {2}, Kappa = [{3},{4},{5}] is [{6},{7},{8}]", cti, bottomlayerAlpha, bottomlayerIOR, 
        bottomlayerKappa[0], bottomlayerKappa[1], bottomlayerKappa[2],
        fgd_ref_x, fgd_ref_y, fgd_ref_z) << std::endl;


    


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

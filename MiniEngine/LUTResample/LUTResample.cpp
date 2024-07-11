// LUTResample.cpp : Resamples LUTs used in the transfer matrix shader to half precision.
//

#include <iostream>
#include "../TransferMatrix/TransferMatrixResources.h"
#include "../TransferMatrix/Texture3D.h"
#include "DirectXMath.h"
#include "DirectXPackedVector.h"
#include "DirectXTex.h"

int main()
{
    std::cout << "Writing TIR LUT" << std::endl;
    //TIR LUT
    using DirectX::PackedVector::HALF;

    constexpr size_t TIR_DIM = 64;
    constexpr size_t TIR_DEPTH = 3;
    constexpr size_t TIR_ELEMENTS = TIR_DIM * TIR_DIM * TIR_DIM;
    std::vector<float> TIR_Data (static_cast<size_t>(std::pow(TIR_DIM, TIR_DEPTH)) );
    //std::vector<HALF> Resampled_TIR_Data(sizeof(HALF) * static_cast<size_t>(std::pow(TIR_DIM, TIR_DEPTH)));
    TransferMatrixResources::LoadLUTFromFile<TIR_DEPTH>("tm_TIR.bin", TIR_Data);

    DirectX::ScratchImage TIR_image{};
    TIR_image.Initialize3D(DXGI_FORMAT_R16_FLOAT, TIR_DIM, TIR_DIM, TIR_DIM, 0);

    auto TIR_begin = reinterpret_cast<HALF*>(TIR_image.GetPixels());

    for (size_t i = 0; i < TIR_ELEMENTS; i++)
    {
        auto sample = TIR_Data[i];
        TIR_begin[i] = DirectX::PackedVector::XMConvertFloatToHalf(sample);
    }

    DirectX::SaveToDDSFile(TIR_image.GetImages(), TIR_image.GetImageCount(), TIR_image.GetMetadata(), DDS_FLAGS_NONE, L"tm_TIR.dds");

    std::cout << "Writing FGD LUT" << std::endl;

    //FGD LUT
    constexpr size_t FGD_DIM = 64;
    constexpr size_t FGD_DEPTH = 4;
    constexpr size_t FGD_ELEMENTS = FGD_DIM * FGD_DIM * FGD_DIM * FGD_DIM;

    std::vector<float> FGD_Data(FGD_ELEMENTS);
    TransferMatrixResources::LoadLUTFromFile<FGD_DEPTH>("tm_FGD.bin", FGD_Data);


    DirectX::ScratchImage FGD_image{};
    FGD_image.Initialize3D(DXGI_FORMAT_R16_FLOAT, FGD_DIM, FGD_DIM, FGD_DIM * (FGD_DIM / 2), 0);

    auto FGD_begin = reinterpret_cast<HALF*>(FGD_image.GetPixels());
    const size_t fgd_image_data_size = FGD_DIM * FGD_DIM * FGD_DIM * (FGD_DIM / 2 );
    constexpr size_t StrideDimension4 = FGD_DIM * FGD_DIM * FGD_DIM;

    //crush 4th dimension to 32. 
    size_t i = 0;
    while (i < fgd_image_data_size)
    {
        const auto src_index = i;
        const auto end_index = i + StrideDimension4;
        const auto dest_index = i;

        //copies per-slice, transforming each element into half precision.
        std::transform(FGD_Data.begin() + src_index, FGD_Data.begin() + std::min(fgd_image_data_size, end_index), FGD_begin + dest_index, DirectX::PackedVector::XMConvertFloatToHalf);
        if (i % 2 == 0 && i != 0)
        {
            i++;
        }
        else {
            i += StrideDimension4;
        }
    }

    DirectX::SaveToDDSFile(FGD_image.GetImages(), FGD_image.GetImageCount(), FGD_image.GetMetadata(), DDS_FLAGS_NONE, L"tm_FGD.dds");


}

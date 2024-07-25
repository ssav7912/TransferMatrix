#pragma once
#include "../Core/PipelineState.h"
#include "../Core/DescriptorHeap.h"
#include "../Core/Texture.h"
#include "Texture3D.h"
#include "../Core/TextureManager.h"

#include <iostream>
#include <filesystem>
#include <fstream>


class TransferMatrixResources
{
public:
	TransferMatrixResources() = default;
	TransferMatrixResources(const std::string & FGD_path, const std::wstring& FGD_4D_path, const std::wstring& GD_path, const std::wstring & TIR_path);

	void Initialise(const std::string& FGD_path, const std::wstring& FGD_4D_path, const std::wstring& GD_path, const std::wstring& TIR_path);

	bool UseTM6 = false;

	GraphicsPSO TM2PSO{ L"2-flux Transfer Matrix PSO" };
	GraphicsPSO TM6PSO{ L"6-flux Transfer Matrix PSO" };

	Texture3D TIR_LUT;
	TextureRef FGD_LUT;
	Texture GD_LUT;
	Texture3D FGD_4D_LUT;


	static constexpr int32_t NUM_LAYERS = 2;
	static constexpr int32_t MAX_LAYERS = 5;
	static constexpr int32_t NUM_ATTRIBUTES = 7;
	static constexpr int32_t MAX_TEXTURES = MAX_LAYERS * NUM_ATTRIBUTES;

	

	Texture3D LoadTIRLutFromFile(const std::wstring& TIR_path);
	Texture3D LoadFGDLUTFromFile(const std::wstring& FGD_path);
	Texture LoadGDLUTFromFile(const std::wstring& GD_path);

	template<size_t LutDimension>
	static void LoadLUTFromFile(const std::wstring& path, std::vector<float>& data)
	{
		struct range {
			float min;
			float max;
		};

		std::ifstream in{ path, std::ios_base::in | std::ios_base::binary };

		

		int32_t size[LutDimension] = { 0 };
		range lut_range[LutDimension] = { 0 };

		for (int i = 0; i < LutDimension; i++)
		{
			in.read(reinterpret_cast<char*>(&size[i]), sizeof(int32_t));
		}

		for (int i = 0; i < LutDimension; i++)
		{
			in.read(reinterpret_cast<char*>(&lut_range[i].min), sizeof(float));
			in.read(reinterpret_cast<char*>(&lut_range[i].max), sizeof(float));
		}

		size_t linear_size = [size]()->size_t {int lsize = 1;
		for (int i = 0; i < LutDimension; i++)
		{
			lsize *= size[i];
		}
		return lsize;  }();

		data.assign(linear_size, float());
		in.read(reinterpret_cast<char*>(data.data()), linear_size * sizeof(float));
	}
};

//float3 padded out to 16 bytes for HLSL cbuffer array alignment.
struct float3_padded
{
	float xyz[3];
	float pad;

};

//float padded out to 16 bytes for HLSL cbuffer array alignment.
__declspec(align(16)) struct float_with_pad
{
	float x;
};

//Constant buffer for layer parameters
struct LayerConstants
{
	__declspec(align(16)) float3_padded IORs[TransferMatrixResources::MAX_LAYERS];
	__declspec(align(16)) float3_padded Kappas[TransferMatrixResources::MAX_LAYERS];

	__declspec(align(16)) float3_padded Sigma_S[TransferMatrixResources::MAX_LAYERS];

	__declspec(align(16)) float3_padded Sigma_K[TransferMatrixResources::MAX_LAYERS];

	__declspec(align(16)) float_with_pad Depths[TransferMatrixResources::MAX_LAYERS];

	__declspec(align(16)) float_with_pad G[TransferMatrixResources::MAX_LAYERS];

	__declspec(align(16)) float_with_pad Roughs[TransferMatrixResources::MAX_LAYERS];
	//uint8_t _pad3[16 * 3];
	int32_t layers;
	int32_t num_samples;
};


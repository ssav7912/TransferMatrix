#pragma once
#include "../Core/PipelineState.h"
#include "../Core/DescriptorHeap.h"
#include "../Core/Texture.h"
#include "Texture3D.h"
#include "../Core/TextureManager.h"

#include <iostream>
#include <filesystem>
#include <fstream>

class TM2Resources
{
public:
	TM2Resources() = default;
	TM2Resources(const std::string & FGD_path, const std::string& FGD_4D_path, const std::string & TIR_path);

	void Initialise(const std::string& FGD_path, const std::string& FGD_4D_path, const std::string& TIR_path);


	GraphicsPSO TM2PSO{ L"2-flux Transfer Matrix PSO" };

	Texture3D TIR_LUT;
	TextureRef FGD_LUT;
	Texture3D FGD_4D_LUT;

	static constexpr int32_t NUM_LAYERS = 2;
	static constexpr int32_t MAX_LAYERS = 5;
	
	//const void* TM2_PS = g_pTM2DielectricPS;
	//const void* TM2_VS = g_pDefaultVS;

	//static constexpr size_t sizeof_TM2PS = sizeof(g_pTM2DielectricPS);
	//static constexpr size_t sizeof_TM2VS = sizeof(g_pDefaultVS);

private:

	Texture3D LoadTIRLutFromFile(const std::string& TIR_path);
	Texture3D LoadFGDLUTFromFile(const std::string& FGD_path);

	template<size_t LutDimension>
	void LoadLUTFromFile(const std::string& path, std::vector<float>& data)
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
	__declspec(align(16)) float3_padded IORs[TM2Resources::MAX_LAYERS];
	__declspec(align(16)) float3_padded Kappas[TM2Resources::MAX_LAYERS];

	__declspec(align(16)) float3_padded Sigma_S[TM2Resources::MAX_LAYERS];

	__declspec(align(16)) float3_padded Sigma_K[TM2Resources::MAX_LAYERS];

	__declspec(align(16)) float_with_pad Depths[TM2Resources::MAX_LAYERS];

	__declspec(align(16)) float_with_pad G[TM2Resources::MAX_LAYERS];

	__declspec(align(16)) float_with_pad Roughs[TM2Resources::MAX_LAYERS];
	//uint8_t _pad3[16 * 3];
	int32_t layers;
	int32_t num_samples;
};


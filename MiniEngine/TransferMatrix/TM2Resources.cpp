#include "TM2Resources.h"

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include "Texture3D.h"

//TODO: fix this stupid shit
#include "../Build/x64/Debug/Output/TransferMatrix/CompiledShaders/DefaultVS.h"
#include "../Build/x64/Debug/Output/TransferMatrix/CompiledShaders/TM2DielectricPS.h"
#include "../Build/x64/Debug/Output/TransferMatrix/CompiledShaders/TM2OpaquePS.h"


TM2Resources::TM2Resources(const std::string & FGD_path, const std::string& FGD_4D_path, const std::string & TIR_path)
{

	//init PSO
	Initialise(FGD_path, FGD_4D_path, TIR_path);

}

void TM2Resources::Initialise(const std::string& FGD_path, const std::string& FGD_4D_path, const std::string& TIR_path)
{

	auto tex = TextureManager::LoadDDSFromFile(FGD_path);
	FGD_LUT = tex;

	//init PSO
	TM2PSO.SetRasterizerState(Graphics::RasterizerDefault);
	TM2PSO.SetBlendState(Graphics::BlendDisable);
	TM2PSO.SetDepthStencilState(Graphics::DepthStateReadWrite);
	TM2PSO.SetInputLayout(0, nullptr);
	TM2PSO.SetPrimitiveTopologyType(D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE);
	TM2PSO.SetVertexShader(g_pDefaultVS, sizeof(g_pDefaultVS));
	TM2PSO.SetPixelShader(g_pTM2OpaquePS, sizeof(g_pTM2OpaquePS));


	TIR_LUT = LoadTIRLutFromFile(TIR_path);
	FGD_4D_LUT = LoadFGDLUTFromFile(FGD_4D_path);

}

/// <summary>
/// Read TIR bin LUT from file and create 3D TIR texture.
/// Adapted from 
/// Randrianandrasana et al. (2021) Transfer Matrix Based Layered Materials Rendering. 
/// </summary>
/// <param name="TIR_path"></param>
/// <returns></returns>
Texture3D TM2Resources::LoadTIRLutFromFile(const std::string& TIR_path)
{
	struct range {
		float min;
		float max;
	};

	std::ifstream in {TIR_path, std::ios_base::in | std::ios_base::binary};

	constexpr int32_t LUT_DIM = 3;

	std::vector<float> data;
	int32_t size[LUT_DIM] = { 0 };
	range lut_range[LUT_DIM] = { 0 };

	for (int i = 0; i < LUT_DIM; i++)
	{
		in.read(reinterpret_cast<char*>(&size[i]), sizeof(int32_t));
	}

	for (int i = 0; i < LUT_DIM; i++)
	{
		in.read(reinterpret_cast<char*>(&lut_range[i].min), sizeof(float));
		in.read(reinterpret_cast<char*>(&lut_range[i].max), sizeof(float)); 
	}

	size_t linear_size = [LUT_DIM, size]()->size_t {int lsize = 1;
	for (int i = 0; i < LUT_DIM; i++)
	{
		lsize *= size[i];
	}
	return lsize;  }();

	data.assign(linear_size, float());
	in.read(reinterpret_cast<char*>(data.data()), linear_size * sizeof(float));



	Texture3D TIR{};
	constexpr size_t Dim = 64;


	ASSERT(data.size() == Dim * Dim * Dim);
	
	const size_t RowPitchBytes = Dim * sizeof(float);
	TIR.Create3D(RowPitchBytes, Dim, Dim, Dim, DXGI_FORMAT_R32_FLOAT, data.data());
	return TIR;
}

Texture3D TM2Resources::LoadFGDLUTFromFile(const std::string& FGD_path)
{
	std::vector<float> data;
	LoadLUTFromFile<4>(FGD_path, data);

	Texture3D FGD{};
	constexpr size_t Dim = 64;

	//crush 4th dimension by dropping every 2nd value to fit into DX12 resource limits (HACK!!!!)
	//for (int i = 0; i < Dim; i++)
	//{
	//	for (int j = 0; j < Dim; j++)
	//	{
	//		for (int k = 0; k < Dim; k++)
	//		{
	//			for (int l = 0; l < Dim; l++)
	//			{
	//				size_t linear_index = i * (Dim * Dim * Dim) + j * (Dim * Dim) + k * Dim + l;

	//				if (l % 2 == 0)
	//				{
	//					data.erase(data.begin() + linear_index);
	//				}


	//			}
	//		}
	//	}
	//}

	//crushed vector
	std::vector<float> new_data(Dim * Dim * Dim * (Dim / 2), std::nanf("NaN"));

	const size_t data_size = data.size();
	constexpr size_t StrideDimension4 = Dim * Dim * Dim; //access every 4th dimensional element. 


	size_t i = 0; 
	while (i < new_data.size())
	{
		//copy from i to next element in stride
		const auto src_index = i;
		const auto end_index = i + StrideDimension4;
		const auto dest_index =  i;
		std::copy(data.begin() + src_index, data.begin() + std::min(new_data.size(), end_index), new_data.begin() + dest_index);
		if (i % 2 == 0 && i != 0)
		{
			i++;
		}
		else {
			i += StrideDimension4;
		}
	}

	ASSERT(std::none_of(new_data.cbegin(), new_data.cend(), [](float v) -> bool {return std::isnan(v); }));

	//4D texture
	ASSERT(new_data.size() == Dim * Dim * Dim * (Dim/2));

	//encode 4D LUT in 3D texture, Z dim encodes Z + W of index.
	//i.e. dim = 64 * 64 * (64 * 64) = 64 * 64 * 4096. 
	const size_t RowPitchbytes = Dim * sizeof(float);

	FGD.Create3D(RowPitchbytes, Dim, Dim, Dim * (Dim/2), DXGI_FORMAT_R32_FLOAT, new_data.data());

	return FGD;
}

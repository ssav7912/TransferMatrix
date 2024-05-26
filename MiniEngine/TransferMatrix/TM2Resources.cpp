#include "TM2Resources.h"

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include "Texture3D.h"

//TODO: fix this stupid shit
#include "../Build/x64/Debug/Output/TransferMatrix/CompiledShaders/DefaultVS.h"
#include "../Build/x64/Debug/Output/TransferMatrix/CompiledShaders/TM2DielectricPS.h"

TM2Resources::TM2Resources(const std::string & FGD_path, const std::string & TIR_path) : FGD_LUT(TextureManager::LoadDDSFromFile(FGD_path))
{

	//init PSO
	TM2PSO.SetRasterizerState(Graphics::RasterizerDefault);
	TM2PSO.SetBlendState(Graphics::BlendDisable);
	TM2PSO.SetDepthStencilState(Graphics::DepthStateReadWrite);
	TM2PSO.SetInputLayout(0, nullptr);
	TM2PSO.SetPrimitiveTopologyType(D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE);
	TM2PSO.SetVertexShader(g_pDefaultVS, sizeof(g_pDefaultVS));
	TM2PSO.SetPixelShader(g_pTM2DielectricPS, sizeof(g_pTM2DielectricPS));

	
	TIR_LUT = LoadTIRLutFromFile(TIR_path);

}

void TM2Resources::Initialise(const std::string& FGD_path, const std::string& TIR_path)
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
	TM2PSO.SetPixelShader(g_pTM2DielectricPS, sizeof(g_pTM2DielectricPS));


	TIR_LUT = LoadTIRLutFromFile(TIR_path);

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

	const size_t RowPitchBytes = Dim * sizeof(float) * LUT_DIM;
	TIR.Create3D(RowPitchBytes, Dim, Dim, Dim, DXGI_FORMAT_R32G32B32_FLOAT, data.data());
	return TIR;
}

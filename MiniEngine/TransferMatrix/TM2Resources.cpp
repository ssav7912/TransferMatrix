#include "TM2Resources.h"

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include "CompiledShaders/DefaultVS.h"
#include "CompiledSahders/TM2DielectricPS.h"

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

	

}

TextureRef TM2Resources::LoadTIRLutFromFile(const std::string& TIR_path)
{
	struct range {
		float min;
		float max;
	};

	std::ifstream in {TIR_path, std::ios_base::in | std::ios_base::binary};

	constexpr int32_t LUT_DIM = 3;

	std::vector<float> data;
	int32_t size[LUT_DIM];
	range lut_range[LUT_DIM];

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


	Texture t{};
	t.CreateDDSFromMemory()
	return TextureRef();
}

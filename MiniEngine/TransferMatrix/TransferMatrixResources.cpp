#include "TransferMatrixResources.h"

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <codecvt>
#include <locale>

#include "Texture3D.h"
#include "DirectXTex.h"

#include "DefaultVS.h"
//#include "TM2DielectricPS.h"
#include "TM2OpaquePS.h"
#include "TM6PS.h"

TransferMatrixResources::TransferMatrixResources(const std::string & FGD_path, const std::wstring& FGD_4D_path, const std::wstring& GD_path, const std::wstring & TIR_path)
{

	//init PSO
	Initialise(FGD_path, FGD_4D_path, GD_path, TIR_path);

}

void TransferMatrixResources::Initialise(const std::string& FGD_path, const std::wstring& FGD_4D_path, const std::wstring& GD_path, const std::wstring& TIR_path)
{

	auto tex = TextureManager::LoadDDSFromFile(FGD_path);
	FGD_LUT = tex;

	GD_LUT = LoadGDLUTFromFile(GD_path);


	//init PSO
	TM2PSO.SetRasterizerState(Graphics::RasterizerDefault);
	TM2PSO.SetBlendState(Graphics::BlendDisable);
	TM2PSO.SetDepthStencilState(Graphics::DepthStateReadWrite);
	TM2PSO.SetInputLayout(0, nullptr);
	TM2PSO.SetPrimitiveTopologyType(D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE);
	TM2PSO.SetVertexShader(g_pDefaultVS, sizeof(g_pDefaultVS));
	TM2PSO.SetPixelShader(g_pTM2OpaquePS, sizeof(g_pTM2OpaquePS));

	TM6PSO = TM2PSO;
	TM6PSO.SetPixelShader(g_pTM6PS, sizeof(g_pTM6PS));


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
Texture3D TransferMatrixResources::LoadTIRLutFromFile(const std::wstring& TIR_path)
{


	Texture3D TIR{};
	constexpr size_t Dim = 64;

#if USEFP16 == 1
	const DXGI_FORMAT pixel_format = DXGI_FORMAT_R16_FLOAT;
	TexMetadata md;
	
	DirectX::ScratchImage out_image{};



	HRESULT res = DirectX::LoadFromDDSFile(TIR_path.c_str(), DDS_FLAGS_NONE, &md, out_image);

	ASSERT_SUCCEEDED(res);
	
	const void* data_ptr = out_image.GetPixels();
	const size_t RowPitchBytes = Dim * sizeof(uint16_t); //half precision

#else
	const DXGI_FORMAT pixel_format = DXGI_FORMAT_R32_FLOAT;

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



	ASSERT(data.size() == Dim * Dim * Dim);

	const void* data_ptr = data.data();


	const size_t RowPitchBytes = Dim * sizeof(float);
#endif

	
	TIR.Create3D(RowPitchBytes, Dim, Dim, Dim, pixel_format, data_ptr);

	return TIR;
}

Texture3D TransferMatrixResources::LoadFGDLUTFromFile(const std::wstring& FGD_path)
{

	Texture3D FGD{};
	constexpr size_t Dim = 64;

#if USEFP16 == 1
	const DXGI_FORMAT pixel_format = DXGI_FORMAT_R16_FLOAT;

	DirectX::ScratchImage img{};

	HRESULT res = LoadFromDDSFile(FGD_path.c_str(), DDS_FLAGS_NONE, nullptr, img);

	ASSERT_SUCCEEDED(res);

	const void* data_ptr = img.GetPixels(); 
	const size_t RowPitchBytes = Dim * sizeof(uint16_t); //half precision.
#else
	const DXGI_FORMAT pixel_format = DXGI_FORMAT_R32_FLOAT;
	std::vector<float> data;
	LoadLUTFromFile<4>(FGD_path, data);


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

#ifdef DEBUG

	const auto max = std::max_element(new_data.cbegin(), new_data.cend());
	const auto min = std::min_element(new_data.cbegin(), new_data.cend());

	std::cout << std::format("Max element is {}, Min element is {}", *max, *min) << std::endl;

#endif


	//4D texture
	ASSERT(new_data.size() == Dim * Dim * Dim * (Dim/2));

	const void* data_ptr = new_data.data();

	//encode 4D LUT in 3D texture, Z dim encodes Z + W of index.
	//i.e. dim = 64 * 64 * (64 * 64) = 64 * 64 * 4096. 
	const size_t RowPitchBytes = Dim * sizeof(float);
#endif


	FGD.Create3D(RowPitchBytes, Dim, Dim, Dim * (Dim/2), pixel_format, data_ptr);

	return FGD;
}

Texture TransferMatrixResources::LoadGDLUTFromFile(const std::wstring& GD_path)
{
	constexpr size_t LUT_DIM = 64;
	std::vector<float> data;
	LoadLUTFromFile<2>(GD_path, data);

	Texture tex{};

	tex.Create2D(LUT_DIM * sizeof(float), LUT_DIM, LUT_DIM, DXGI_FORMAT_R32_FLOAT, data.data());
	

	return tex;
}

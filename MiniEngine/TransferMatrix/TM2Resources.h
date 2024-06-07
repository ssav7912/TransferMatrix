#pragma once
#include "../Core/PipelineState.h"
#include "../Core/DescriptorHeap.h"
#include "../Core/Texture.h"
#include "Texture3D.h"
#include "../Core/TextureManager.h"


class TM2Resources
{
public:
	TM2Resources() = default;
	TM2Resources(const std::string & FGD_path, const std::string & TIR_path);

	void Initialise(const std::string& FGD_path, const std::string& TIR_path);


	GraphicsPSO TM2PSO{ L"2-flux Transfer Matrix PSO" };

	Texture3D TIR_LUT;
	TextureRef FGD_LUT;

	static constexpr int32_t NUM_LAYERS = 2;
	static constexpr int32_t MAX_LAYERS = 5;
	
	//const void* TM2_PS = g_pTM2DielectricPS;
	//const void* TM2_VS = g_pDefaultVS;

	//static constexpr size_t sizeof_TM2PS = sizeof(g_pTM2DielectricPS);
	//static constexpr size_t sizeof_TM2VS = sizeof(g_pDefaultVS);

private:

	Texture3D LoadTIRLutFromFile(const std::string& TIR_path);

};

//Constant buffer for layer parameters
struct LayerConstants
{
	Math::Vector3 IORs[TM2Resources::MAX_LAYERS];
	Math::Vector3 Kappas[TM2Resources::MAX_LAYERS];
	float Roughs[TM2Resources::MAX_LAYERS];
	int32_t layers;
};


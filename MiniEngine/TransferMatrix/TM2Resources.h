#pragma once
#include "../Core/PipelineState.h"
#include "../Core/DescriptorHeap.h"
#include "../Core/Texture.h"
#include "../Core/TextureManager.h"

class TM2Resources
{
	TM2Resources(const std::string & FGD_path, const std::string & TIR_path);

public:
	GraphicsPSO TM2PSO{ L"2-flux Transfer Matrix PSO" };


private:

	TextureRef LoadTIRLutFromFile(const std::string& TIR_path);

	TextureRef TIR_LUT;
	TextureRef FGD_LUT;
};


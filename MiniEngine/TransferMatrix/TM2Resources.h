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

private:

	Texture3D LoadTIRLutFromFile(const std::string& TIR_path);

};


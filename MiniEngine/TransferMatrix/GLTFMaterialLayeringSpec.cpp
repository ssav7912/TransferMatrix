#include "GLTFMaterialLayeringSpec.h"

TransferMatrixMaterial GLTFMaterialLayeringSpec::ProcessMaterialLayer(json& material)
{
	TransferMatrixMaterial outMaterial{};

	if (material.find("TransferMatrix") != material.end())
	{
		for (int32_t i = 0; i < TransferMatrixResources::MAX_LAYERS; i++)
		{
			glTF::Asset::ReadTextureInfo
		}
	}

	return outMaterial;

}

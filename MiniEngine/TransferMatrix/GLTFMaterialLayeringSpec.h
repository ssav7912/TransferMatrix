#pragma once
#include "../Model/glTF.h"
#include "TransferMatrixResources.h"
#include "../Model/json.hpp"

struct TransferMatrixMaterial
{
	float IORs[TransferMatrixResources::MAX_LAYERS];
	float Kappas[TransferMatrixResources::MAX_LAYERS];
	float Sigma_S[TransferMatrixResources::MAX_LAYERS];
	float Sigma_K[TransferMatrixResources::MAX_LAYERS];
	float Depths[TransferMatrixResources::MAX_LAYERS];
	float Phase[TransferMatrixResources::MAX_LAYERS];
	float Roughs[TransferMatrixResources::MAX_LAYERS];

	enum {kIORs = 0,
		  kKappas = TransferMatrixResources::MAX_LAYERS, 
		  kSigma_S = TransferMatrixResources::MAX_LAYERS * 2, 
		  kSigma_K = TransferMatrixResources::MAX_LAYERS * 3, 
		  kDepths = TransferMatrixResources::MAX_LAYERS * 4, 
		  kPhase = TransferMatrixResources::MAX_LAYERS * 5, 
		  kRoughs = TransferMatrixResources::MAX_LAYERS * 6};

	static constexpr size_t kNumAttributes = 7;

	glTF::Texture* textures[TransferMatrixResources::NUM_LAYERS * kNumAttributes];
};

class GLTFMaterialLayeringSpec
{

	TransferMatrixMaterial ProcessMaterialLayer(json& material);

};


#pragma once
#include "Math/Vector.h"
#include "../ImGUI/imgui.h"
#include "GraphicsCore.h"
#include "d3d12.h"
#include "DescriptorHeap.h"
#include <format>
#include <vector>
#include <wrl.h>


class GUI
{
public:
	void Initialise(struct ID3D12Device* device, DescriptorHeap SRVDescriptorHeap);
	void Teardown();

	void LayerUI(int32_t num_layers, int32_t max_layers, bool UseTM6 = false);

	

	std::vector<Math::Vector3> IORs = {};
	std::vector<Math::Vector3> Kappas = {};
	std::vector<float> Roughs = {0.01f, 0.1f};

	DescriptorHeap FontHeap;

	static constexpr size_t VECTOR_DIM = 3;
	float IOR_DEFAULT[5][VECTOR_DIM] = { {1.0,1.0,1.0}, {1.5,1.5,1.5}, {1.0,1.0,1.0}, {1.0,1.0,1.0} };
	float KAPPA_DEFAULT[5][VECTOR_DIM] = { {1.0f, 0.1f, 0.1f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f } };
	float ROUGH_DEFAULT[5] = { 0.01f, 0.1f, 0.0f, 0.0f };

	float SIGMA_S_DEFAULT[5][VECTOR_DIM] = { 0 };
	float SIGMA_K_DEFAULT[5][VECTOR_DIM] = { 0 };
	float DEPTH_DEFAULT[5] = { 0 };
	float G_DEFAULT[5] = { 0 };


	int32_t NumSamples = 5;
	int32_t NumLayers = 2;


};


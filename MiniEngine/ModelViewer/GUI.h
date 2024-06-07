#pragma once
#include "Math/Vector.h"
#include "../ImGUI/imgui.h"
#include "GraphicsCore.h"
#include "d3d12.h"
#include "DescriptorHeap.h"
#include <format>
#include <vector>
#include <wrl.h>

static constexpr float IOR_DEFAULT[2][3] = { {1.0,1.0,1.0}, {1.5,1.5,1.5} };
static constexpr float KAPPA_DEFAULT[2][3] = { {1.0f, 0.1f, 0.1f }, { 0.0f, 0.0f, 0.0f } };
static constexpr float ROUGH_DEFAULT[2] = { 0.01f, 0.1f };

class GUI
{
public:
	void Initialise(struct ID3D12Device* device, DescriptorHeap SRVDescriptorHeap);
	void Teardown();

	void LayerUI(uint32_t num_layers, uint32_t max_layers);

	

	std::vector<Math::Vector3> IORs = {};
	std::vector<Math::Vector3> Kappas = {};
	std::vector<float> Roughs = {0.01f, 0.1f};

	DescriptorHeap FontHeap;

};


#include "GUI.h"
#include "../ImGUI/backends/imgui_impl_dx12.h"
#include "../ImGUI/backends/imgui_impl_win32.h"
#include "GraphicsCore.h"
#include "GameCore.h"
#include "d3d12.h"
#include "dxgi.h"


void GUI::Initialise(ID3D12Device* device, DescriptorHeap SRVDescriptorHeap)
{
	FontHeap = SRVDescriptorHeap;

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
	io.MouseDrawCursor = true;

	ImGui_ImplWin32_Init(GameCore::g_hWnd);
	ImGui_ImplDX12_Init(device, 1, DXGI_FORMAT_R11G11B10_FLOAT, FontHeap.GetHeapPointer(),
		FontHeap.GetHeapPointer()->GetCPUDescriptorHandleForHeapStart(), FontHeap.GetHeapPointer()->GetGPUDescriptorHandleForHeapStart());



}

void GUI::Teardown()
{
	ImGui_ImplWin32_Shutdown();
	ImGui_ImplDX12_Shutdown();
	ImGui::DestroyContext();

}

void GUI::LayerUI(int32_t num_layers, int32_t max_layers)
{
	num_layers = std::max(1,std::min(num_layers, max_layers));

	ImGui::SetNextWindowSize({ 600, 768 });
	ImGui::Begin("Layers");
	for (int32_t i = num_layers - 1; i >= 0; --i)
	{
		ImGui::Text(std::format("Layer {0} Params", i).c_str());



		ImGui::DragFloat3(std::format("IOR Layer {0}", i).c_str(), IOR_DEFAULT[i],0.1, 0.0f);
		ImGui::DragFloat3(std::format("Kappa Layer {0}", i).c_str(), KAPPA_DEFAULT[i], 0.1f, 0.0f);
		ImGui::SliderFloat(std::format("Roughness Layer {0}", i).c_str(), &ROUGH_DEFAULT[i], 0.0f, 1.0f);

		if (UseTM6)
		{
			ImGui::DragFloat3(std::format("Sigma S Layer {0}", i).c_str(), SIGMA_S_DEFAULT[i], 0.1, 0.0f);
			ImGui::DragFloat3(std::format("Sigma K Layer {0}", i).c_str(), SIGMA_K_DEFAULT[i], 0.1, 0.0f);
			ImGui::DragFloat(std::format("Depth Layer {0}", i).c_str(), &DEPTH_DEFAULT[i], 0.1, 0.0);
			ImGui::SliderFloat(std::format("Phase Media {0}", i).c_str(), &G_DEFAULT[i], 0.0f, 1.0f);
		}

		//IORs[i] = Math::Vector3(IOR_DEFAULT[i][0], IOR_DEFAULT[i][1], IOR_DEFAULT[i][2]);
		//Kappas[i] = Math::Vector3(KAPPA_DEFAULT[i][0], KAPPA_DEFAULT[i][1], KAPPA_DEFAULT[i][2]);
		IOR_DEFAULT[i][0] = std::max(0.0f, IOR_DEFAULT[i][0]);
		IOR_DEFAULT[i][1] = std::max(0.0f, IOR_DEFAULT[i][1]);
		IOR_DEFAULT[i][2] = std::max(0.0f, IOR_DEFAULT[i][2]);

		KAPPA_DEFAULT[i][0] = std::max(0.0f, KAPPA_DEFAULT[i][0]);
		KAPPA_DEFAULT[i][1] = std::max(0.0f, KAPPA_DEFAULT[i][1]);
		KAPPA_DEFAULT[i][2] = std::max(0.0f, KAPPA_DEFAULT[i][2]);

	}

	ImGui::DragInt("Number of Layers", &NumLayers, 1, max_layers);
	NumLayers = std::max(1, std::min(NumLayers, max_layers));
	ImGui::Checkbox("Use 6-flux Matrix", &UseTM6);

	ImGui::End();

}

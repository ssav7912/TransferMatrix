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

void GUI::LayerUI(uint32_t num_layers, uint32_t max_layers)
{
	IORs.assign(max_layers, { 1.0f, 1.0f, 1.0f });
	Kappas.assign(max_layers, { 0.0f, 0.0f, 0.0f });
	Roughs.assign(max_layers, 1.0f);

	ImGui::SetNextWindowSize({ 600, 768 });
	ImGui::Begin("Layers");
	for (int32_t i = num_layers - 1u; i >= 0; --i)
	{
		ImGui::Text(std::format("Layer {0} Params", i).c_str());

		constexpr size_t VECTOR_DIM = 3;
		float IOR[VECTOR_DIM] = { IOR_DEFAULT[i][0], IOR_DEFAULT[i][1], IOR_DEFAULT[i][2] };
		float Kappa[VECTOR_DIM] = { KAPPA_DEFAULT[i][0], KAPPA_DEFAULT[i][1], KAPPA_DEFAULT[i][2] };
		Roughs[i] = ROUGH_DEFAULT[i];

		ImGui::InputFloat3("IOR", IOR);
		ImGui::InputFloat3("Kappa", Kappa);
		ImGui::SliderFloat("Roughness", &Roughs[i], 0.0f, 1.0f);

		IORs[i] = Math::Vector3(IOR[0], IOR[1], IOR[2]);
		Kappas[i] = Math::Vector3(Kappa[0], Kappa[1], Kappa[2]);
		
	}

	ImGui::End();

}

#include "Texture3D.h"
#include "../Core/GraphicsCore.h"
#include "../Core/CommandContext.h"

void Texture3D::Create3D(size_t RowPitchBytes, size_t Width, size_t Height, size_t Depth, DXGI_FORMAT Format, const void* InitData)
{
	this->m_UsageState = D3D12_RESOURCE_STATE_COPY_DEST;
	
	m_Width = Width;
	m_Height = Height;
	m_Depth = Depth;

	D3D12_RESOURCE_DESC texDesc = {};
	texDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE3D;
	texDesc.Width = Width;
	texDesc.Height = Height;
	texDesc.DepthOrArraySize = Depth;
	texDesc.MipLevels = 1;
	texDesc.Format = Format;
	texDesc.SampleDesc.Count = 1;
	texDesc.SampleDesc.Quality = 0;
	texDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
	texDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

	D3D12_HEAP_PROPERTIES HeapProps;
	HeapProps.Type = D3D12_HEAP_TYPE_DEFAULT;
	HeapProps.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
	HeapProps.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
	HeapProps.CreationNodeMask = 1;
	HeapProps.VisibleNodeMask = 1;


	ASSERT_SUCCEEDED(g_Device->CreateCommittedResource(&HeapProps, D3D12_HEAP_FLAG_NONE, &texDesc,
		m_UsageState, nullptr, MY_IID_PPV_ARGS(m_pResource.ReleaseAndGetAddressOf())));

	m_pResource->SetName(L"3DTexture");

	D3D12_SUBRESOURCE_DATA texResource;
	texResource.pData = InitData;
	texResource.RowPitch = RowPitchBytes;
	texResource.SlicePitch = texResource.RowPitch * Height;

	CommandContext::InitializeTexture(*this, 1, &texResource);

	if (m_hCpuDescriptorHandle.ptr == D3D12_GPU_VIRTUAL_ADDRESS_UNKNOWN)
		m_hCpuDescriptorHandle = AllocateDescriptor(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
	g_Device->CreateShaderResourceView(m_pResource.Get(), nullptr, m_hCpuDescriptorHandle);

	
}

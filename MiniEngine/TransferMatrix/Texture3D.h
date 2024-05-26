#pragma once
#include "../Core/Texture.h"

class Texture3D : public Texture
{
public:
	void Create3D(size_t RowPitchBytes, size_t Width, size_t Height, size_t Depth, DXGI_FORMAT Format, const void* InitData);
};
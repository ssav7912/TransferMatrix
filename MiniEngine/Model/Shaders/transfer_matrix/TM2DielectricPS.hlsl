#include "HenyeyGreenstein.hlsli"


///@param incident incoming light direction in local coords
///@param outgoing outgoing light direction in local coords.
float3 sample(float3 incident, float3 outgoing, float3 iors[2], float alphas[2])
{
	//interface between air and arbitrary surface.
    float3 iors[2] = { float3(1.0, 1.0, 1.0f), float3(1.4f, 1.4, 1.4f)};
    float alphas[2] = { 0.0f, 1.0f };
	
    henyey_greenstein lobes[2];
    float3 lobes_incoming_direction[2];
	
	//TODO: Outgoing lobes
	
    float w[2];
    float w_sum = 0.0f;
	
    w[0] = (lobes[0].norm.x + lobes[0].norm.y + lobes[0].norm.z) / 3.0f;
    w_sum += w[0];
    w[1] = (lobes[1].norm.x + lobes[1].norm.y + lobes[1].norm.z) / 3.0f;
    w_sum += w[1];
	
	//TODO: random sampling?
    float sel_w = 1.0f * w_sum - w[0];
    int sel_i = 0;
    const int lobe_count = 2;
    for (sel_i = 0; sel_w > 0.f && sel_i < lobe_count; ++sel_i)
        sel_w -= w[sel_i + 1];
    
    const bool reflection = incident.z * lobes_incoming_direction[sel_i].z > 0.f;
    
    const float sel_a = hg_to_ggx(lobes[sel_i].asymmetry);
    
    //TODO: Microfacet distribution.
	
    //upper hemisphere 
    float3 sel_wi = lobes_incoming_direction[sel_i];
    sel_wi.z = abs(sel_wi.z);
    
    float H;  //ndf 
    
    float3 sel_wo = reflect(H, sel_wi);
    
	//Throughput
    
}

float4 main() : SV_TARGET
{
	
	
	
	return float4(1.0f, 1.0f, 1.0f, 1.0f);
}
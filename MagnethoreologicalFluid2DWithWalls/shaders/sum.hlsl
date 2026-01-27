[[vk::binding(0, 0)]] RWStructuredBuffer<double> x_0;
[[vk::binding(1, 0)]] RWStructuredBuffer<double> y_0;
[[vk::binding(2, 0)]] RWStructuredBuffer<double> z_0;

[[vk::binding(3, 0)]] RWStructuredBuffer<double> x_1;
[[vk::binding(4, 0)]] RWStructuredBuffer<double> y_1;
[[vk::binding(5, 0)]] RWStructuredBuffer<double> z_1;

[numthreads(1, 1, 1)]
void Main(uint3 DTid : SV_DispatchThreadID)
{
    x_1[DTid.x] = x_0[DTid.x] + 3.0;
    x_0[DTid.x] = x_1[DTid.x];
    y_1[DTid.x] = y_0[DTid.x];
    z_1[DTid.x] = z_0[DTid.x];
}
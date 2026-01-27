[[vk::binding(0, 0)]] RWStructuredBuffer<double> x_0;
[[vk::binding(1, 0)]] RWStructuredBuffer<double> y_0;
[[vk::binding(2, 0)]] RWStructuredBuffer<double> z_0;

[[vk::binding(12, 0)]] RWStructuredBuffer<int> dimensions;
[[vk::binding(13, 0)]] RWStructuredBuffer<int> length;
[[vk::binding(21, 0)]] RWStructuredBuffer<int> particle_0_array;
[[vk::binding(22, 0)]] RWStructuredBuffer<int> particle_1_array;
[[vk::binding(26, 0)]] RWStructuredBuffer<int> valid;

double sqrt_d(double x)
{
    double y = 1.0 / rsqrt((float) x);
    y = 0.5 * (y + x / y);
    y = 0.5 * (y + x / y);
    return y;
}

[numthreads(1, 1, 1)]
void Main(uint3 DTid : SV_DispatchThreadID)
{
    double r;
    double A = 1.0;
    double B = 100.0;
    double r_min = 1 - log(100.0) / B;
    double dot_product;
    int idx = DTid.x;
    int particle_0 = particle_0_array[idx];
    int particle_1 = particle_1_array[idx];
    
    if (dimensions[0] == 2)
    {
        double distances[3];
        distances[0] = sqrt_d((x_0[particle_1] - x_0[particle_0]) * (x_0[particle_1] - x_0[particle_0]) + (y_0[particle_1] - y_0[particle_0]) * (y_0[particle_1] - y_0[particle_0]));
        distances[1] = sqrt_d((x_0[particle_1] - x_0[particle_0] + (length[0])) * (x_0[particle_1] - x_0[particle_0] + (length[0])) + (y_0[particle_1] - y_0[particle_0]) * (y_0[particle_1] - y_0[particle_0]));
        distances[2] = sqrt_d((x_0[particle_1] - x_0[particle_0] - (length[0])) * (x_0[particle_1] - x_0[particle_0] - (length[0])) + (y_0[particle_1] - y_0[particle_0]) * (y_0[particle_1] - y_0[particle_0]));
        int index = 0;
        for (int j = 0; j < 3; j++)
        {
            if (distances[index] > distances[j])
            {
                index = j;
            }
        }
        r = distances[index];
    }
    else
    {
        double distances[9];
        distances[0] = sqrt_d((x_0[particle_1] - x_0[particle_0]) * (x_0[particle_1] - x_0[particle_0]) + (y_0[particle_1] - y_0[particle_0]) * (y_0[particle_1] - y_0[particle_0]) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[1] = sqrt_d((x_0[particle_1] - x_0[particle_0]) * (x_0[particle_1] - x_0[particle_0]) + (y_0[particle_1] - y_0[particle_0] - (length[0])) * (y_0[particle_1] - y_0[particle_0] - (length[0])) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[2] = sqrt_d((x_0[particle_1] - x_0[particle_0] + (length[0])) * (x_0[particle_1] - x_0[particle_0] + (length[0])) + (y_0[particle_1] - y_0[particle_0] - (length[0])) * (y_0[particle_1] - y_0[particle_0] - (length[0])) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[3] = sqrt_d((x_0[particle_1] - x_0[particle_0] + (length[0])) * (x_0[particle_1] - x_0[particle_0] + (length[0])) + (y_0[particle_1] - y_0[particle_0]) * (y_0[particle_1] - y_0[particle_0]) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[4] = sqrt_d((x_0[particle_1] - x_0[particle_0] + (length[0])) * (x_0[particle_1] - x_0[particle_0] + (length[0])) + (y_0[particle_1] - y_0[particle_0] + (length[0])) * (y_0[particle_1] - y_0[particle_0] + (length[0])) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[5] = sqrt_d((x_0[particle_1] - x_0[particle_0]) * (x_0[particle_1] - x_0[particle_0]) + (y_0[particle_1] - y_0[particle_0] + (length[0])) * (y_0[particle_1] - y_0[particle_0] + (length[0])) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[6] = sqrt_d((x_0[particle_1] - x_0[particle_0] - (length[0])) * (x_0[particle_1] - x_0[particle_0] - (length[0])) + (y_0[particle_1] - y_0[particle_0] + (length[0])) * (y_0[particle_1] - y_0[particle_0] + (length[0])) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[7] = sqrt_d((x_0[particle_1] - x_0[particle_0] - (length[0])) * (x_0[particle_1] - x_0[particle_0] - (length[0])) + (y_0[particle_1] - y_0[particle_0]) * (y_0[particle_1] - y_0[particle_0]) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        distances[8] = sqrt_d((x_0[particle_1] - x_0[particle_0] - (length[0])) * (x_0[particle_1] - x_0[particle_0] - (length[0])) + (y_0[particle_1] - y_0[particle_0] - (length[0])) * (y_0[particle_1] - y_0[particle_0] - (length[0])) + (z_0[particle_1] - z_0[particle_0]) * (z_0[particle_1] - z_0[particle_0]));
        int index = 0;
        for (int j = 0; j < 9; j++)
        {
            if (distances[index] > distances[j])
            {
                index = j;
            }
        }
        r = distances[index];
    }
    
    if (r < r_min)
    { valid[0] = 0;
    }
}
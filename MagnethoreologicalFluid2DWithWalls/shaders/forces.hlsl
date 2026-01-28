[[vk::binding(0, 0)]] RWStructuredBuffer<double> x_0;
[[vk::binding(1, 0)]] RWStructuredBuffer<double> y_0;
[[vk::binding(2, 0)]] RWStructuredBuffer<double> z_0;

[[vk::binding(6, 0)]] RWStructuredBuffer<double> forces_x;
[[vk::binding(7, 0)]] RWStructuredBuffer<double> forces_y;
[[vk::binding(8, 0)]] RWStructuredBuffer<double> forces_z;

[[vk::binding(9, 0)]] RWStructuredBuffer<double> magnetic_field;
[[vk::binding(12, 0)]] RWStructuredBuffer<int> dimensions;
[[vk::binding(13, 0)]] RWStructuredBuffer<int> length;
[[vk::binding(15, 0)]] RWStructuredBuffer<int> particles;
[[vk::binding(21, 0)]] RWStructuredBuffer<int> particle_0_array;
[[vk::binding(22, 0)]] RWStructuredBuffer<int> particle_1_array;
[[vk::binding(23, 0)]] RWStructuredBuffer<int> matrix_size;

double sqrt_d(double x){
    double y = 1.0 / rsqrt((float) x);
    y = 0.5 * (y + x / y);
    y = 0.5 * (y + x / y);
    return y;
}

[numthreads(32, 1, 1)]
void Main(uint3 DTid : SV_DispatchThreadID)
{
    double r;
    double A = 1.0;
    double B = 100.0;
    double x;
    double y;
    double z;
    double dot_product;
    int idx = DTid.x;
    int particle_0 = particle_0_array[idx];
    int particle_1 = particle_1_array[idx];
    
    if ( dimensions[0] == 2)
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

        x = (length[0]) * ((index == 1) - (index == 2)) + x_0[particle_1] - x_0[particle_0];
        y = y_0[particle_1] - y_0[particle_0];
        x = x / r;
        y = y / r;

        dot_product = magnetic_field[0] * x + magnetic_field[1] * y;

        forces_x[idx] = ((1 - 5 * (dot_product * dot_product)) / (r * r * r * r) + A * (double) exp(-B * (r - 1.0))) * x + (2.0 * dot_product * magnetic_field[0]) / (r * r * r * r);
        forces_y[idx] = ((1 - 5 * (dot_product * dot_product)) / (r * r * r * r) + A * (double) exp(-B * (r - 1.0))) * y + (2.0 * dot_product * magnetic_field[1]) / (r * r * r * r);
        forces_z[idx] = 0.0;
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

        x = ( length[0]) * ((index == 2) + (index == 3) + (index == 4) - (index == 6) - (index == 7) - (index == 8)) + x_0[particle_1] - x_0[particle_0];
        y = ( length[0]) * ((index == 4) + (index == 5) + (index == 6) - (index == 1) - (index == 2) - (index == 8)) + y_0[particle_1] - y_0[particle_0];
        z = z_0[particle_1] - z_0[particle_0];
        x = x / r;
        y = y / r;
        z = z / r;

        dot_product = magnetic_field[0] * x + magnetic_field[1] * y + magnetic_field[2] * z;

        forces_x[idx] = ((1 - 5 * (dot_product * dot_product)) / (r * r * r * r) + A * (double) exp(-B * (r - 1.0))) * x + (2.0 * dot_product * magnetic_field[0]) / (r * r * r * r);
        forces_y[idx] = ((1 - 5 * (dot_product * dot_product)) / (r * r * r * r) + A * (double) exp(-B * (r - 1.0))) * y + (2.0 * dot_product * magnetic_field[1]) / (r * r * r * r);
        forces_z[idx] = ((1 - 5 * (dot_product * dot_product)) / (r * r * r * r) + A * (double) exp(-B * (r - 1.0))) * z + (2.0 * dot_product * magnetic_field[2]) / (r * r * r * r);
    }
}
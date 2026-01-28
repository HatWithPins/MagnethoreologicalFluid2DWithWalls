[[vk::binding(0, 0)]] RWStructuredBuffer<double> x_0;
[[vk::binding(1, 0)]] RWStructuredBuffer<double> y_0;
[[vk::binding(2, 0)]] RWStructuredBuffer<double> z_0;

[[vk::binding(3, 0)]] RWStructuredBuffer<double> x_1;
[[vk::binding(4, 0)]] RWStructuredBuffer<double> y_1;
[[vk::binding(5, 0)]] RWStructuredBuffer<double> z_1;

[[vk::binding(6, 0)]] RWStructuredBuffer<double> forces_x;
[[vk::binding(7, 0)]] RWStructuredBuffer<double> forces_y;
[[vk::binding(8, 0)]] RWStructuredBuffer<double> forces_z;

[[vk::binding(10, 0)]] RWStructuredBuffer<int> mode;
[[vk::binding(11, 0)]] RWStructuredBuffer<int> phase;
[[vk::binding(12, 0)]] RWStructuredBuffer<int> dimensions;
[[vk::binding(13, 0)]] RWStructuredBuffer<int> length;
[[vk::binding(15, 0)]] RWStructuredBuffer<int> particles;
[[vk::binding(16, 0)]] RWStructuredBuffer<double> delta_t;
[[vk::binding(19, 0)]] RWStructuredBuffer<int> initial_indices_sum;
[[vk::binding(20, 0)]] RWStructuredBuffer<int> last_indices_sum;
[[vk::binding(23, 0)]] RWStructuredBuffer<int> matrix_size;
[[vk::binding(26, 0)]] RWStructuredBuffer<int> valid;
[[vk::binding(28, 0)]] RWStructuredBuffer<double> stress_array;
[[vk::binding(29, 0)]] RWStructuredBuffer<double> wall_velocity;

double sqrt_d(double x)
{
    double y = 1.0 / rsqrt((float) x);
    y = 0.5 * (y + x / y);
    y = 0.5 * (y + x / y);
    return y;
}

[numthreads(32, 1, 1)]
void Main(uint3 DTid : SV_DispatchThreadID)
{
    double r;
    double A = 1;
    double B = 100;
    double top_separation = 0.5;
    double top_repulsion = 0;
    double bottom_separation = 0.5;
    double bottom_repulsion = 0;
    double stress = 0;
    int particle_0 = DTid.x;
    int remainder_x;
    int initial_index_sum = initial_indices_sum[particle_0];
    int last_index_sum = last_indices_sum[particle_0];
    int index_sub = particle_0 - 1;

    x_1[particle_0] = 0;
    y_1[particle_0] = 0;
    z_1[particle_0] = 0;
    if ( dimensions[0] == 2)
    {
        top_separation = (length[0]) - y_0[particle_0] + 0.5;
        top_repulsion = -A * exp(-B * (top_separation - 1));
        bottom_separation = y_0[particle_0] + 0.5;
        bottom_repulsion = A * exp(-B * (bottom_separation - 1));
        y_1[particle_0] += top_repulsion + bottom_repulsion;
    }
    else
    {
        top_separation = (length[0]) - z_0[particle_0] + 0.5;
        top_repulsion = -A * exp(-B * (top_separation - 1));
        bottom_separation = z_0[particle_0] + 0.5;
        bottom_repulsion = A * exp(-B * (bottom_separation - 1));
        z_1[particle_0] += top_repulsion + bottom_repulsion;

        if ( mode[0] == 1)
        {
            stress = z_0[particle_0] * (wall_velocity[0]) / (length[0]);
        }
    }


    if (initial_index_sum < (matrix_size[0]))
    {
        for (int i = initial_index_sum; i < last_index_sum; i++)
        {
            x_1[particle_0] += forces_x[i];
            y_1[particle_0] += forces_y[i];
            z_1[particle_0] += forces_z[i];
        }
    }
    if (particle_0 != 0)
    {
        for (int i = 0; i < particle_0; i++)
        {
            x_1[particle_0] -= forces_x[index_sub];
            y_1[particle_0] -= forces_y[index_sub];
            z_1[particle_0] -= forces_z[index_sub];

            index_sub += (particles[0]) - 2 - i;
        }
    }

    stress_array[particle_0] = x_1[particle_0] * z_0[particle_0];
    x_1[particle_0] += stress;
    x_1[particle_0] = (delta_t[0]) * x_1[particle_0];
    y_1[particle_0] = (delta_t[0]) * y_1[particle_0];
    z_1[particle_0] = (delta_t[0]) * z_1[particle_0];

    r = sqrt_d(x_1[particle_0] * x_1[particle_0] + y_1[particle_0] * y_1[particle_0] + z_1[particle_0] * z_1[particle_0]);
    if (r > 0.03)
    {
        valid[0] = 0;
    }

    x_1[particle_0] += x_0[particle_0];
    y_1[particle_0] += y_0[particle_0];
    z_1[particle_0] += z_0[particle_0];

    remainder_x = x_1[particle_0] / (length[0]);
    x_1[particle_0] = x_1[particle_0] - remainder_x * (length[0]) + (length[0]) * (x_1[particle_0] < 0);
    if (dimensions[0] == 3)
    {
        int remainder_y = y_1[particle_0] / (length[0]);
        y_1[particle_0] = y_1[particle_0] - remainder_y * (length[0]) + (length[0]) * (y_1[particle_0] < 0);
    }
}
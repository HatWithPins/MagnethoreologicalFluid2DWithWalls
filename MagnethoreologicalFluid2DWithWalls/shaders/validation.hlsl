[[vk::binding(0, 0)]] RWStructuredBuffer<double> x_0;
[[vk::binding(1, 0)]] RWStructuredBuffer<double> y_0;
[[vk::binding(2, 0)]] RWStructuredBuffer<double> z_0;

[[vk::binding(3, 0)]] RWStructuredBuffer<double> x_1;
[[vk::binding(4, 0)]] RWStructuredBuffer<double> y_1;
[[vk::binding(5, 0)]] RWStructuredBuffer<double> z_1;

[[vk::binding(9, 0)]] RWStructuredBuffer<double> magnetic_field;
[[vk::binding(10, 0)]] RWStructuredBuffer<int> mode;
[[vk::binding(11, 0)]] RWStructuredBuffer<int> phase;
[[vk::binding(12, 0)]] RWStructuredBuffer<int> dimensions;
[[vk::binding(13, 0)]] RWStructuredBuffer<int> length;
[[vk::binding(14, 0)]] RWStructuredBuffer<double> field_direction;
[[vk::binding(15, 0)]] RWStructuredBuffer<int> particles;
[[vk::binding(16, 0)]] RWStructuredBuffer<double> delta_t;
[[vk::binding(17, 0)]] RWStructuredBuffer<double> original_delta_t;
[[vk::binding(18, 0)]] RWStructuredBuffer<double> t;
[[vk::binding(24, 0)]] RWStructuredBuffer<double> mason;
[[vk::binding(25, 0)]] RWStructuredBuffer<double> amplitude_relationship;
[[vk::binding(26, 0)]] RWStructuredBuffer<int> valid;
[[vk::binding(27, 0)]] RWStructuredBuffer<double> stress;
[[vk::binding(28, 0)]] RWStructuredBuffer<double> stress_array;

static const double PI = 3.14159265358979323846;
static const double TWO_PI = 6.28318530717958647692;
static const double HALF_PI = 1.57079632679489661923;

double reduce_angle(double x)
{
    double k = floor(x / TWO_PI + 0.5);
    return x - k * TWO_PI;
}

double sin_minimax(double x)
{
    double x2 = x * x;
    return x * (
        1.0
        + x2 * (-1.66666666666666324348e-1
        + x2 * (8.33333333332248946124e-3
        + x2 * (-1.98412698298579493134e-4
        + x2 * (2.75573137070700676789e-6
        + x2 * (-2.50507602534068634195e-8)))))
    );
}

double cos_minimax(double x)
{
    double x2 = x * x;
    return 1.0
        + x2 * (-0.5
        + x2 * (4.16666666666666019037e-2
        + x2 * (-1.38888888888741095749e-3
        + x2 * (2.48015872894767294178e-5
        + x2 * (-2.75573143513906633035e-7)))));
}

double sin_double(double x)
{
    x = reduce_angle(x);

    if (x > HALF_PI)
        return sin_minimax(PI - x);
    if (x < -HALF_PI)
        return -sin_minimax(PI + x);

    return sin_minimax(x);
}

double cos_double(double x)
{
    x = reduce_angle(x);

    if (x > HALF_PI)
        return -cos_minimax(PI - x);
    if (x < -HALF_PI)
        return -cos_minimax(PI + x);

    return cos_minimax(x);
}

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
    double norm = 1.0;
    double pi = 3.14159265359;
    int particle_0 = DTid.x;
    if ( valid[0] == 1)
    {
        x_0[particle_0] = x_1[particle_0];
        y_0[particle_0] = y_1[particle_0];
        z_0[particle_0] = z_1[particle_0];
        double volume = (length[0]) * (length[0]) * (length[0]);

        if (particle_0 == 0)
        {
            t[0] +=  delta_t[0];
            delta_t[0] =  original_delta_t[0];

            if (dimensions[0] == 2)
            {
                norm = sqrt_d((amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) * (amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) + 1);
                magnetic_field[2] = 0;
                magnetic_field[1] = 1.0 / norm;
                magnetic_field[0] = (amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) * magnetic_field[1];
            }
            else
            {
                norm = sqrt_d((amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) * (amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) + 1);
                magnetic_field[2] = 1.0 / norm * (mode[0] != 1) + (mode[0] == 1);
                magnetic_field[1] = (amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) * magnetic_field[2] * sin_double(pi * (field_direction[0]) / 360) * (mode[0] != 1);
                magnetic_field[0] = (amplitude_relationship[0]) * sin_double((mason[0]) * (t[0])) * magnetic_field[2] * cos_double(pi * (field_direction[0]) / 360) * (mode[0] != 1);

                if (mode[0] == 1)
                {
                    stress[0] = 0;
                    for (int i = 0; i <  particles[0]; i++)
                    {
                        stress[0] += stress_array[i];
                    }
                    stress[0] = -( stress[0]) / volume;
                }
            }
        }
    }
    else if (particle_0 == 0)
    {
        delta_t[0] =  delta_t[0] / 2;
        valid[0] = 1;
    }
}
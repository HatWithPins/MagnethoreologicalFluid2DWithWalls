# MagnethoreologicalFluid2DWithWalls

I will expand this later, but for now, let's just say that this code is meant for simulating magnetorheological
fluids using OpenCL. Code still needs some cleaning for making it more easy to use.

I coded this in Visual Studio, using Windows and CMake. 
Don't know exactly if code needs something else for running in Linux or Mac systems.

And one quick note about OpenCL, if you have a NVidia card, just install CUDA and you are fine to go.
If using an AMD card, like I did, you need [this library](https://github.com/GPUOpen-LibrariesAndSDKs/OCL-SDK/releases).
I'm writing this on July 2021, if you found this much later in the future, it's possible that a different library 
is needed for compile and run OpenCL programs in AMD cards.

## Usage
``MagnethorologicalFluid2dWithWalls repetitions=2 particles=400 concentration=0.07 ma=0.4 ar=1.5 dimensions=3 field_direction=0``

**repetitions** can be any number between 1 and 5. It's the number of parallel simulations.\
**particles** is the number of particles per simulation, minimum 2. This variable is a 32 bit int.\
**concentration** is the desired volumetric concentration.\
**ma** is Mason number. If it's 0, then simulation only runs with DC field.\
**ar** is the amplitude relationship between DC field and perturbation. Must be greater than 0.\
**dimensions** can be 2 or 3.\
**field_direction** can be 0, for DC field oriented in the y axis, or 1, for DC field oriented in the z axis.
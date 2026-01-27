# CMakeList.txt : CMake project for MagnethoreologicalFluid2DWithWalls, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

find_package(Vulkan)
find_package (OpenCL)
include_directories( ${OpenCL_INCLUDE_DIR} )
LINK_DIRECTORIES(${OpenCL_LIBRARY})
if( NOT OpenCL_FOUND )
	MESSAGE( SEND_ERROR "OpenCL not found. Install SDK from NVidia (CUDA), AMD or Intel" )
endif( NOT OpenCL_FOUND )

# Add source to this project's executable.
add_executable (MagnethoreologicalFluid2DWithWalls "MagnethoreologicalFluid2DWithWalls.cpp" "MagnethoreologicalFluid2DWithWalls.h" "box.h" "dynamics.h" "analysis.h" "dynamics.cpp" "analysis.cpp" "box.cpp"  "VulkanContext.h" "VulkanContext.cpp"   "SimulationContext.h" "SimulationContext.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET MagnethoreologicalFluid2DWithWalls PROPERTY CXX_STANDARD 20)
endif()

 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/forces.cl
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sum.cl
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/distances.cl
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/validation.cl
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/positions)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/analysis)

add_custom_target(posiciones ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory positions)
add_custom_target(analisis ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory analysis)

add_custom_command(
	OUTPUT "${CMAKE_BINARY_DIR}/shaders/forces.spv"
	COMMAND $ENV{VULKAN_SDK}/bin/dxc.exe -T cs_6_2 -E "Main" -spirv -fvk-use-dx-layout -fspv-target-env=vulkan1.1 -Fo "${CMAKE_BINARY_DIR}/shaders/forces.spv" "shaders/forces.hlsl"
	DEPENDS "shaders/forces.hlsl"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Buiding Force Shader"
)
add_custom_command(
	OUTPUT "${CMAKE_BINARY_DIR}/shaders/distances.spv"
	COMMAND $ENV{VULKAN_SDK}/bin/dxc.exe -T cs_6_2 -E "Main" -spirv -fvk-use-dx-layout -fspv-target-env=vulkan1.1 -Fo "${CMAKE_BINARY_DIR}/shaders/distances.spv" "shaders/distances.hlsl"
	DEPENDS "shaders/distances.hlsl"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Buiding Distances Shader"
)
add_custom_command(
	OUTPUT "${CMAKE_BINARY_DIR}/shaders/sum.spv"
	COMMAND $ENV{VULKAN_SDK}/bin/dxc.exe -T cs_6_2 -E "Main" -spirv -fvk-use-dx-layout -fspv-target-env=vulkan1.1 -Fo "${CMAKE_BINARY_DIR}/shaders/sum.spv" "shaders/sum.hlsl"
	DEPENDS "shaders/sum.hlsl"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Buiding Sum Shader"
)
add_custom_command(
	OUTPUT "${CMAKE_BINARY_DIR}/shaders/validation.spv"
	COMMAND $ENV{VULKAN_SDK}/bin/dxc.exe -T cs_6_2 -E "Main" -spirv -fvk-use-dx-layout -fspv-target-env=vulkan1.1 -Fo "${CMAKE_BINARY_DIR}/shaders/validation.spv" "shaders/validation.hlsl"
	DEPENDS "shaders/validation.hlsl"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Buiding Validation Shader"
)
add_custom_target(ComputeShader DEPENDS "${CMAKE_BINARY_DIR}/shaders/forces.spv" "${CMAKE_BINARY_DIR}/shaders/sum.spv" "${CMAKE_BINARY_DIR}/shaders/validation.spv" "${CMAKE_BINARY_DIR}/shaders/distances.spv")
# TODO: Add tests and install targets if needed.
target_link_libraries(MagnethoreologicalFluid2DWithWalls OpenCL::OpenCL Vulkan::Vulkan)
add_dependencies(MagnethoreologicalFluid2DWithWalls ComputeShader)
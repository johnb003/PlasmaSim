cmake_minimum_required(VERSION 3.4)

project (PlasmaSim)

add_subdirectory(adMath)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PlasmaSim_SOURCE_DIR}/cmake")
set(SDL2_PATH ${PlasmaSim_SOURCE_DIR}/../sdl ${PlasmaSim_SOURCE_DIR}/../sdl/make)
set(CUDA_TOOLKIT_ROOT_DIR /Developer/NVIDIA/CUDA-9.0)

find_package(SDL2 REQUIRED)
find_package(CUDA REQUIRED)
find_package(OpenGL REQUIRED)

message( STATUS "CUDA_VERSION_MAJOR: " ${CUDA_VERSION_MAJOR} )
message( STATUS "CUDA_VERSION_MINOR: " ${CUDA_VERSION_MINOR} )
message( STATUS "CUDA_VERSION: " ${CUDA_VERSION} )
message( STATUS "CUDA_VERSION_STRING: " ${CUDA_VERSION_STRING} )
message( STATUS "CUDA_TOOLKIT_ROOT_DIR: " ${CUDA_TOOLKIT_ROOT_DIR} )
message( STATUS "CUDA_SDK_ROOT_DIR: " ${CUDA_SDK_ROOT_DIR} )
message( STATUS "CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS} )
message( STATUS "CUDA_LIBRARIES: " ${CUDA_LIBRARIES} )
message( STATUS "CUDA_CUFFT_LIBRARIES: " ${CUDA_CUFFT_LIBRARIES} )
message( STATUS "CUDA_CUBLAS_LIBRARIES: " ${CUDA_CUBLAS_LIBRARIES} )
message( STATUS "CUDA_cupti_LIBRARY: " ${CUDA_cupti_LIBRARY} )
message( STATUS "CUDA_curand_LIBRARY: " ${CUDA_curand_LIBRARY} )
message( STATUS "CUDA_cusparse_LIBRARY: " ${CUDA_cusparse_LIBRARY} )
message( STATUS "CUDA_npp_LIBRARY: " ${CUDA_npp_LIBRARY} )
message( STATUS "CUDA_nppc_LIBRARY: " ${CUDA_nppc_LIBRARY} )
message( STATUS "CUDA_nppi_LIBRARY: " ${CUDA_nppi_LIBRARY} )
message( STATUS "CUDA_npps_LIBRARY: " ${CUDA_npps_LIBRARY} )



include_directories("${SDL2_INCLUDE_DIR}")
include_directories("${CUDA_INCLUDE_DIRS}")
include_directories(${PlasmaSim_SOURCE_DIR}/../adMath/include)

cuda_add_executable(PlasmaSim
	src/CurrentLineField.cpp
	src/CurrentLineVectorField.cpp
	src/FieldAccumulator.cpp
	src/FieldRenderer.cpp
	src/ParticleSystem.cpp
	src/ParticleSimulationCuda.cu
	src/ParticleSimulation.cpp
	src/sdl_gl_texture.cpp
	src/engine.cpp
	src/font.cpp
	src/main.cpp
	src/tweaks.cpp
)

target_link_libraries (PlasmaSim
	"${SDL2_LIBRARY}"
	"${CUDA_LIBRARIES}"
	"${OPENGL_LIBRARIES}"
	adMath)



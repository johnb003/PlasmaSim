cmake_minimum_required(VERSION 3.4)

# Create a library called adMath which includes the source files.
# The extension is already found. Any number of sources could be listed here.

add_library (adMath
	src/adMatrix.cpp
	src/adVector.cpp
	src/adQuaternion.cpp
	)

# Make sure the compiler can find include files for our library
# when other libraries or executables link to it.
target_include_directories (adMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
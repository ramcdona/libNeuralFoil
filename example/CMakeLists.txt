cmake_minimum_required( VERSION 3.23 )

project( lnfExampleProject )

set( CMAKE_CXX_STANDARD 17 )

include( FetchContent )

# Since this example is bundled inside the libNeuralFoil repository, point FetchContent at the enclosing
# directory.  However in normal circumstances, you should point FetchContent at the libNeuralFoil repository
# or another source.
FetchContent_Declare( libneuralfoil
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..
        # GIT_REPOSITORY https://github.com/ramcdona/libNeuralFoil
        # GIT_TAG origin/main
)

FetchContent_MakeAvailable( libneuralfoil )

add_executable( lnfexample
        lnfexample.cpp
)

target_link_libraries( lnfexample PUBLIC
        neuralfoil
)

get_target_property( neuralfoil_RESOURCE_DIR neuralfoil RESOURCE_DIR )

# Copy the NeuralFoil data files from the included NeuralFoil source to the build directory.  You will
# need to include the desired data files with your program.
add_custom_command ( TARGET lnfexample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${neuralfoil_RESOURCE_DIR}/
        ${CMAKE_CURRENT_BINARY_DIR}
)

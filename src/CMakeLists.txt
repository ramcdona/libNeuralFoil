cmake_minimum_required( VERSION 3.23 )

set( CMAKE_CXX_STANDARD 17 )

add_library( neuralfoil
        neuralfoil.cpp
        airfoilnormalizer.cpp
        kulfancst.cpp
        matsolve.cpp
        matsolve.h
)

set_target_properties( neuralfoil PROPERTIES
        PUBLIC_HEADER "../include/neuralfoil.h"
)

target_include_directories( neuralfoil PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/external
)

target_include_directories( neuralfoil PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

# Define resource destination variable for when project is installed, this will result
# in neuralfoil_RESOURCE_DIR being defined in a project including this library via
# FindPackage
set_target_properties( neuralfoil PROPERTIES
        RESOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external/NeuralFoil/neuralfoil/nn_weights_and_biases
)

set_property( TARGET neuralfoil APPEND PROPERTY
        EXPORT_PROPERTIES RESOURCE_DIR
)

target_link_libraries( neuralfoil PUBLIC npy::npy )

install( TARGETS neuralfoil
        EXPORT neuralfoiltargets
)

install( DIRECTORY
        ${CMAKE_SOURCE_DIR}/external/NeuralFoil/neuralfoil/nn_weights_and_biases/
        DESTINATION share/neuralfoil/resources
)

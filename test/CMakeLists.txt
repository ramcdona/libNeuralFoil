cmake_minimum_required( VERSION 3.23 )

set( CMAKE_CXX_STANDARD 17 )

add_executable( test
        test.cpp
        neuralfoil_test.cpp
        neuralfoil_test.h
)

target_link_libraries( test PUBLIC
        neuralfoil
)

get_target_property( neuralfoil_RESOURCE_DIR neuralfoil RESOURCE_DIR )

ADD_CUSTOM_COMMAND ( TARGET test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${neuralfoil_RESOURCE_DIR}/
        ${CMAKE_CURRENT_BINARY_DIR}
)

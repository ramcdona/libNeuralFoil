CMAKE_MINIMUM_REQUIRED(VERSION 3.24)

PROJECT( libneuralfoil
        VERSION 0.3.0
        LANGUAGES CXX
)

option( LNF_BUILD_DOCUMENTATION "Build API documentation with Doxygen" OFF )

if( LNF_BUILD_DOCUMENTATION )
    find_package( Doxygen )
    if( DOXYGEN_FOUND )
        # Set Doxygen input and output directories
        set( DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in )
        set( DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile )

        # Configure Doxyfile from template
        configure_file( ${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY )

        # Add target for building documentation
        add_custom_target( docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
        )

        message( STATUS "Doxygen found. Documentation target 'docs' available." )
        message( STATUS "Build documentation with: cmake --build . --target docs" )
    else()
        message( WARNING "Doxygen not found. Documentation target will not be available." )
    endif()
else()
    message( WARNING "Documentation not enabled." )
endif()

ADD_SUBDIRECTORY( external )

ADD_SUBDIRECTORY( src )

IF( LNF_BUILD_TESTS )
    ADD_SUBDIRECTORY( test )
ENDIF()

install( EXPORT neuralfoiltargets
        FILE NeuralFoilTargets.cmake
        DESTINATION lib/cmake/NeuralFoil
)

include( CMakePackageConfigHelpers )
# generate the config file that includes the exports
configure_package_config_file( ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/NeuralFoilConfig.cmake"
        INSTALL_DESTINATION "lib/cmake/NeuralFoil"
)

write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/NeuralFoilConfigVersion.cmake"
        VERSION "${libneuralfoil_VERSION_MAJOR}.${libneuralfoil_VERSION_MINOR}"
        COMPATIBILITY AnyNewerVersion
)

install( FILES
        ${CMAKE_CURRENT_BINARY_DIR}/NeuralFoilConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/NeuralFoilConfigVersion.cmake
        DESTINATION lib/cmake/NeuralFoil
)

export( EXPORT neuralfoiltargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/NeuralFoilTargets.cmake"
)

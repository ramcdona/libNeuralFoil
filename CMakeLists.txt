CMAKE_MINIMUM_REQUIRED(VERSION 3.24)

PROJECT( libneuralfoil
        VERSION 0.1.0
        LANGUAGES CXX
)

ADD_SUBDIRECTORY( external )

ADD_SUBDIRECTORY( src )

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

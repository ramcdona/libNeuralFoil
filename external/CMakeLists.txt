include( FetchContent )

FetchContent_Declare( libnpy
    URL ${CMAKE_CURRENT_SOURCE_DIR}/libnpy-e76b1dae9fe7.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable( libnpy )

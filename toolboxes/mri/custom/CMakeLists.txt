    if (WIN32)
        add_definitions(-D__BUILD_GADGETRON_CUSTOM__)
    endif ()

    add_library(gadgetron_toolbox_custom  INTERFACE)

    target_link_libraries(gadgetron_toolbox_custom INTERFACE gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_log )

    target_include_directories(gadgetron_toolbox_custom
            INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR})

    install(FILES
            CustomExport.h
            CustomEPIReconXObjectTrapezoid.h
            DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

    # install(TARGETS epi DESTINATION lib)

    #if (GTEST_FOUND)
    #    add_subdirectory(ut)
    #endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)
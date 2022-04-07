if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CUSTOM__)
endif ()


add_library(gadgetron_custom SHARED
		gadgetron_custom_export.h
#	    CustomEPIReconXGadget.h
#		CustomEPIReconXGadget.cpp
		CustomGenericReconCartesianGrappaGadget.h
		CustomGenericReconCartesianGrappaGadget.cpp
		SMSReconCartesianGrappaGadget.h
		SMSReconCartesianGrappaGadget.cpp
)

set_target_properties(gadgetron_custom PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_custom
#  gadgetron_epi
  gadgetron_core
  gadgetron_mricore
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_log
#  gadgetron_toolbox_epi
  gadgetron_toolbox_custom
  gadgetron_core_writers
  gadgetron_toolbox_image_analyze_io
  gadgetron_toolbox_denoise
  gadgetron_toolbox_hostutils
  gadgetron_toolbox_cpuoperator
  gadgetron_toolbox_cpuklt
  gadgetron_toolbox_mri_core
)


install(FILES
#  CustomEPIReconXGadget.h
  CustomGenericReconCartesianGrappaGadget.h
  SMSReconCartesianGrappaGadget.h
  gadgetron_custom_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_custom
#	EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


target_include_directories(gadgetron_custom
		PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

#if (ZFP_FOUND)
#	message("ZFP Found")
#	add_definitions(-DGADGETRON_COMPRESSION_ZFP)
#	target_include_directories(
#			gadgetron_mricore
#			PUBLIC
#			${ZFP_INCLUDE_DIR})
#
#	target_link_libraries(gadgetron_custom ${ZFP_LIBRARIES})
#else ()
#	message("ZFP NOT Found")
#endif ()

#if (Boost_PYTHON3_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
#	target_link_libraries(gadgetron_custom
#			gadgetron_toolbox_python
#			)
#
#	install(FILES
#			${gadgetron_python_models_files}
#			DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)
#endif ()


set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)

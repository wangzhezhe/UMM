include(CMakeFindDependencyMacro)

set(REQUIRED_VARS)

if(TRUE)
	find_dependency(ZLIB)
	list(APPEND REQUIRED_VARS ZLIB_FOUND)
endif()

if(TRUE)
	find_dependency(ZSTD)
	list(APPEND REQUIRED_VARS ZSTD_FOUND)
endif()

if()
	find_dependency(NVCOMP)
	list(APPEND REQUIRED_VARS NVCOMP_FOUND)
endif()

if()
	find_dependency(CUDAToolkit)
	list(APPEND REQUIRED_VARS CUDAToolkit_FOUND)
endif()

if(FALSE)
	find_dependency(MOAB)
	list(APPEND REQUIRED_VARS MOAB_FOUND)
endif()

if(TRUE)
	find_dependency(OpenMP)
	list(APPEND REQUIRED_VARS OpenMP_FOUND)
endif()

include(FindPackageHandleStandardArgs)
set(${CMAKE_FIND_PACKAGE_NAME}_CONFIG ${CMAKE_CURRENT_LIST_FILE})
find_package_handle_standard_args(
	${CMAKE_FIND_PACKAGE_NAME}
	REQUIRED_VARS ${REQUIRED_VARS}
	CONFIG_MODE
)

if(NOT TARGET mgard::mgard)
  include("${CMAKE_CURRENT_LIST_DIR}/mgard-targets.cmake")
endif()

set(MGARD_LIBRARIES mgard::mgard)
set(MGARD_INCLUDE_DIRS
  $<TARGET_PROPERTY:mgard::mgard,INTERFACE_INCLUDE_DIRECTORIES>
)

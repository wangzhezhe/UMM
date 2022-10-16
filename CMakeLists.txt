cmake_minimum_required (VERSION 3.2)
project(UMM LANGUAGES CXX)

find_package(VTK COMPONENTS
IOLegacy
IOXML
IOImage
vtkFiltersCore
)

add_executable(addVar addVar.cpp)
target_link_libraries(addVar ${VTK_LIBRARIES})

add_executable(interp interp.cpp)
target_link_libraries(interp ${VTK_LIBRARIES})

find_package(mgard)
add_executable(resampleCompress resampleCompress.cpp)
target_link_libraries(resampleCompress ${VTK_LIBRARIES} mgard::mgard)

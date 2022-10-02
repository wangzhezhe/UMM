cmake_minimum_required (VERSION 3.2)
project(UMM LANGUAGES CXX)

find_package(VTK COMPONENTS
IOLegacy
IOXML
IOImage
vtkFiltersCore
)

add_executable(DataGenerator DataGenerator.cpp)
target_link_libraries(DataGenerator ${VTK_LIBRARIES})

find_package(mgard)
add_executable(NaiveMgard NaiveMgard.cpp)
target_link_libraries(NaiveMgard ${VTK_LIBRARIES} mgard::mgard)

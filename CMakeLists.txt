cmake_minimum_required(VERSION 3.9)
project (memgraph)

add_subdirectory(vendor/osmpbf osmpbf)

set(SOURCES_CPP
	src/Graph.cpp
)

add_library(${PROJECT_NAME} STATIC ${SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} PUBLIC osmpbf)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

#Build the example
get_directory_property(hasParent PARENT_DIRECTORY)
if(NOT hasParent)
	add_executable(${PROJECT_NAME}_example example/main.cpp)
	target_link_libraries(${PROJECT_NAME}_example ${PROJECT_NAME})
endif()

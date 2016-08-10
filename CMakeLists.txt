cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
	${CMAKE_SOURCE_DIR}/cmake
)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(MEMGRAPH_SOURCES_CPP
		${CMAKE_CURRENT_SOURCE_DIR}/src/Graph.cpp
	PARENT_SCOPE)
	set(MEMGRAPH_INCLUDE_DIR
		${CMAKE_CURRENT_SOURCE_DIR}/include
	PARENT_SCOPE)
else()#hasParent


	project (memgraph)
	find_package(Protobuf REQUIRED)
	find_package(ZLIB REQUIRED)
	# find_package(LIBRT)

	include_directories(
		${CMAKE_SOURCE_DIR}
		${CMAKE_CURRENT_BINARY_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/include
		${CMAKE_CURRENT_SOURCE_DIR}/vendor/osmpbf
		${CMAKE_CURRENT_SOURCE_DIR}/vendor/osmpbf/osmpbf/include
	)

	set(MY_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4")
	set(MY_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4")
	set(DEBUG_FLAGS "${DEBUG_FLAGS}")#-Wconversion

	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
	set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")

	add_subdirectory(vendor/osmpbf osmpbf)

	SET(MY_LINK_LIBS
		${osmpbf_BINARY_DIR}/libosmpbf.a
		${PROTOBUF_LIBRARIES}
		${ZLIB_LIBRARIES}
		${LIBRT_LIBRARIES}
	)

	set(SOURCES_CPP
		src/Graph.cpp
		src/main.cpp
	)

	# The executable itself.
	add_executable(memgraph ${SOURCES_CPP} ${SOURCES_MOC_CPP})
	target_link_libraries(memgraph ${MY_LINK_LIBS})
	add_dependencies(memgraph osmpbf)
	
endif()##hasParent

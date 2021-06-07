cmake_minimum_required(VERSION 3.14)
project(CluES)

#set(CMAKE_CXX_STANDARD 17)

#******************************************************************************************  Google-tests
#configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#    message(FATAL_ERROR "Build step for googletest failed: ${result}")
#endif()
#
## Add googletest directly to our build. This defines
## the gtest and gtest_main targets.
#add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
#        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
#        EXCLUDE_FROM_ALL)
#******************************************************************************************  Google-tests



SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_DEBUG -std=c++17 -Wreturn-type -O0") # compilation for debugging
#SET(CMAKE_CXX_FLAGS "-std=c++17 -O0") # compilation for profiling

#remove -pthread
#SET(CMAKE_CXX_FLAGS "-static -std=c++17 -O3") # turbo fast
#SET(CMAKE_CXX_FLAGS "-static -std=c++17 -lm -O3 -DNDEBUG") # turbo fast - disabling assertions
#SET(CMAKE_CXX_FLAGS "-m64 -lm -O0 -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_SANITIZE_VECTOR -D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG_ASSERT -D_GLIBCXX_DEBUG_PEDASSERT -D_GLIBCXX_DEBUG_VERIFY -pedantic -march=native -ffast-math -Wextra") # compilation for faster execution -O2 uses compiler optimizations, -O0 uses no optimizations at all


#***************************** Cluster Editing

# QUICK_BUILD is a makro that excludes some dependancy-heavy code, such as GraphUtils::writeBasicGraphStatistics()
#add_compile_definitions(QUICK_BUILD)

include_directories( include
        include/graphs
        include/combinatorics
        include/utils
        include/CONTESTS/PACE21/*
        include/CONTESTS/PACE21/*/*
        include/CONTESTS/PACE21/*/*/*
        )

file(GLOB SOURCES
        src/graphs/GraphReader.cpp
        src/graphs/generators/GraphGenerator.cpp
        src/graphs/GraphUtils.cpp
        src/graphs/GraphInducer.cpp
        src/graphs/flow/UnitFlow2.cpp
        src/graphs/flow/MaxFlow.cpp
        src/graphs/graphtraversals/BFS.cpp

#        below section is for GraphUtils::writeBasicGraphStatistics
        src/graphs/graphtraversals/DFS.cpp
        src/graphs/trees/Tree.cpp
        src/graphs/components/BridgesAndArtPoints.cpp
        src/graphs/components/ConnctedComponents.cpp
        src/graphs/cliques/CliqueExtension.cpp
        src/graphs/cliques/CliqueUtils.cpp
        src/graphs/GraphTrimmer.cpp

        src/combinatorics/CombinatoricUtils.cpp

        src/datastructures/FAU.cpp

        src/utils/RandomNumberGenerators.cpp
        src/utils/StandardUtils.cpp
        src/utils/TimeMeasurer.cpp

        src/CONTESTS/PACE21/*.cpp
        src/CONTESTS/PACE21/kernelization/*.cpp
        src/CONTESTS/PACE21/heur/*.cpp
        src/CONTESTS/PACE21/*/EOCreators/*.cpp
        src/CONTESTS/PACE21/*/SwapCandidates/*.cpp
        src/CONTESTS/PACE21/*/SwpCndCreators/*.cpp
        src/CONTESTS/PACE21/*/StateImprovers/*.cpp

        ) # include all source files

#***************************** Cluster Editing




#*************** adding tests
#enable_testing()
#file(GLOB TESTS ${SOURCES}
##        "src/graphs/GraphUtils.cpp"
#        "src/graphs/unit_tests/test_GraphTrimmer.cpp"
#        "src/graphs/GraphTrimmer.cpp"
##        "src/CONTESTS/PACE21/kernelization/CEKernelizer.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_CEKernelizer.cpp"
##        "src/CONTESTS/PACE21/heur/EOCreators/ComponentExpansion.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_ComponentExpansion.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_State.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_SwapCandidate.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_SwpCndNode.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_SwpCndEdge.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_SwpCndTriangle.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_ConvexHullTrickDynamic.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_SwpCndEO.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_ComponentExpansionRepulsion.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_ComponentExpansionAttraction.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_Solver.cpp"
#        "src/CONTESTS/PACE21/unit_tests/test_NodeEdgeGreedy.cpp"
#        )
#
#add_executable(Tests ${TESTS})
#target_link_libraries(Tests gtest_main)
#target_link_libraries(Tests gmock_main)
#*************** adding tests

add_executable(CluES "main.cpp" ${SOURCES})
#target_link_libraries(CluES gtest_main) # necessary to write tests in the same .cpp file as source code

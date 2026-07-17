string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()

#set(MPICC "cc")
#set(MPICXX "CC")
#set(MPIFC "ftn")
#set(SCC "icx")
set(SCXX "icpx")
#set(SFC "ifx")

# fp-model source does not work on in mpicxx, this seems to be a compiler bug, need to manually switch to precise
set(CMAKE_CXX_FLAGS " ") # hardcode it here to blank, then try to do same things as in intel.cmake
if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=precise") # and manually add precise
#message(STATUS "ndk CXXFLAGS=${CXXFLAGS}")

string(APPEND CMAKE_Fortran_FLAGS " -fp-model consistent -fimf-use-svml")
string(APPEND CMAKE_CXX_FLAGS " -fp-model consistent")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -qno-opt-dynamic-align")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lpthread")

string(APPEND SPIO_CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

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

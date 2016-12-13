include(ExternalProject)
set(radiation_LIBRARIES "")

# Find Eigen.
find_package( Eigen3 REQUIRED )
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
list(APPEND radiation_LIBRARIES ${EIGEN3_LIBRARIES})

# Find Ceres.
find_package( Ceres REQUIRED )
include_directories(SYSTEM ${CERES_INCLUDE_DIRS})
list(APPEND radiation_LIBRARIES ${CERES_LIBRARIES})

# Find OpenGL.
find_package( OpenGL REQUIRED )
include_directories(SYSTEM ${OPENGL_INCLUDE_DIRS})
list(APPEND radiation_LIBRARIES ${OPENGL_LIBRARIES})

# Find GLUT.
find_package( GLUT REQUIRED )
include_directories(SYSTEM ${GLUT_INCLUDE_DIRS})
list(APPEND radiation_LIBRARIES ${GLUT_LIBRARIES})

# Find Google-gflags.
#include("cmake/External/gflags.cmake")
include("cmake/Modules/FindGflags.cmake")
include_directories(SYSTEM ${GFLAGS_INCLUDE_DIRS})
list(APPEND radiation_LIBRARIES ${GFLAGS_LIBRARIES})

# Find Google-glog.
#include("cmake/External/glog.cmake")
include("cmake/Modules/FindGlog.cmake")
include_directories(SYSTEM ${GLOG_INCLUDE_DIRS})
list(APPEND radiation_LIBRARIES ${GLOG_LIBRARIES})

# Find Gurobi. Right now Travis is only Linux machine and it
# does not have Gurobi, so this will filter out Travis.
# FIX THIS LATER!
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  include("cmake/Modules/FindGurobi.cmake")
  include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
  list(APPEND radiation_LIBRARIES ${GUROBI_LIBRARIES})
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Linux")

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

# Find Gurobi.
include("cmake/Modules/FindGurobi.cmake")
include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
list(APPEND radiation_LIBRARIES ${GUROBI_LIBRARIES})

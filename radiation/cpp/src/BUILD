 cc_library(
    name = "main",
    srcs = glob(["*.cc"]),
    hdrs = glob(["*.h"]),
    copts = [
        "-Iexternal/gtest/include",
        "-Iexternal/eigen/Eigen/src"
    ],
    deps = [
        #"@glog//:main",
        "@eigen//:main",
    ],
    visibility = ["//test:__pkg__"],
)
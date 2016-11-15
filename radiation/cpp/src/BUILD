 cc_library(
    name = "radiation",
    srcs = glob(["*.cc])
    hdrs = glob(["*.h"])
    copts = [
        "-Iexternal/gtest/include",
        "-Iexternal/glog/src/windows/glog",
        "-Iexternal/eigen/Eigen/src"
    ],
    deps = [
        "@glog//:main",
        "@eigen//:main",
    ]
    visibility = ["//visibility:public"],
)
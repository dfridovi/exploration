cc_test(
    name = "tests",
    srcs = glob(["*.cc"]),
    copts = [
        "-Iexternal/gtest/include",
        "-Iexternal/glog/src/windows/glog",
        "-Iexternal/eigen/Eigen/src"
    ],
    deps = ["@gtest//:main"],
    size = "small",
)
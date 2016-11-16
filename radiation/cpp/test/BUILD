cc_test(
    name = "tests",
    srcs = glob(["*.cc"]),
    copts = [
        "-Iexternal/gtest/include",
    ],
    deps = [
        "@gtest//:main",
        "//src:main",
    ],
    size = "small",
)
cc_library(
    name = "main",
    srcs = glob(
        ["src/**/*.cc"],
    ),
    hdrs = glob([
        "src/**/*.h",
    ]),
    copts = ["-Iexternal/glog/include"],
    linkopts = ["-pthread"],
    visibility = ["//visibility:public"],
)
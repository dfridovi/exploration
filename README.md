# exploration

[![Build Status](https://travis-ci.org/dfridovi/exploration.svg?branch=master)](https://travis-ci.org/dfridovi/exploration)
[![License](https://img.shields.io/badge/license-BSD-blue.svg)](https://github.com/dfridovi/exploration/blob/master/LICENSE)

Simulations of autonomous robot navigation, using information-theoretic metrics. For a detailed description of this work, please review the document [here](http://people.eecs.berkeley.edu/~dfk/pdfs/ee290s_writeup.pdf), which gives a brief overview of related work and presents some interesting theoretical results.

## Status
Although we hope to build **exploration** into a full-fledged, generic information-maximization solver for applications in robotics, we have decided to begin by tackling a toy problem inspired by _radiation detection_. 

## Structure
The _radiation detection_ problem is contained in the `radiation/` directory, and is subdivided into Python and C++ implementations. In the Python subdirectory (`radiation/python/`), the source code is readily visible, unit tests are located in the `radiation/python/test/` directory, and scripts are located in the `radiation/python/scripts/` directory. In the C++ implementation (`radiation/cpp/`), source code and header files are located in `radiation/cpp/src/` and `radiation/cpp/include/`, respectively. Unit tests are in `radiation/cpp/test/`, and executables are in `radiation/cpp/exec/`.

## Dependencies
I may miss a few here, but for the Python implementation the dependencies are as follows:

* nose (unit testing)
* scipy (for optimization)
* numpy (for general numerics)

And for the C++ implementation:

* [Google Ceres](http://ceres-solver.org) (_fast_ nonlinear least squares solver)
* [Gurobi](http://www.gurobi.com) (general commercial convex optimization solver)
* [Eigen3](http://eigen.tuxfamily.org/dox/) (header-only linear algebra library)
* Gflags (Google's command-line flag manager)
* Glog (Google's logging tool)

All of these may be installed very easily. If you run into any trouble, though, we are more than happy to help you figure out what's going on. Just post an [issue](https://github.com/dfridovi/exploration/issues) on this repository and we will reply as soon as possible.

## Usage
### Python implementation
Begin by adding the directory `exploration/radiation/python/` to your PYTHONPATH, e.g. by appending the following line to your `~/.bashrc` (Linux) or `~/.bash_profile` (Mac) or (sorry, Windows):

```
export PYTHONPATH=$PYTHONPATH:(path-to-this-repo)/exploration/radiation/python
```

To run unit tests, you'll need to invoke the following command from `exploration/radiation/python/`:

```
nosetests -v
```

And to run any of the scripts within `exploration/radiation/python/scripts/`, you'll need to run the command `python (path-to-script)/(name-of-script)`, e.g. to run the LP-based explorer, you can run

```
python run_explorer_lp.py
```

from within the `scripts/` directory.

### C++ implementation
You'll need to begin by building the repository. From the top directory (`exploration/radiation/cpp/`), type the following sequence of commands:

```
mkdir bin
mkdir build
cd build
cmake ..
make -j4
```

This should build all tests and executables. In order to run tests, you can run the following command:

```
./run_tests
```

from within the `build/` directory you just made. All the tests should pass, and none should take more than a second or so to run.

Executables are automatically placed within the `bin/` directory that you created. To run them, just type `./(name-of-executable)`, e.g. to run the LP explorer:

```
./run_explorer_lp
```

To the extent that it makes sense, all parameters are accessible from the command line via Gflags. So, for example, in order to run the LP explorer for 100 iterations instead of the default 10, you could run the following:

```
./run_explorer_lp --num_iterations=100
```

## API documentation
We use Doxygen to auto-generate web-based [documentation](https://dfridovi.github.io/exploration/radiation/cpp/documentation/html/). Although we do not follow the Doxygen guidelines for writing comments, auto-generation still seems to do a fairly reasonable job.

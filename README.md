# RObustly Complete control Synthesis (ROCS)

This is an improved version of **ROCS** the robust symbolic control synthesis tool for nonlinear systems, which provides two control synthesis methods from linear temporal logic (LTL) formulas:
- Abstraction-based control synthesis: state and control space are uniformly discretized.
- Specification-guided control synthesis: state space are partitioned adaptively with respect to given specifications and local dynamics of the systems.



## New features
- Supports the standard **abstraction-based control synthesis**.
- Supports control synthesis for sampled-data systems by incoporating Taylor models of ODE solutions. No analytical bounds, e.g., growth functions, are needed.
- Supports goal and avoiding area representation in terms of both intervals and nonlinear inequilities.



## Installation
### Prerequisite packages and softwares
 - [Boost c++ libraries](http://www.boost.org) >= 1.69.0: libraries for graph searching, ode solvers;
 - [Armadillo](http://arma.sourceforge.net): a c++ linear algebra library;
 - matlab: control simulation and graphics.
 
### Configuration
 - Update the environment variable DYLD_LIBRARY_PATH in terminal:
```
export DYLD_LIBRARY_PATH=“/Applications/MATLAB_R2016a.app/bin/maci64:/Applications/MATLAB_R2016a.app/sys/os/maci64:$DYLD_LIBRARY_PATH”
```


## Usage
- Clone this repositary. Source code will be downloaded and ready for use.
- Make sure the prerequisites are installed.
- Prepare the configuration file for build:
  + The default configuration file *cflags* in the root path is for Mac OSX, the user only need to modify the paths for the external packages, e.g.,
  ```
  BSTDIR := /usr/local/opt/boost
  MATDIR := /Applications/MATLAB_R2016a.app
  LADIR := /usr/local/opt/armadillo
  ```
  + for Linux users, rename the file *cflags-linux* to *cflags*. If all the prerequisites are installed in the default paths, only modify the matlab path, e.g., 
  ```
  MATDIR = /usr/local/MATLAB/R2016b
  ```
- To run examples, go to the folder *examples*:
  + To build all the examples: type `$ make examples` in terminal.
  + To build one of the examples: go to the corresponding subfolder and type `$ make` in terminal.
- Execute compiled examples by typing the name of the generated executable, e.g. `./dcdc`. **Data will be generated in each subfolder corresponding to the specific example**.
- Run *simulation.m* in the same example folder using matlab.



## Examples
Simulations of the examples can be found under the folder *examples*:
- *examples/dcdc*: DCDC converter invariance control.
- *examples/car*: Control a unicycle to a goal area while avoiding obstacles as well as automatic parallel parking.
- *examples/ipdl*: Regulate an inverted pendulum to the upright position.
- *examples/temp*: Control the room temperature (4-mode system) to a desired temperature (a setpoint) and keep the temperature around the setpoint.
- *examples/vdp*: Estimation of the region of attraction (ROA) for Van der Pol equations.



## Create your own examples
To create your own examples, you have to provide a *main file* specifying
- the original discrete-time or continuous-time dynamics, and
- all necessary data for problem settings, including the state space X, the set of control values U, sampling time, and target area and so on.

For the DCDC converter example (go to the folder *examples/dcdc*), the user can find a main file *dcdc.cpp*.



## Features under developement
- Control synthesis strategies to deal with richer fragments of LTL, e.g., GR(1), co-safe LTL (or scLTL).



## Note
Tested on
- Mac OSX 10.11/12.
- Ubuntu 16.04.3 LTS.



## References
[1][Y. Li, J. Liu, ROCS: A Robustly Complete Control Synthesis Tool for Nonlinear Dynamical Systems. *Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control (part of CPS Week) - HSCC '18*.](http://dl.acm.org/citation.cfm?doid=3178126.3178153)


## Contact
Email: <yinan.li@uwaterloo.ca>
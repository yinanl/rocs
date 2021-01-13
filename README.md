# RObustly Complete control Synthesis (ROCS)

This is an improved version of **ROCS** [1] the robust symbolic control synthesis tool for nonlinear systems, which provides two control synthesis methods from linear temporal logic (LTL) formulas:
- Abstraction-based control synthesis: state and control space are uniformly discretized.
- Specification-guided control synthesis: state space are partitioned adaptively with respect to given specifications and local dynamics of the systems.



## New features
- Supports control synthesis with respect to a general class of LTL formulas that can be translated to deterministic B\"uchi automata (DBA) or **DBA-recognizable** formulas.
- Supports the standard **abstraction-based control synthesis**.
- Supports control synthesis for sampled-data systems by incoporating Taylor models of ODE solutions. No analytical bounds, e.g., growth functions, are needed.
- Provide python simulation files for the examples and python simulation interfaces in the folder *python*.
- Supports goal and avoiding area representation in terms of both intervals and nonlinear inequilities.



## Installation
### Prerequisite packages and softwares
 - [Boost c++ libraries](http://www.boost.org) >= 1.69.0: libraries for graph searching, ode solvers;
 - [Armadillo](http://arma.sourceforge.net) (tested on 9.900.2 and 9.850.1): a c++ linear algebra library;
 - [HDF5](https://www.hdfgroup.org/downloads/hdf5/) (1.10.0): a high-performance data management and storage suite.
 - MATLAB: control simulation and graphics. **Installation of MATLAB is no longer mandatory for saving control synthesis results. Data saved in HDF5 format can also be loaded in MATLAB for simulation.** Update the environment variable DYLD_LIBRARY_PATH in terminal if you have MATLAB installed and wish to save data in MATLAB data format:
```
export DYLD_LIBRARY_PATH=“Your_MATLAB_install_path/bin/maci64:Your_MATLAB_install_path/sys/os/maci64:$DYLD_LIBRARY_PATH”
```


## Usage
- Clone this repositary. Source code will be downloaded and ready for use.
- Make sure the prerequisites are installed.
- Prepare the configuration file for build:
  + The default configuration file *cflags* in the root path is for Mac OSX, the user only need to modify the paths for the external packages, e.g.,
  ```
  BSTDIR := /usr/local/opt/boost
  LADIR := /usr/local/opt/armadillo
  ```
  If you choose to save control synthesis results in MATLAB data format, also edit:
  ```
  MATDIR := Your_MATLAB_install_path
  ```
  
  + for Linux users, rename the file *cflags-linux* to *cflags*. Do not need to modify *cflags* if all the prerequisites are installed in the default paths. Same as in Mac OSX, you may need to modify the matlab path, e.g., 
  ```
  MATDIR = /usr/local/MATLAB/XXX(your version)
  ```
  
- To run examples, go to the folder *examples* and the subfolder corresponding to the specific case. Then type `$ make XXX` in terminal, where `XXX` is the target executable. There is a target executable for performing every control synthesis problem, and there are several control synthesis problems provided for some of the dynamical systems. Check in `Makefile` for target names.
- Execute compiled examples by typing the name of the generated executable, e.g. `./dcdcInv` for the invariance control of DCDC converter example. **Data will be generated in each subfolder corresponding to the specific example**.
- To check the control synthesis results, run MATLAB or Python simulation files in the same example folder using matlab.



## Examples
Simulations of the examples can be found under the folder *examples*:
- *examples/aircraft*: reachability control of an aircraft longitudinal model for safe landing [2].
- *examples/dcdc*: DCDC converter invariance control [3].
- *examples/car*: control a unicycle to a goal area while avoiding obstacles as well as automatic parallel parking [4].
- *examples/integrator*: reach control of a double integrator.
- *examples/locomotion*: mode switching for a bipedal robot (see instructions in the related folder).
- *examples/ipdl*: regulate an inverted pendulum to the upright position [3].
- *examples/temp*: control the room temperature (4-mode system) to a desired temperature (a setpoint) and keep the temperature around the setpoint [5].
- *examples/vdp*: estimation of the region of attraction (ROA) for Van der Pol equations [6].
- *examples/Moore-Greitzer*: control of operation points of Moore-Greitzer jet engine.
- *examples/scara*: motion planning of a two-link SCARA manipulator.



## Create your own examples
To create your own examples, you have to provide
- a *main file* specifying
  + the original discrete-time or continuous-time dynamics, and
  + all necessary data for problem settings, including the state space X, the set of control values U, sampling time, and target area and so on.

For the DCDC converter example (go to the folder *examples/dcdc*), the user can find a main file *dcdc.cpp*.

- a specification file if you wish to perform DBA control synthesis. Examples of such files can be found in some of the examples such as *examples/car*.



## Features under developement
- Control synthesis strategies to deal with $`\omega`$-regular languages.
- Interface between LTL to DRA translator and ROCS
- Parallelization



## Note
Tested on
- Mac OSX 10.11/12.
- Ubuntu 16.04.3 LTS.



## References
[1][Y. Li, J. Liu, ROCS: A Robustly Complete Control Synthesis Tool for Nonlinear Dynamical Systems. *Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control (part of CPS Week) - HSCC '18*.](http://dl.acm.org/citation.cfm?doid=3178126.3178153)

[2][Reissig, G., Weber, A., and Rungger, M.(2017). Feed-back refinement relations for the synthesis of symboliccontrollers. IEEE Trans. Automat. Contr., 62(4), 1781–1796.](http://ieeexplore.ieee.org/document/7519063/)

[3][Y. Li, J. Liu, Invariance Control Synthesis for Switched Nonlinear Systems: An Interval Analysis Approach. IEEE Trans. Automat. Contr., 63(7), 2206-2211.](http://ieeexplore.ieee.org/document/8062786/)

[4][Y. Li, J. Liu, Robustly Complete Synthesis of Memoryless Controllers for Nonlinear Systems with Reach-and-Stay Specifications, IEEE Trans. Automat. Contr.](https://ieeexplore.ieee.org/document/9067073)

[5][J. Liu, N. Ozay, U. Topcu, R.M. Murray, Synthesis of Reactive Switching Protocols From Temporal Logic Specifications. IEEE Trans. Automat. Contr., 58(7), 1771-1785.](http://ieeexplore.ieee.org/document/6457409/)

[6][Y. Li, J. Liu, Robustly Complete Synthesis of Memoryless Controllers for Nonlinear Systems with Reach-and-Stay Specifications. ArXiv: 1802.09082v2](http://arxiv.org/abs/1802.09082v2)


## Contact
Email: <yinan.li@uwaterloo.ca>
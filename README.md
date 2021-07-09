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



## System requirements

- Operating system:
  - Mac OSX (tested on >=10.11), or
  - Ubuntu (tested on 20.04 LTS), or
  - Manjaro (tested on 5.9.16).
- Compilers: gcc (clang), g++ (clang++)



## Installation
### Prerequisite packages and softwares

Note: To simplify the editing of the configuration file *cflags*, we suggest install **Boost**, **Armadillo** and **HDF5** by using *Homebrew* (for Mac OSX) or *apt* (for Ubuntu).

 - [Boost c++ libraries](http://www.boost.org) >= 1.69.0: libraries for graph searching, ode solvers;

 - [Armadillo](http://arma.sourceforge.net) (tested on 9.900.2, 9.850.1, and 10.2.2): a c++ linear algebra library;

 - [HDF5](https://www.hdfgroup.org/downloads/hdf5/) (1.10.0/1.12.0): a high-performance data management and storage suite. Install the package from source code may take a while.

      + For Mac OSX users, using [Homebrew](https://formulae.brew.sh/formula/hdf5) for installation is stongly recommended.
      + For Ubuntu users, it is also suggested to install the [package](https://launchpad.net/ubuntu/focal/+source/hdf5) by running
      ```
      sudo apt install libhdf5-dev
      ```
      + Manjaro users can also use the [package](https://archlinux.org/packages/community/x86_64/hdf5)

 - MATLAB: control simulation and graphics. **Installation of MATLAB is no longer mandatory for saving control synthesis results. Data saved in HDF5 format can also be loaded in MATLAB for simulation.** Update the environment variable `DYLD_LIBRARY_PATH` in terminal if you have MATLAB installed and wish to save data in MATLAB data format:
```
export DYLD_LIBRARY_PATH=“Your_MATLAB_install_path/bin/maci64:Your_MATLAB_install_path/sys/os/maci64:$DYLD_LIBRARY_PATH”
```
 - Python 3: an alternative of MATLAB. Make sure the packages **scipy**, **numpy**, **matplotlib**, and **h5py** are installed. The following command can be used:

   ```
   pip3 install scipy numpy matplotlib h5py
   ```



### Configuration of makefiles

In the configuration file *cflags* in the root path, the paths assigned to `INCS` and `LDFLAGS` might need to be changed accordingly.

- For Mac OSX users, modify the path for the external packages to your actual install path, e.g.,
  ```
  EXTDIR := Your_package_install_path
  ```
  If the prerequisites are installed in different paths, make sure the paths are given correctly to `INCS` and `LDFLAGS`.
  
- For Linux users, rename the file *cflags-linux* to *cflags*.
  * On Ubuntu, if HDF5 is installed via `apt` and `mpi` is not available, HDF5 may be installed only for serial compuation. In this case, the following changes may apply to the *cflags*:
    + Add to `INCS` the following path:
    ```
    -I$(INCDIR)/hdf5/serial
    ```
    + Change `-lhdf5_hl` and `-lhdf5` to `-lhdf5_serial_hl` and `-lhdf5_serial`, respectively.
  * If HDF5 is installed by compiling the source code from HDF5 official website, make sure that *zlib* and *szip* are installed as in their instruction (*release_docs/INSTALL*). Use
  ```
  ./configure --prefix=Your_install_path --enable-cxx --with-szip=Path_to_szip
  ```
  to configure *cmake*. After installation, add to `INCS` the following paths:
  ```
  -I$(Your_install_path)/include -I$(Path_to_szip)/include
  ```
  and add to `LDFLAGS` the following paths:
  ```
  -I$(Your_install_path)/lib -I$(Path_to_szip)/lib
  ```
  
- For all users, if you choose to save control synthesis results in MATLAB data format, make sure that MATLAB is installed and edit:
  ```
  MATDIR := Your_MATLAB_install_path
  ```
  Otherwise, **comment out or delete** the line starting with `MATDIR`.



## Usage

- Clone this repositary. Source code will be downloaded and ready for use.
- Make sure the prerequisites are installed.
- Prepare the configuration file *cflags* for build according to the installation instruction.
- To run examples, go to the folder *examples* and the subfolder corresponding to the specific case. Then type `$ make XXX` in terminal, where `XXX` is the target executable. There is a target executable for performing every control synthesis problem, and there are several control synthesis problems provided for some of the dynamical systems. Check in `Makefile` for target names.
- Execute compiled examples by typing the name of the generated executable, e.g. `./dcdcInv` for the invariance control of DCDC converter example. **Data will be generated in each subfolder corresponding to the specific example**.
- To check the control synthesis results, run MATLAB or Python simulation files in the same example folder.



## Examples
Simulations of the examples can be found under the folder *examples*:
- *examples/aircraft*: reachability control of an aircraft longitudinal model for safe landing [3].
- *examples/dcdc*: DCDC converter invariance control [4].
- *examples/car*: control a unicycle to a goal area while avoiding obstacles as well as automatic parallel parking [5].
- *examples/integrator*: reach control of a double integrator.
- *examples/locomotion*: mode switching for a bipedal robot (see instructions in the related folder).
- *examples/ipdl*: regulate an inverted pendulum to the upright position [4].
- *examples/temp*: control the room temperature (4-mode system) to a desired temperature (a setpoint) and keep the temperature around the setpoint [6].
- *examples/vdp*: estimation of the region of attraction (ROA) for Van der Pol equations [7].
- *examples/Moore-Greitzer*: control of operation points of Moore-Greitzer jet engine [8].
- *examples/scara*: motion planning of a two-link SCARA manipulator [8].



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



## References
[1][Y. Li, J. Liu, ROCS: A Robustly Complete Control Synthesis Tool for Nonlinear Dynamical Systems. *Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control (part of CPS Week) - HSCC '18*.](http://dl.acm.org/citation.cfm?doid=3178126.3178153)

[2]Y. Li, Z. Sun, and J. Liu, ROCS 2.0: An Integrated Temporal Logic Control Synthesis Tool for Nonlinear Dynamical Systems. *7th IFAC Conference on Analysis and Design of Hybrid Systems (ADHS 2021)*.

[3][Reissig, G., Weber, A., and Rungger, M.(2017). Feed-back refinement relations for the synthesis of symboliccontrollers. IEEE Trans. Automat. Contr., 62(4), 1781–1796.](http://ieeexplore.ieee.org/document/7519063/)

[4][Y. Li, J. Liu (2018). Invariance Control Synthesis for Switched Nonlinear Systems: An Interval Analysis Approach. IEEE Trans. Automat. Contr., 63(7), 2206-2211.](http://ieeexplore.ieee.org/document/8062786/)

[5][Y. Li, J. Liu (2021). Robustly Complete Synthesis of Memoryless Controllers for Nonlinear Systems with Reach-and-Stay Specifications, IEEE Trans. Automat. Contr., 66(3), 1199-1206](https://ieeexplore.ieee.org/document/9067073)

[6][J. Liu, N. Ozay, U. Topcu, R.M. Murray (2013). Synthesis of Reactive Switching Protocols From Temporal Logic Specifications. IEEE Trans. Automat. Contr., 58(7), 1771-1785.](http://ieeexplore.ieee.org/document/6457409/)

[7][Y. Li, J. Liu (2018). Robustly Complete Synthesis of Memoryless Controllers for Nonlinear Systems with Reach-and-Stay Specifications. ArXiv: 1802.09082v2](http://arxiv.org/abs/1802.09082v2)

[8][Y. Li, Z. Sun, and J. Liu (2021). A Specification-Guided Framework for Temporal Logic Control of Nonlinear Systems. ArXiv:2104.01385](https://arxiv.org/abs/2104.01385)


## Contact
Email: <yinan.li@uwaterloo.ca>
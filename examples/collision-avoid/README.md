- The number in file names, e.g. "safety-x1-x2-x3-x4.eps", "controller_safety_x1-x2-x3-x4.h5":
  + x1: disturbance range
  + x2: least safety distance
  + x3: sampling time
  + x4: state space discretisation precision

- Compile and execute static planner
  + `make plan`
  + `./plan $(SPEC)` (default precision 0.2 0.2 0.2) or `./plan $(SPEC) $(eta[0] eta[1] eta[2])`

- Compile and execute local collision avoidance
  + `make safesetAbst`
  + `/.safesetAbst` or `/.safesetAbst $(samplingtime) $(eta_r[0] eta_r[1] eta_r[2])`

- Simulate online collision avoidance
  + With RH module only: `make sim`
  ```
  ./sim $(dbafile ctlrfile cafile graphfile eta[0] eta_r[0])
  ```
  + With Detouring: `make simRe`
  ```
  ./simRe $(dbafile ctlrfile cafile graphfile eta[0] eta_r[0])
  ```
  + Handle multiple obstacles: `make simMulti`
  ```
  ./simMulti $(dbafile ctlrfile cafile graphfile eta[0] eta_r[0])
  ```

- Visualize simulation
  + `python replay.py $(case) $(simulationtime)`
  + `python replay_multi.py $(simulationtime)`
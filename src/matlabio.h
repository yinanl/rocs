/**
 *  matlabio.h
 *
 *  Input/output classes to matlab files.
 *
 *  Created by Yinan Li on May 10, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _matlabio_h
#define _matlabio_h

#include <iostream>

#include <mat.h>
#include <matrix.h>

#include "grid.h"
#include "csolver.h"
#include "transition.h"
#include "dsolver.h"


namespace rocs {

  /**
   * A mat file writer class.
   */
  class matWriter {
  private:
    const char* _filename;
    MATFile *_pmat;

  public:
    matWriter(const char* f) : _filename(f) {}

    /**
     * Open a mat file.
     **/
    bool open() {
      _pmat = matOpen(_filename, "w");
      return _pmat == NULL;
    }
    bool openfile(const char* f) {
      _filename = f;
      return open();
    }

    /**
     * Close a mat file.
     **/
    bool close() {return EOF == matClose(_pmat);}


    /**
     * Write a real number to a mat file.
     * @param x the real number.
     * @param varname the variable name to be saved as.
     **/
    void write_real_number(const double x, const char *varname);

    /**
     * Write an array of real numbers to a mat file.
     * @param x the array of real numbers.
     * @param varname the variable name to be saved as.
     **/
    template<typename T>
    void write_real_array(const std::vector<T> x, const char *varname);

    /**
     * Write the state space to a mat file.
     * @param ws the system state space (a ivec).
     * @param varname the variable name to be saved as.
     **/
    void write_state_space(const ivec &ws, const char *varname);

    /**
     * Write the input values to a mat file.
     * @param ugrid the set of input values (a grid).
     * @param varname the variable name to be saved as.
     **/
    void write_input_values(const grid &ugrid, const char *varname);

    /**
     * Write the goal area to a mat file.
     * @param goal the goal area (in the form of ivec)
     * @param varname the variable name to be saved as.
     **/
    void write_goal_area(const ivec &goal, const char *varname);

    /**
     * Write the areas-to-avoid to a mat file.
     * @param obs the areas-to-avoid (a vector of ivecs)
     * @param varname the variable name to be saved as.
     **/
    void write_unsafe_area(const std::vector<ivec> &obs, const char *varname);

    /**
     * Write leaves under a node to a .mat file (Depth-first).
     * @param ctlr the controller in the form of SPtree.
     * @param ptrn pointer to the given node.
     */
    void write_sptree_leaves(const SPtree &ctlr, SPnode *ptrn);

    /**
     * Extract the boundary of a given grid.
     *
     * @param g the grid data.
     * @param win the list indicating inside/outside of the target region.
     * @param varname the variable name to be saved as.
     **/
    void write_boundary(grid &g, const std::vector<bool> &win, const char *varname);
    
    /**
     * Extract the boundary of the winning set.
     *
     * Discretize the state space and collect the ones lie on the boundary.
     *
     * @param eta the grid size of the state space.
     * @param sol the solver which contains the spec info.
     * @param varname the variable name to be saved as.
     **/
    void write_winset_boundary(const double eta[], const CSolver &sol, const char *varname);

    /**
     * Write the problem setting into a mat file.
     * @param sys the defined control system (use template type S).
     * @param sol the solver that contains the spec info.
     **/
    template<typename S>
      void write_problem_setting(const S &sys, const CSolver &sol) {
      if (_pmat == NULL) {open();}
      write_real_number(sys._tau, "ts");
      write_state_space(sys._workspace, "X");
      write_input_values(sys._ugrid, "U");
      write_goal_area(sol._goal, "G");
      write_unsafe_area(sol._obs, "xobs");
    }

    /**
     * Write the controller into a mat file.
     *
     * The controller is saved in 3 variables:
     * + pavings: the interval boxes,
     * + ctlr: a list of binary values indicating the admissible control inputs, and
     * + tags: the tags (0,1,or 2) of the boxes.
     *
     * @param sol the solver which contains the controller info.
     **/
    void write_sptree_controller(const CSolver &sol) {
	write_sptree_leaves(sol._ctlr, sol._ctlr._root);
    }

    /**
     * Write the controller into a mat file as arrays (keep tree structure).
     *
     * - ctree: index-searching array (2^H-1 x [split axis(1), intervals(2dim)]).
     * - cindex: comparison array (# of nodes x intervals(2dim)).
     * - cvalue: control array (# of leaves x [index(1), # of control values]).
     *
     * @param sol the solver which contains the controller info.
     */
    void write_controller_serialized(const CSolver &sol);


    /**
     * Save all transitions to a .mat file.
     *
     * The file is rewritten each time.
     * @param trans the transitions to be saved: 
     * trans_post, post, postptr,
     * trans_pre, pre, preptr
     */
    void write_transitions(const fts &trans,
			   const char* vtranspost, const char* vpost, const char* vptrpost,
			   const char* vtranspre, const char* vpre, const char* vptrpre);


    /**
     * Save uniform grids to a .mat file
     *
     * @param varname 
     */
    void write_uniform_grids(const grid &g, const char* varname);


    /**
     * Write the discrete controller into a mat file.
     *
     * The controller is saved in 2 variables:
     * + optctlr: the optimal controller (a vector of control indices).
     * + leastctlr: a list of binary values indicating the admissible control inputs.
     *
     * @param sol the solver which contains the controller info.
     **/
    void write_discrete_controller(const DSolver &dsol,
				   const char* vleastctlr, const char* voptctlr);
    
  };//matWriter class

  
  template<typename T>
    void matWriter::write_real_array(const std::vector<T> x, const char *varname) {
    if (_pmat == NULL) {open();}
    mxArray *mx = mxCreateDoubleMatrix(x.size(), 1, mxREAL);
    double *ptrx = mxGetPr(mx);
    for (int i = 0; i < x.size(); ++i) {
      ptrx[i] = double(x[i]);
    }
    if (matPutVariable(_pmat, varname, mx)!=0) {
      std::cout << "matWriter::write_real_array: Error writing the real number.\n";
    }
    mxDestroyArray(mx);
  }//matWriter::write_real_array


  /**
   * A mat file reader class.
   */
  class matReader {
  private:
    const char* _filename;
    MATFile *_pmat;

  public:
  matReader(const char* f) : _filename(f) {}

    /**
     * Open a mat file.
     **/
    bool open() {
      _pmat = matOpen(_filename, "r");
      return _pmat == NULL;
    }
    bool openfile(const char* f) {
      _filename = f;
      return open();
    }

    /**
     * Close a mat file.
     **/
    bool close() {return EOF == matClose(_pmat);}

    /**
     * Read a transition system from .mat file.
     * @param trans a transition system to be constructed from the file.
     */
    void read_transitions(fts &trans);

  };


}//namespace rocs

#endif

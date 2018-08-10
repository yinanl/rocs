/**
 *  grid.h
 *
 *  A uniform grid class.
 *
 *  Created by Yinan Li on Jan. 21, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _grid_h_
#define _grid_h_


#include <cstdlib>

#include <cassert>
#include <iostream>
#include <fstream>

#include <cmath>
#include <vector>

#include "interval_vector.h"


namespace rocs {
  
class grid {

 public:

  int _dim; /**< the dimension */

  ivec _bds; /**< an interval of the space */

  std::vector<double> _gw;  /**< grid width */

  size_t _nv; /**< overall number of grids */

  std::vector<size_t> _size; /**< number of grids in each dimension */

  std::vector<double> _valmin; /**< minimal grid center in each dimension */

  std::vector< std::vector<double> > _data; /**< grid data */


  /**
   * A default constructor.
   */
  grid() {}
  /**
   * Construct from lower and upper bound arrays.
   * @param n dimension
   * @param eta[n] an array of grid size.
   * @param lb[n] an array of lower bounds.
   * @param ub[n] an array of upper bounds.
   */
  grid(const int n, const double eta[], const double lb[], const double ub[]) :
  _dim(n), _bds(n), _gw(eta, eta + n), _nv(1), _size(n), _valmin(n) {

    for (int i = 0; i < _dim; ++i) {

      _bds[i] = interval(lb[i], ub[i]);

      _size[i] = floor((ub[i] - lb[i]) / _gw[i]) + 1;

      _nv *= _size[i];

      _valmin[i] = (ub[i] - lb[i] - floor((ub[i]-lb[i])/eta[i]) * eta[i])/2. + lb[i];
      
    }

  }
  /**
   * Construct from an interval vector.
   * @param n dimension
   * @param eta[n] an array of grid size.
   * @param x an interval vector.
   */
  grid(const int n, const double eta[], const ivec &x) :
  _dim(n), _bds(x), _gw(eta, eta + n), _nv(1), _size(n), _valmin(n) {

    std::vector<double> lb = x.getinf();
    std::vector<double> ub = x.getsup();
    
    for (int i = 0; i < _dim; ++i) {

	_size[i] = floor((ub[i] - lb[i]) / _gw[i]) + 1;

	_nv *= _size[i];

	_valmin[i] = ((ub[i] - lb[i]) - floor((ub[i] - lb[i]) / eta[i]) * eta[i])/2. + lb[i];
    }

  }

  
  /**
   * Initialization from bound arrays.
   * @see grid()
   */
  void init(const int n, const double eta[], const double lb[], const double ub[]);
  /**
   * Initialization from an interval vector.
   * @see grid().
   */
  void init(const int n, const double eta[], const ivec &x);


  /**
   * \brief Perform gridding.
   *
   * Generate center point of the grids by the formula:
   *
   * \f$i= i_1 + N_1 i_2 + N_1 N_2 i_3 + \cdots + N_1\cdots N_{n-1} i_n\f$
   */
  void gridding();
  
  /**
   * Perform gridding from bound arrays without initialization aprioi.
   * @see grid().
   */
  void gridding(const int n, const double eta[],
		const double lb[], const double ub[]);

  /**
   * Perform gridding an interval vector without initialization aprioi.
   * @see grid().
   */
  void gridding(const int n, const double eta[], const ivec &x);

  /**
   * Subgrid a grided interval.
   * @param xc grid centers.
   * @param rp relative portion of subgridding, e.g. 0.2.
   * @return a list of subgrid centers.
   */
  std::vector<std::vector<double> >
    subgridding(std::vector<double> &xc, std::vector<double> &rp);

  /**
   * Perform the actual job of gridding.
   * @param[inout] data the grid centers to be output.
   * @param[in] xmin the list of minimum values in each dimension.
   * @param[in] number the list of numbers of sub grid points in each dimension.
   * @param[in] nv the overall number of sub grid points.
   */
  void griddingHelper(std::vector<std::vector<double> > &data,
		      const std::vector<double> &xmin,
		      const std::vector<double> &gw,
		      const std::vector<size_t> number,
		      const size_t nv);
  

  /**
   * \brief Get the ID number in the grid for a real state.
   *
   * \verbatim val[k] --> ik \endverbatim
   * \f$i = i_1 + N_1 i_2 + \cdots + N_1 N_2\cdots N_{n-1} i_n\f$
   * @param val the real state.
   * @return id the corresponding grid ID.
   */
  size_t val_to_id(std::vector<double> val);

  /**
   * \brief Get the grids covered by an interval area.
   *
   * \f$i = i_1 + N_1 i_2 + \cdots + N_1 N_2 \cdots N_{n-1} i_n \f$
   * @param box the interval area.
   * @param boxin 0: allow box out of range, 1: don't allow.
   * @param strictin 1: collect the grids that entirely inside the box,
   *             0: collect the grids that intersect the box.
   * @return a list of grid IDs.
   */
  std::vector<size_t> subset(ivec &box, bool boxin, bool strictin);

  /**
   * \brief Get the neighbour grids of the given grid point.
   *
   * \f$i = i_1 + N_1 i_2 + \cdots + N_1 N_2 \cdots N_{n-1} i_n \f$
   * @param id the id of the given grid point.
   * @return a list of neighbour grids ID.
   */
  std::vector<size_t> neighbours(size_t id);
  

  /**
   * Print grid information.
   */
  void print_info();
  
  /**
   * Print grid centers. 
   */
  void print_data();

};

} // namespace rocs

#endif

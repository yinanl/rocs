/**
 *  grid.h
 *  A uniform grid class.
 *
 *  Created by Yinan Li on Jan. 21, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _grid_h_
#define _grid_h_


#include <cstdlib>
#include <exception>
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
	std::vector<size_t> _base; /**< _base[k]=_size[k-1]*_size[k-2]*...*1 */
	std::vector<double> _valmin; /**< minimal grid center in each dimension */
	std::vector< std::vector<double> > _data; /**< grid data */

	
	/**
	 * No default constructor.
	 */
	grid() = delete;
	/**
	 * Construct from lower and upper bound arrays.
	 * @param n dimension
	 */
	grid(const int n) : _dim(n), _bds(n), _gw(n), _nv(1), _size(n), _base(n), _valmin(n) {}
	/**
	 * Construct from lower and upper bound arrays.
	 * @param n dimension
	 * @param eta[n] an array of grid size.
	 * @param lb[n] an array of lower bounds.
	 * @param ub[n] an array of upper bounds.
	 */
	grid(const int n, const double eta[], const double lb[], const double ub[]) :
	    _dim(n), _bds(n) , _gw(n), _nv(1), _size(n), _base(n), _valmin(n) {
	    
	    init(eta, lb, ub);
	}
	/**
	 * Construct from an interval vector.
	 * @param n dimension
	 * @param eta[n] an array of grid size.
	 * @param x an interval vector.
	 */
	grid(const int n, const double eta[], const ivec &x) :
	    _dim(n), _bds(x), _gw(n), _nv(1), _size(n), _base(n), _valmin(n) {
	    
	    init(eta, x);
	}

  
	/**
	 * Initialize all member variables from bound arrays.
	 * @param eta an array of grid width.
	 * @param lb an array of lower bound of the domain.
	 * @param ub an array of upper bound of the domain.
	 */
	// void init(const double eta[], const double lb[], const double ub[]);
	template<typename T>
	void init(const double eta[], const T lb, const T ub) {
	    _gw.assign(eta, eta + _dim);
	    for (int i = 0; i < _dim; ++i)
		_bds[i] = interval(lb[i], ub[i]);
	    size_t b = 1;
	    for (int i = 0; i < _dim; ++i) {
		_size[i] = floor((ub[i] - lb[i]) / eta[i]) + 1;
		_nv *= _size[i];
		_valmin[i] = ((ub[i] - lb[i]) - floor((ub[i] - lb[i]) / eta[i]) * eta[i])/2. + lb[i];
		_base[i] = b;
		b *= _size[i];
	    }
	}
	/**
	 * Initialization from an interval vector.
	 * @see grid().
	 */
	void init(const double eta[], const ivec &x) {
	    std::vector<double> lb(_dim), ub(_dim);
	    x.getinf(lb);
	    x.getsup(ub);
	    init(eta, lb, ub);
	}


	/**
	 * \brief Write center point of the grids to _data by the formula:
	 * \f$i= i_1 + N_1 i_2 + N_1 N_2 i_3 + \cdots + N_1\cdots N_{n-1} i_n\f$
	 */
	void gridding();
	void gridding(const double eta[], const double lb[], const double ub[]); /* used after the constructor grid(n) */
	void gridding(const double eta[], const ivec &x); /* used after the constructor grid(n) */

	/**
	 * Subgrid a grided interval.
	 * @param xc grid centers.
	 * @param rp relative portion of subgridding, e.g. 0.2.
	 * @return a list of subgrid centers.
	 */
	std::vector<std::vector<double> >
	subgridding(std::vector<double> &xc, std::vector<double> &rp);

	/**
	 * Compute the value of every index (id_to_val).
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
	 * \brief Get the center value x of a grid with ID number i.
	 * 
	 * @param[inout] val the real state.
	 * @param[in] id the corresponding grid ID.
	 */
	void id_to_val(std::vector<double> &val, size_t id) const;

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

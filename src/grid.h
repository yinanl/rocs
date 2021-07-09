/**
 *  The header file of a uniform grid class.
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

    /**
     * \brief A uniform grid class.
     */
    class grid {

    public:

	int _dim; /**< The dimension */
	ivec _bds; /**< An interval of the space */
	std::vector<double> _gw;  /**< Grid width */
    
	size_t _nv; /**< Overall number of grids */
	std::vector<size_t> _size; /**< Number of grids in each dimension */
	std::vector<size_t> _base; /**< _base[k]=_size[k-1]*_size[k-2]*...*1 */
	std::vector<double> _valmin; /**< Minimal grid center in each dimension */
	std::vector< std::vector<double> > _data; /**< Grid data */

	
	/**
	 * No default constructor.
	 */
	grid() = delete;
	
	/**
	 * \brief Construct a grid from lower and upper bound arrays.
	 *
	 * @param[in] n The system state dimension
	 */
	grid(const int n) : _dim(n), _bds(n), _gw(n), _nv(1), _size(n), _base(n), _valmin(n) {}
	
	/**
	 * \brief Construct a grid from lower and upper bound arrays.
	 *
	 * @param[in] n The system state dimension
	 * @param[in] eta An array of grid size
	 * @param[in] lb An array of lower bounds
	 * @param[in] ub An array of upper bounds
	 */
	grid(const int n, const double eta[], const double lb[], const double ub[]) :
	    _dim(n), _bds(n) , _gw(n), _nv(1), _size(n), _base(n), _valmin(n) {
	    
	    init(eta, lb, ub);
	}
	
	/**
	 * \brief Construct a grid from an interval vector.
	 *
	 * @param[in] n The system state dimension
	 * @param[in] eta An array of grid size
	 * @param[in] x An interval vector
	 */
	grid(const int n, const double eta[], const ivec &x) :
	    _dim(n), _bds(x), _gw(n), _nv(1), _size(n), _base(n), _valmin(n) {
	    
	    init(eta, x);
	}
  
	/**
	 * \brief Initialize all member variables from arrays of bounds.
	 *
	 * @param[in] eta An array of grid width
	 * @param[in] lb An array of lower bound of the domain
	 * @param[in] ub An array of upper bound of the domain
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
	 * \brief Initialization from an interval vector.
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

	/**
	 * \brief Write center point of the grids to _data.
	 *
	 * Use it only after using the constructor grid(n)
	 */
	void gridding(const double eta[], const double lb[], const double ub[]);

	/**
	 * \brief Write center point of the grids to _data.
	 *
	 * Use it only after using the constructor grid(n)
	 */
	void gridding(const double eta[], const ivec &x); 

	/**
	 * \brief Subgrid a grided interval.
	 *
	 * @param[in] xc The grid centers
	 * @param[in] rp The relative precision w.r.t. a cell
	 * @return a list of subgrid centers.
	 */
	std::vector<std::vector<double> >
	subgridding(std::vector<double> &xc, std::vector<double> &rp);

	/**
	 * \brief Compute the value of every grid index (id_to_val).
	 *
	 * @param[in,out] data The grid centers to be output
	 * @param[in] xmin The list of minimum values in each dimension
	 * @param[in] number The list of numbers of sub grid points in each dimension
	 * @param[in] nv The overall number of sub grid points
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
	 * @param val The real state
	 * @return the corresponding grid ID.
	 */
	size_t val_to_id(std::vector<double> val);

	/**
	 * \brief Get the center value x of a grid with ID number i.
	 * 
	 * @param[in,out] val The real state
	 * @param[in] the corresponding grid ID.
	 */
	void id_to_val(std::vector<double> &val, size_t id) const;

	/**
	 * \brief Get the grids covered by an interval area.
	 *
	 * \f$i = i_1 + N_1 i_2 + \cdots + N_1 N_2 \cdots N_{n-1} i_n \f$
	 *
	 * @param[in] box The interval
	 * @param[in] boxin 0: allow box out of range, 1: don't allow
	 * @param[in] strictin 1: collect the grids that entirely inside the box,
	 *             0: collect the grids that intersect the box
	 * @return a list of grid IDs.
	 */
	std::vector<size_t> subset(ivec &box, bool boxin, bool strictin);

	/**
	 * \brief Get the neighbour grids of the given grid point.
	 *
	 * \f$i = i_1 + N_1 i_2 + \cdots + N_1 N_2 \cdots N_{n-1} i_n \f$
	 *
	 * @param[in] id The id of the given grid point
	 * @return a list of neighbour grids ID.
	 */
	std::vector<size_t> neighbours(size_t id);

	/**
	 * \brief Print grid information.
	 */
	void print_info();
  
	/**
	 * \brief Print grid centers. 
	 */
	void print_data();

    };

} // namespace rocs

#endif

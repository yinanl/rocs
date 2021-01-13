/**
 *  safeset.cpp
 *
 *  Compute the safety (collision-free) set of two agents.
 *
 *  Created by Yinan Li on July 13, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"
#include "src/hdf5io.h"


struct twoagent {
	static const int n = 3;  // system dimension
	static const int nu = 2;  // control dimension
	rocs::ivec d{rocs::interval(-0.8, 0.8),
		     rocs::interval(-0.8, 0.8)};

	/* template constructor
	 * @param[out] dx
	 * @param[in] x = [xr, yr, psir]
	 * @param u = [v, w]
	 * @param d = [v', w']
	 */
	template<typename S>
	twoagent(S *dx, const S *x, rocs::Rn u) {
		dx[0] = -u[0] + d[0]*cos(x[2]) + u[1]*x[1];
		dx[1] = d[0]*sin(x[2]) - u[1]*x[0];
		dx[2] = d[1] - u[1];
	}
};


/** safety constraint f(x)<=0
 * x^2+y^2>=d^2 <==> -(x^2+y^2)+d^2<=0
 *
 **/
template<typename T>
T collision_free(const T &x) {
  const double d = 1.2;
  T y(1);
  y[0] = d*d - x[0]*x[0] - x[1]*x[1];
  return y;
}
int main()
{
    /**
     * Set the state and control space.
     */
    double xlb[] = {-3.0, -3.0, -M_PI};
    double xub[] = {3.0, 3.0, M_PI};
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /**
     * Define the two-agent system
     */
    double t = 0.3;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    rocs::CTCntlSys<twoagent> safety("collision-free", t,
				     twoagent::n, twoagent::nu,
				     delta, &controlparams);

    safety.init_workspace(xlb, xub);
    safety.init_inputset(mu, ulb, uub);
    safety.allocate_flows();

    rocs::CSolver solver(&safety, 0, rocs::RELMAX);
    double e[]{0.2, 0.2, 15.0};
    solver.init(rocs::GOAL, &collision_free<rocs::ivec>, e);
    solver.init_goal_area();
    double eta[]{0.2, 0.2, 0.2};
    solver.invariance_control(&safety, eta);
    solver.print_controller_info();

    // rocs::matWriter wtr("data_safe_0.8.mat");
    // wtr.open();
    // wtr.write_problem_setting(safety, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string filename = "controller_safety_new_0.8-1.2-0.3-0.2.h5";
    rocs::h5FileHandler wtr(filename, H5F_ACC_TRUNC);
    wtr.write_problem_setting< rocs::CTCntlSys<twoagent> >(safety);
    wtr.write_ivec_array(solver._goal, "G");
    wtr.write_ivec_array(solver._obs, "xobs");
    wtr.write_sptree_controller(solver);

    safety.release_flows();
    return 0;
}

/*
 *  Simulate the real boundary of the region of attraction (ROA) 
 *  for the reversed VDP example.
 *
 *  Created by Yinan Li on May 19, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */



#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include "src/definitions.h"
#include "src/matlabio.h"


/* user defined dynamics */
void vdp(const rocs::Rn &x, rocs::Rn &dxdt, double t) {
    dxdt[0] = -x[1];
    dxdt[1] = x[0] + x[1]*(x[0]*x[0]-1);
}



int main()
{
    double xlb[] = {-4, -4};
    double xub[] = {4, 4};
    double eta[] = {0.001, 0.001};
    double T = 10.0;
    double dt = 0.01;

    rocs::grid xg(2, eta, xlb, xub);
    xg.gridding();
    
    /* ROA computation by running simulations */
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;
    std::vector<bool> roa(xg._nv);
    rocs::Rn x;
    bool stable;
    for (int i = 0; i < xg._nv; ++i) {
	
	x = xg._data[i];
	
	// boost::numeric::odeint::integrate(vdp, x, 0.0, T);
	// boost::numeric::odeint::integrate_adaptive(make_controlled(1E-12, 1E-12, dopri5()), vdp, x, 0.0, T, dt);
	boost::numeric::odeint::integrate_const(rk45, vdp, x, 0.0, T, dt);

	stable = true;
	for (int j = 0; j < xg._dim; ++j) {
	    // std::cout << x[j] << ", "; 
	    stable = stable & (x[j] < 0.5) & (x[j] > -0.5);
	}
	// std::cout <<'\n';

	if (stable) 
	    roa[i] = true;
    }
    
    

    /* write roa boundary to a mat file */
    rocs::matWriter wtr("data_roavdp_real.mat");
    wtr.open();

    wtr.write_boundary(xg, roa, "bd0");
    wtr.close();
    
    return 0;
}

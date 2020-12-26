/**
 *  scara.h
 *
 *  The dynamics of a two-link SCARA manipulator
 *
 *  Created by Yinan Li on July 07, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _scara_h
#define _scara_h

const double m1 = 0.1;
const double m2 = 0.1;
const double l1 = 0.15;
const double l2 = 0.15;
const double r1 = 0.5*l1;
const double r2 = 0.5*l2;
const double I1 = 1.33e-5;
const double I2 = 1.33e-5;


/**
 * The full dynamics
 */
struct scaraode {

    static const int n = 4;  // system dimension
    static const int nu = 2;  // control dimension

    double z1; // = I1 + I2 +m1*r1*r1 + m2*(l1*l1+r2*r2);
    double z2; // = m2*l1*r2;
    double z3; // = I2 + m2*r2*r2;

    /* template constructor
     * @param[out] dx
     * @param[in] x = [theta1, theta2, w1, w2]
     * @param u = [tau1, tau2]
     */
    // scaraode(T<rocs::interval> *dx, const T<rocs::interval> *x, rocs::Rn u) :
	// z1(I1+I2+m1*r1*r1+m2*(l1*l1+r2*r2)), z2(m2*l1*r2), z3(I2+m2*r2*r2) {
	// auto detM = z3*(z1-z3) - z2*z2*cos(x[1])*cos(x[1]);
	// auto a = z2*sin(x[1])*(2*x[2]+x[3])*x[3];
	// auto b = z2*cos(x[1]);
	// auto c = z2*x[2]*x[2]*sin(x[1])-u[1];
	// dx[0] = x[2];
	// dx[1] = x[3];
	// dx[2] = (z3*u[0] + z3*a + (z3+b)*c) / detM;
	// dx[3] = ((z1+2*b)*(-c) - (z3+b)*(u[0]+a)) / detM;
	// // dx[2] = (z3*u[0] + z3*(z2*sin(x[1])*(2*x[2]+x[3])*x[3]) +
	// // 	 (z3+z2*cos(x[1]))*(z2*x[2]*sin(x[1])-u[1])) /
	// //     (z3*(z1-z3) - z2*z2*cos(x[1])*cos(x[1]));
	// // dx[3] = ((z1+2*(z2*cos(x[1])))*(u[1]-z2*x[2]*sin(x[1])) -
	// // 	 (z3+z2*cos(x[1]))*(u[0]+z2*sin(x[1])*(2*x[2]+x[3])*x[3])) /
	// //     (z3*(z1-z3) - z2*z2*cos(x[1])*cos(x[1]));
    // }
    // scaraode(double *dx, const double *x, rocs::Rn u) :
    // 	z1(I1+I2+m1*r1*r1+m2*(l1*l1+r2*r2)), z2(m2*l1*r2), z3(I2+m2*r2*r2) {
    // 	double detM = z3*(z1-z3) - z2*z2*std::cos(x[1])*std::cos(x[1]);
    // 	double a = z2*std::sin(x[1])*(2*x[2]+x[3])*x[3];
    // 	double b = z2*std::cos(x[1]);
    // 	double c = z2*x[2]*x[2]*std::sin(x[1])-u[1];

    // 	dx[0] = x[2];
    // 	dx[1] = x[3];
    // 	dx[2] = (z3*u[0] + z3*a + (z3+b)*c) / detM;
    // 	dx[3] = ((z1+2*b)*(-c) - (z3+b)*(u[0]+a)) / detM;
    // }

    template<typename S>
    scaraode(S *dx, const S *x, rocs::Rn u) :
    	z1(I1+I2+m1*r1*r1+m2*(l1*l1+r2*r2)), z2(m2*l1*r2), z3(I2+m2*r2*r2) {
    	dx[0] = x[2];
    	dx[1] = x[3];
    	dx[2] = (z3*u[0] + z3*(z2*sin(x[1])*(2*x[2]+x[3])*x[3]) +
    		 (z3+z2*cos(x[1]))*(z2*x[2]*x[2]*sin(x[1])-u[1])) /
    	    (z3*(z1-z3) - z2*z2*cos(x[1])*cos(x[1]));
    	dx[3] = ((z1+2*(z2*cos(x[1])))*(u[1]-z2*x[2]*x[2]*sin(x[1])) -
    		 (z3+z2*cos(x[1]))*(u[0]+z2*sin(x[1])*(2*x[2]+x[3])*x[3])) /
    	    (z3*(z1-z3) - z2*z2*cos(x[1])*cos(x[1]));
    }
};


/**
 * The discrete-time 2 double integrators:
 * theta1_dot = w1, theta2_dot = w2, w1_dot = u1, w2_dot = u2
 */
const double ts = 0.1; // the sampling time
struct integrator {
    static const int n = 4;  // system dimension
    static const int m = 2;  // control dimension

    /* template constructor
     * @param[out] dx
     * @param[in] x = [theta1, theta2, w1, w2]
     * @param u = [acc]
     */
    template<typename S>
    integrator(S &dx, const S &x, rocs::Rn u) {
		dx[0] = x[0] + ts*x[2] + 0.5*ts*ts*u[0];
		dx[2] = x[2] + ts*u[0];

		dx[1] = x[1] + ts*x[3] + 0.5*ts*ts*u[1];
		dx[3] = x[3] + ts*u[1];
    }
};


struct simple {
    static const int n = 2;  // system dimension
    static const int m = 1;  // control dimension
    template<typename S>
    simple(S &dx, const S &x, rocs::Rn u) {
		dx[0] = x[0] + ts*x[1] + 0.5*ts*ts*u[0];
		dx[1] = x[1] + ts*u[0];
    }
};


/*
 * Convert (x, y) in the operational space to (theta1, theta2) in the joint space
 */
void xy2theta(rocs::ivec &theta, rocs::ivec &x) {
	rocs::interval l3 = sqrt(x[0]*x[0] + x[1]*x[1]);
	rocs::interval a = acos((l3*l3 +l1*l1 - l2*l2)/(2*l1*l3));
	rocs::interval b = acos((l3*l3 +l2*l2 - l1*l1)/(2*l2*l3));
	theta[0] = atan(x[1]/x[0]) + a;
	theta[1] = -(a+b);
	theta[2] = atan(x[1]/x[0]) - a;
	theta[3] = a+b;
}

/*
 * Convert the joint space to the operational space
 */
void theta2xy(rocs::ivec &x, rocs::ivec &theta) {
	x[0] = l1 * cos(theta[0]) + l2 * cos(theta[0] + theta[1]);
	x[1] = l1 * sin(theta[0]) + l2 * sin(theta[0] + theta[1]);
}

#endif

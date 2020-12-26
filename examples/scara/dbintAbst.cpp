/**
 *  dbintAbst.cpp
 *
 *  Abstraction-based DBA control of the simplified SCARA manipulator dynamics (the double integrator model).
 *
 *  Created by Yinan Li on August 8, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>
#include <sys/stat.h>

#include "src/DBAparser.h"
#include "src/system.hpp"
#include "src/abstraction.hpp"
#include "src/bsolver.hpp"

#include "src/hdf5io.h"

#include "scara.hpp"


const double h = 0.8*l1;
const double r = 0.5*l1;
double a1 = atan(h / r); // the upper bound for theta1
double a2 = asin(h / l1);



int main(int argc, char *argv[])
{
    /* Input arguments: 
     * carAbst dbafile 
     */
    if (argc != 2) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }


    clock_t tb, te;
    /**
     * Setup the control system.
     */
    double xlb[] = {0, -M_PI, -1, -1};
    double xub[] = {M_PI/2.0, M_PI, 1, 1};
    double ulb[] = {-5.0, -5.0};
    double uub[] = {5.0, 5.0};
    double mu[] = {0.5, 0.5};

    rocs::DTCntlSys<integrator> scara("dba", ts, integrator::n, integrator::m);
    scara.init_workspace(xlb, xub);
    scara.init_inputset(mu, ulb, uub);


    /**
     * Abstraction 
     */
    rocs::abstraction< rocs::DTCntlSys<integrator> > abst(&scara);
    const double eta[]{0.05, 0.05, 0.1, 0.1};
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';

    /* Assign labels to states: has to be consistent with the dba file */
    double glb[][4] = {{0.4980, 1.5739, -0.1, -0.1},
		       {0.5103, -0.9363, -0.1, -0.1}};
    double gub[][4] = {{0.5772, 1.7055, 0.1, 0.1},
		       {0.5769, -0.8363, 0.1, 0.1}};
    rocs::UintSmall nGoal = 2;
    auto label_target = [&glb, &gub, &nGoal, &abst, &eta](size_t i) {
			    std::vector<double> x(abst._x._dim);
			    abst._x.id_to_val(x, i);
			    double c[]{eta[0]/2.0,eta[1]/2.0,eta[2]/2.0,eta[3]/2.0};
			    boost::dynamic_bitset<> label(nGoal, false); // n is the number of goals
			    for(rocs::UintSmall i = 0; i < nGoal; ++i) {
				label[i] = (glb[i][0] <= (x[0]-c[0]) && (x[0]+c[0]) <= gub[i][0] && 
					    glb[i][1] <= (x[1]-c[1]) && (x[1]+c[1]) <= gub[i][1])
				    ? true: false;
			    }
			    return label.to_ulong();
			};
    abst.assign_labels(label_target);
    std::cout << "Specification assignment is done.\n";
    std::vector<size_t> targetIDs;
    std::vector<rocs::Rn> targetPts;   //initial invariant set
    rocs::Rn x(abst._x._dim);
    for(size_t i = 0; i < abst._labels.size(); ++i) {
    	if(abst._labels[i] > 0) {
    	    targetIDs.push_back(i);
    	    abst._x.id_to_val(x, i);
    	    targetPts.push_back(x);
    	}
    }
    /* Mark obstacles */
    double obs[2][4] = {{a1, xlb[1], xlb[2], xlb[3]},
    			{M_PI/2.0, xub[1], xub[2], xub[3]}};
    const double e[]{0.01, 0.01, 3, 3};  // only bisect x[0] and x[1].
    auto avoid = [&obs, abst, eta](size_t& id) {
		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, id);
    		     double c[]{eta[0]/2.0,eta[1]/2.0,eta[2]/2.0,eta[3]/2.0}; //+1e-10;
		     if((obs[0][0]-c[0]) <= x[0] && x[0] <= (obs[1][0]+c[0]) &&
			(obs[0][1]-c[1]) <= x[1] && x[1] <= (obs[1][1]+c[1]) &&
			(obs[0][2]-c[2]) <= x[2] && x[2] <= (obs[1][2]+c[2]) &&
			(obs[0][3]-c[3]) <= x[3] && x[3] <= (obs[1][3]+c[3]))
			 return -1;
		     if(x[0]<=a2 &&
			x[0]+x[1]+atan2(h-l1*sin(x[0]),l1*cos(x[0])-r)>=M_PI )
			 return -1;
		     if(x[0]<=a1 && x[0]>=a2 &&
			x[0]+x[1]-atan2(l1*sin(x[0])-h,l1*cos(x[0]))>=M_PI )
			 return -1;
		     return 0;
		 };
    abst.assign_labels(avoid);
    abst.assign_label_outofdomain(0);
    std::vector<size_t> obstacles;
    for (size_t i = 0; i < abst._x._nv; ++i) {
    	if (abst._labels[i] < 0)
    	    obstacles.push_back(i);
    }

    /* Compute abstraction */
    std::string transfile = "abst_0.05-0.01.h5";
    struct stat buffer;
    float tabst;
    if(stat(transfile.c_str(), &buffer) == 0) {
	/* Read from a file */
	std::cout << "Transitions have been computed. Reading transitions...\n";
	rocs::h5FileHandler transRdr(transfile, H5F_ACC_RDONLY);
	tb = clock();
	transRdr.read_transitions(abst._ts);
	te = clock();
	tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of reading abstraction: " << tabst << '\n';
    } else {
	/* Robustness margins */
	double e1[] = {0.0, 0.0, 0.0, 0.0};
	double e2[] = {0.0, 0.0, 0.0, 0.0};
	tb = clock();
	abst.assign_transitions(e1, e2);
	te = clock();
	float tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of computing abstraction: " << tabst << '\n';
	/* Write abstraction to file */
	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
	transWtr.write_transitions(abst._ts);
    }
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
    

    /**
     * Read specification file
     **/
    std::cout << "Loading the specification...\n";
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    std::string specfile = std::string(argv[1]);
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc))
	std::exit(1);


    /**
     * Solve a Buchi game on the product of NTS and DBA.
     */
    std::cout << "Start solving a Buchi game on the product of the abstraction and DBA...\n";
    rocs::BSolver solver; // memories will be allocated for psolver
    solver.construct_dba((int)nAP, (int)nNodes, (int)q0, acc, arrayM);
    tb = clock();
    solver.load_abstraction(abst);
    solver.generate_product(abst);
    solver.solve_buchigame_on_product();
    te = clock();
    float tsyn = (float)(te - tb)/CLOCKS_PER_SEC;
    std::cout << "Time of synthesizing controller: " << tsyn << '\n';

    
    /**
     * Display and save memoryless controllers.
     */
    std::cout << "Writing the controller...\n";
    std::string datafile = "controller_abst_0.05-0.01.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<integrator> >(scara);
    ctlrWtr.write_array<size_t>(targetIDs, "G");
    ctlrWtr.write_array<size_t>(obstacles, "A");
    ctlrWtr.write_array<double>(eta, integrator::n, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(&(solver._sol));
    std::cout << "Controller writing is done.\n";

    std::cout << "Total time of used (abstraction+synthesis): " << tabst+tsyn << '\n';

    return 0;
}

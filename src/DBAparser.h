/**
 * DBAparser.h
 * The class that parse a input file to get information of the specification.
 * Created by Yinan Li on Feb. 01, 2020.
 * Hybrid Systems Group, University of Waterloo
 */

#ifndef _dbaparser_h
#define _dbaparser_h

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "definitions.h"

namespace rocs {

    class DBAparser {
    private:
	std::ifstream _file;
	std::string _line;

    public:
	bool open(const std::string filename);
	void close();

	void back_to_the_first_line();

	UintSmall read_number_of_nodes();
	UintSmall read_number_of_atomic_propsitions();
	UintSmall read_number_of_propositions();
	UintSmall read_initial_state();
	void read_accepting_nodes(std::vector<UintSmall> &acc);
	
	void read_spec_name();
	void read_propositions();
	void read_transition_matrix(std::vector< std::vector<UintSmall> > &arrayM);
	
    };

    bool read_spec(std::string specfile, UintSmall &nNodes, UintSmall &nProps,
			      std::vector< std::vector<UintSmall> > &arrayM,
			      std::vector<rocs::UintSmall> &acc);
    bool read_spec(std::string specfile,
		    UintSmall &nNodes, UintSmall &nAP, UintSmall &q0,
		    std::vector< std::vector<UintSmall> > &arrayM,
		    std::vector<rocs::UintSmall> &acc);
}


#endif

/**
 * DBAparser.h
 * The source file for class DBAparser.
 * Created by Yinan Li on Feb. 01, 2020.
 * Hybrid Systems Group, University of Waterloo
*/


#include "DBAparser.h"

namespace rocs {

    bool DBAparser::open(const std::string filename) {
	    _file.open(filename);
	    return _file.good();
	}
    void DBAparser::close() {
	_file.close();
    }

    void DBAparser::back_to_the_first_line() {
	_file.clear();
	_file.seekg(0,std::ios::beg);
    }

    UintSmall DBAparser::read_number_of_nodes() {
	back_to_the_first_line();
	while (std::getline(_file, _line)) {
	    if (_line.find("nNodes") != std::string::npos) {
		std::vector<std::string> tokens;
		boost::split(tokens, _line, boost::is_any_of("="));
		return std::stoi(tokens[1]);
	    }
	}
	return 0;
    }
    UintSmall DBAparser::read_number_of_atomic_propsitions() {
	back_to_the_first_line();
	while (std::getline(_file, _line)) {
	    if (_line.find("nAP") != std::string::npos) {
		std::vector<std::string> tokens;
		boost::split(tokens, _line, boost::is_any_of("="));
		return std::stoi(tokens[1]);
	    }
	}
	return 0;
    }
    UintSmall DBAparser::read_number_of_propositions() {
	back_to_the_first_line();
	while (std::getline(_file, _line)) {
	    if (_line.find("nProps") != std::string::npos) {
		std::vector<std::string> tokens;
		boost::split(tokens, _line, boost::is_any_of("="));
		return std::stoi(tokens[1]);
	    }
	}
	return 0;
    }
    UintSmall DBAparser::read_initial_state() {
	back_to_the_first_line();
	while (std::getline(_file, _line)) {
	    if (_line.find("init") != std::string::npos) {
		std::vector<std::string> tokens;
		boost::split(tokens, _line, boost::is_any_of("="));
		return std::stoi(tokens[1]);
	    }
	}
	return 0;
    }
    void DBAparser::read_spec_name() {
	back_to_the_first_line();
	while (std::getline(_file, _line)) {
	    if (_line.find("Name") != std::string::npos) {
		std::cout << "Formula: " << _line.substr(_line.find("=") + 1) << '\n';
		break;
	    }
	}
    }
    void DBAparser::read_propositions() {
	back_to_the_first_line();
	while (std::getline(_file, _line)) {
	    if (_line.find("AP") != std::string::npos) {
		std::cout << "Propositions: " << _line.substr(_line.find("=") + 1) << '\n';
		break;
	    }
	}
    }
    void DBAparser::read_transition_matrix(std::vector< std::vector<UintSmall> > &arrayM) {
	assert(arrayM.size()>0);
	back_to_the_first_line();
	bool M = false;
	while (std::getline(_file, _line)) {
	    if (M) {
		if (_line.find("]")!=std::string::npos) {
		    break;
		} else {
		    std::vector<std::string> tokens;
		    boost::split(tokens, _line, boost::is_any_of(","));
		    arrayM[std::stoi(tokens[0])][std::stoi(tokens[1])] = std::stoi(tokens[2]);
		}
	    } else {
		if (_line.find("M=[") != std::string::npos) {
		    M = true;
		    continue;
		}
	    }
	}
	if (!M)
	    std::cout << "No transition matrix found.\n";
    }
    void DBAparser::read_accepting_nodes(std::vector<UintSmall> &acc) {
	back_to_the_first_line();
	acc.clear();
	while (std::getline(_file, _line)) {
	    if (_line.find("acc") != std::string::npos) {
		std::vector<std::string> tokens;
		std::string s = _line.substr(_line.find("=") + 1);
		boost::split(tokens, s, boost::is_any_of(","));
		for (UintSmall i = 0; i < tokens.size(); ++i)
		    acc.push_back(std::stoi(tokens[i]));
		break;
	    }
	}
    }

    /* friend function */
    bool read_spec(std::string specfile,
		   UintSmall &nNodes, UintSmall &nProps,
		   std::vector< std::vector<UintSmall> > &arrayM,
		   std::vector<rocs::UintSmall> &acc) {

	DBAparser dba;
	if(!dba.open(specfile)) {/* open the file of specification */
	    std::cout << "Error in opening " << specfile << ".\n";
	    return 0;
	}
    
	dba.read_spec_name(); /* Display the specification name. */
	dba.read_propositions(); /* Display the propositions. */
	nProps = dba.read_number_of_propositions(); /* The number of DBA nodes */
	nNodes = dba.read_number_of_nodes(); /* The number of propositions */
	/* Setup the 2d array by the transition matrix M */
	arrayM.resize(nNodes, std::vector<UintSmall>(nProps, nNodes));
	dba.read_transition_matrix(arrayM);
	/* Get the accepting nodes */
	dba.read_accepting_nodes(acc);
	
	dba.close();

	std::cout << "Number of propositions: " << nProps << '\n';
	std::cout << "Number of nodes: " << nNodes << '\n';
	std::cout << "The transition matrix is:\n";
	for (rocs::UintSmall i = 0; i < nNodes; ++i) {
	    for (rocs::UintSmall j = 0; j < nProps; ++j)
		std::cout << arrayM[i][j] << ' ';
	    std::cout << '\n';
	}
	std::cout << "The accepting nodes are: ";
	for (rocs::UintSmall i = 0; i < acc.size(); ++i) {
	    std::cout << acc[i] << ' ';
	}
	std::cout << '\n';

	return 1;
    }

    bool read_spec(std::string specfile,
		    UintSmall &nNodes, UintSmall &nAP, UintSmall &q0,
		    std::vector< std::vector<UintSmall> > &arrayM,
		    std::vector<rocs::UintSmall> &acc) {

	DBAparser dba;
	if(!dba.open(specfile)) {/* open the file of specification */
	    std::cout << "Error in opening " << specfile << ".\n";
	    return 0;
	}
    
	dba.read_spec_name(); /* Display the specification name. */
	dba.read_propositions(); /* Display the propositions. */
	nAP = dba.read_number_of_atomic_propsitions(); /* The number of atomic propositions */
	nNodes = dba.read_number_of_nodes(); /* The number of DBA nodes */
	q0 = dba.read_initial_state();
	/* Setup the 2d array by the transition matrix M */
	UintSmall nProps = 1; // nProps = 2^nAP
	for(UintSmall i = 0; i < nAP; ++i)
	    nProps *= 2;
	arrayM.resize(nNodes, std::vector<UintSmall>(nProps, nNodes));
	dba.read_transition_matrix(arrayM);
	/* Get the accepting nodes */
	dba.read_accepting_nodes(acc);
	
	dba.close();

	std::cout << "Number of propositions: " << nProps << '\n';
	std::cout << "Number of nodes: " << nNodes << '\n';
	std::cout << "The transition matrix is:\n";
	for (rocs::UintSmall i = 0; i < nNodes; ++i) {
	    for (rocs::UintSmall j = 0; j < nProps; ++j)
		std::cout << arrayM[i][j] << ' ';
	    std::cout << '\n';
	}
	std::cout << "The accepting nodes are: ";
	for (rocs::UintSmall i = 0; i < acc.size(); ++i) {
	    std::cout << acc[i] << ' ';
	}
	std::cout << '\n';

	return 1;
    }
    
} //namespace rocs

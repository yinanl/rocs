/**
 * The header file of the class DBAparser.
 *
 * Created by Yinan Li on Feb. 01, 2020.
 *
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

    /**
     * \brief A class that parses a DBA given in a file.
     *
     * The info extracted by parsing:
     * - # of %DBA states (nNodes)
     * - # of atomic propositions (nAP)
     * - # of propositions
     * - the initial state (q0)
     * - the accepting states (acc)
     */
    class DBAparser {
    private:
	std::ifstream _file; /**< An input file stream object */
	std::string _line; /**< A string object */

    public:
	/**
	 * \brief Open a file.
	 *
	 * @param[in] filename The name of the file
	 * @return whether or not the file is opened.
	 */
	bool open(const std::string filename);
	
	/**
	 * \brief Close the file that has been opened.
	 */
	void close();
	
	/**
	 * \brief Go back to the first line of the file.
	 */
	void back_to_the_first_line();
	
	/**
	 * \brief Read the number of states in the given %DBA file.
	 *
	 * @return the nNodes (an UintSmall variable).
	 */
	UintSmall read_number_of_nodes();
	
	/**
	 * \brief Read the number of atomic propositions in the given %DBA file.
	 *
	 * @return the nAP (an UintSmall variable).
	 */
	UintSmall read_number_of_atomic_propsitions();
	
	// /**
	//  * \brief Read the number of propositions in the given %DBA file.
	//  * @return the # of propositions (an UintSmall variable).
	//  */
	// UintSmall read_number_of_propositions();
	
	/**
	 * \brief Read the initial state in the given %DBA file.
	 *
	 * @return the q0 (an UintSmall variable).
	 */
	UintSmall read_initial_state();
	
	/**
	 * \brief Read the accepting states in the given %DBA file.
	 *
	 * @param[in,out] acc the lvalue reference to a vector of small integers.
	 */
	void read_accepting_nodes(std::vector<UintSmall> &acc);
	
	/**
	 * \brief Read and display the LTL specification.
	 */
	void read_spec_name();
	
	/**
	 * \brief Read and display the propositions.
	 */
	void read_propositions();
	
	/**
	 * \brief Read the transition matrix.
	 *
	 * @param[in,out] arrayM the lvalue reference to the transition matrix (a 2d array of UintSmall variables)
	 */
	void read_transition_matrix(std::vector< std::vector<UintSmall> > &arrayM);
	
    };

    
    // bool read_spec(std::string specfile, UintSmall &nNodes, UintSmall &nProps,
    // 		   std::vector< std::vector<UintSmall> > &arrayM,
    // 		   std::vector<rocs::UintSmall> &acc);
    
    /**
     * \brief The non-member function of DBAparser that returns the DBA info.
     *
     * @param[in] specfile The name of the specification file
     * @param[in,out] nNodes The # of states (an UintSmall variable)
     * @param[in,out] nAP The # of atomic propositions (an UintSmall variable)
     * @param[in,out] q0 The initial state (UintSmall type)
     * @param[in,out] arrayM The transition matrix (a 2d array of UintSmall variables)
     * @param[in,out] acc The accepting states (a vector of UintSmall variables)
     */
    bool read_spec(std::string specfile,
		   UintSmall &nNodes, UintSmall &nAP, UintSmall &q0,
		   std::vector< std::vector<UintSmall> > &arrayM,
		   std::vector<rocs::UintSmall> &acc);
}


#endif

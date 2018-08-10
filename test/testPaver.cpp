//
//  testPaver.cpp
//  
//  test interval subpaving class
//
//  Created by yinan li on 31/08/2016.
//  Copyright (c) 2016 UW. All rights reserved.
// ------------------------------------------------


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PaverClass


#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>
#include <cmath>
#include "src/interval_paver.h"


/*
  test case 1: initialization of SPnode
*/
BOOST_AUTO_TEST_CASE(test_init_SPnode)
{
    rocs::ivec a(3);
    a[0] = rocs::interval(-1, 1);
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);
    
    rocs::ivec e(2);
    e[0] = rocs::interval(1, 0);
    e[1] = rocs::interval(rocs::NINF, rocs::PINF);

    /* test direct init */
    rocs::SPnode node0(a, 5); //init
    boost::dynamic_bitset<> cntl0(5, 0);
    BOOST_CHECK(node0._box == a);  // box=a
    BOOST_CHECK(node0._cntl == cntl0);  // cntl={}
    BOOST_CHECK(node0._tag == 0);  // tag=0
    BOOST_CHECK(node0._b0);  // b0=1
    BOOST_CHECK(!node0._b1);  // b1=0
    BOOST_CHECK(node0._split == -1);  // no splitting dimension
    node0._tag = 1;
    BOOST_CHECK(node0._tag == 1);

    /* test default init */
    rocs::SPnode node1;
    BOOST_CHECK(node1.isempty());

    /* test 2nd init */
    boost::dynamic_bitset<> cntl2(5, 0);
    cntl2[3] = 1;
    rocs::SPnode node2 = rocs::SPnode(a, cntl2, 2, true, false);
    BOOST_CHECK(node2._box == a);  // box=a
    BOOST_CHECK(2 == node2._tag);  // tag=2
    BOOST_CHECK(node2._cntl == cntl2);  // cntl={0,0,1}
    BOOST_CHECK(node2._b0);  // b0=1
    BOOST_CHECK(!node2._b1);  // b1=0
    BOOST_CHECK(node2._split == -1);  // no splitting

    /* test copy init */
    rocs::SPnode node3 = rocs::SPnode(node2);
    BOOST_CHECK(node3._box == a);  // box=a
    BOOST_CHECK(node3._tag == 2);  // tag=2
    BOOST_CHECK(node3._cntl == cntl2);  // cntl={0,0,1}
    BOOST_CHECK(node3._split == -1);  // no splitting

    
    /* test printing */
    std::cout << node0 << std::endl;
    std::cout << node1 << std::endl;
    std::cout << node2 << std::endl;  // node3=node2
}


/*
  test case 2: initialization of SPtree
*/
BOOST_AUTO_TEST_CASE(test_SPtree)
{
    rocs::ivec a(3);
    a[0] = rocs::interval(-1, 1);
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 3.2);

    rocs::SPnode root(a, 5);  // box=a, tag=0, cntl={0,0,0,0,0}


    rocs::SPtree paver1; // default init (empty)
    rocs::SPtree paver2(&root); // init by root pointer
    rocs::SPtree paver3 = rocs::SPtree(paver2); // copy init


    /* expand 2 levels */
    paver3.expand(paver3._root, 2);  // split=2, split along axis 2
    paver3.expand(paver3._root->_left, 2);  // split=2
    paver3.expand(paver3._root->_right, 0);  // split=0

    /* test depth 1 (root--0,leaf--2) */
    rocs::ivec al(3), ar(3);
    al[0] = rocs::interval(-1, 1);
    al[1] = rocs::interval(-0.1, 0.01);
    al[2] = rocs::interval(0, 1.6);
    ar[0] = rocs::interval(-1, 1);
    ar[1] = rocs::interval(-0.1, 0.01);
    ar[2] = rocs::interval(1.6, 3.2);
    rocs::ivec lbox = paver3._root->_left->_box;
    rocs::ivec rbox = paver3._root->_right->_box;
    BOOST_CHECK(lbox == al);
    BOOST_CHECK(rbox == ar);
    BOOST_CHECK(paver3._root->_left->_split == 2);
    BOOST_CHECK(paver3._root->_right->_split == 0);
    BOOST_CHECK(paver3.height(paver3._root) == 3);
    BOOST_CHECK(paver3.leafcount(paver3._root) == 4);
    BOOST_CHECK(paver3.nodecount(paver3._root) == 7);


    /* test printing */
    std::cout << "print paver3:" << std::endl;
    paver3.print_subtree(paver3._root); //printing
    std::cout << std::endl;
    std::cout << "print leaves of paver3:" << std::endl;
    paver3.print_leaves(paver3._root);
    std::cout << std::endl;
    

    /* test expand */
    boost::dynamic_bitset<> cntl2(7, false);
    cntl2[1] = true;
    cntl2[6] = true;
    paver2._root->_cntl = cntl2;
    paver2._root->_tag = 2;

    paver2.expand(paver2._root, 1);

    std::cout << "print paver2:" << std::endl;
    paver2.print_subtree(paver2._root); //printing
    std::cout << std::endl;
    std::cout << "print leaves of paver2:" << std::endl;
    paver2.print_leaves(paver2._root);
    std::cout << std::endl;
    paver1 = paver2;  // copy assignment
    std::cout << "print leaves of paver1:" << std::endl;
    paver1.print_leaves(paver1._root);
    std::cout << std::endl;
}

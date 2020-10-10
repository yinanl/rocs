/**
 *  carAbst_ref.cpp
 *  Abstraction-based DBA control synthesis for vehicle motion planning.
 *
 *  This is a reference file for testing the integrated study case in carAbst.cpp.
 *
 *  Created by Zhibing Sun on May 26, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <math.h>

#include "src/abstraction.hpp"
#include "src/txtfileio.h"
#include "TicToc.hh"
#include "carAbst_ref.h"
#include "car.hpp"


void input_labels_and_pre(HEAD *head)
{
    NTS_PRE *pre;
    NODE_PRE *nts_pre;
    int action,a,*outdeg;
    long long x1,x2,x,n0,m0;
    EDGE *in_edge,*p5;
    DBA *dba=&head->dba;
    BOOL *labels;
    int k;
    
    pre=&head->nts_pre;
    
    /* to measure time */
    TicToc tt;
    
     clock_t tb, te;
    
    /* set the state space */
    const double theta = 3.5;
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};
    const double eta[] = {.2,.2,.2}; /* set precision */
    
    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* define the control system */
    rocs::DTCntlSys<carde> car("reach goal", h, carde::n, carde::m);
    car.init_workspace(xlb, xub);
    car.init_inputset(mu, ulb, uub);

    /* Abstraction */
    // std::string scale_eta ("1");
    // if (argc>1) {
    // 	    scale_eta = argv[1];
    // }
    // double set_eta = 0.2/pow(2,std::stod(scale_eta)-1);
    
    rocs::abstraction<rocs::DTCntlSys<carde>> abst(&car);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';

    double obs[3][4] = {{0.0, 3.2, 4.0, 5.0},
			{5.4, 6.0, 5.0, 10.0},
			{4.5, 5.2, 0.0, 2.5}};
    
    /* avoid function returns 1 if x is in avoid set  */
    rocs::UintSmall nAvoid = 3;
    auto label_avoid = [&obs, &nAvoid, &abst, &eta](size_t i) {
    		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, i);
    		     double c1= eta[0]/2.0+1e-10;
    		     double c2= eta[1]/2.0+1e-10;
    		     for(size_t i = 0; i < nAvoid; ++i) {
    			 if ((obs[i][0]-c1) <= x[0] && x[0] <= (obs[i][1]+c1) && 
    			     (obs[i][2]-c2) <= x[1] && x[1] <= (obs[i][3]+c2))
    			     return -1;
    		     }
    		     return 0;
    		 };
    abst.assign_labels(label_avoid);
    
    /* define target set */
    dba->k=k=6; //  number of labels
    
    double targets[6][4] = {
        
        {1.0, 2.0, 0.5, 2.0},
        {7.1, 9.1, 1.9, 2.9},
        {0.0, 4.0, 6.0, 10.0},
        {6.0, 10.0, 0.0, 4.0},
        {0.5, 2.5, 7.5, 8.5},
        {3.8, 4.6, 3.1, 4.5}
        
    };
    
    
    /* define target set */
    auto target = [&targets,&abst,&eta](int i,const size_t idx) {
        std::vector<double> x(abst._x._dim);
	abst._x.id_to_val(x, idx);
	double c1= eta[0]/2.0; //+1e-10;
	double c2= eta[1]/2.0; //+1e-10;
        /* function returns 1 if cell associated with x is in target set  */
        return (targets[i][0] <= (x[0]-c1) && (x[0]+c1) <= targets[i][1] &&
                targets[i][2] <= (x[1]-c2) && (x[1]+c2) <= targets[i][3])?true:false;
    };
    
    // MARK: input labels
    printf("Input labels:\n");
    /*
     "labelsX_X.txt" gives the info of the labels of nts,
     with a format:
     1-st line: k
     where k is the # of atomic propositions/ labels, of type int,
     the indexes are from 0 to k-1;
     for the next k lines, from 2 to k+1 line:
     each line represents the states with that label, ends with -1;
     the indexes of the states are from 0 to n0-1,
     these states are of type long long int;
     there is no requirement for the order of states for each label.
     */
    
    tt.tic();
    // dba->k=k=3;
    n0 = abst._x._nv;
    dba->labels=labels=(BOOL*)calloc(n0*k,sizeof(BOOL));
    for(int cnt1=0;cnt1<k;cnt1++)
        for(auto cnt2=0;cnt2<n0;cnt2++)
            if(target(cnt1,cnt2))*(labels+cnt2*k+cnt1)=1;
    printf("Time used to input labels = \n");
    tt.toc();
    printf("********\n\n");
    
    std::cout << "Computing the transition function: " << std::endl;
    /* transition function of symbolic model */
    double e1[] = {0,0,0};
    double e2[] = {0,0,0};
    tt.tic();
    abst.assign_transitions(e1, e2);
    tt.toc();
    
    
    // MARK: input nts in pres
    tt.tic();
    printf("Input nts in pres:\n");
    /*
     "ntsX.txt" gives the info of a non-deterministic transition system in its pres,
     with a format:
     1-st line: n0 action m0
     where n0 is the # of states, of type long long int,
     the indexes are from 0 to n0-1;
     action is the # of actions, of type int,
     the indexes are from 0 to action-1;
     m0 is the # of transitions, of type long long int;
     from 2-nd line to the end: each line has format x1 a x2,
     which represents an outgoing edge from state x1
     taking action a to state x2.
     All of the edges are sorted in ascending order w.r.t x2;
     There is no requirement for the order of x1 or a.
     */
    
    pre->n=n0=abst._ts._nx;
    head->action=action=abst._ts._nu;
    pre->m=m0=abst._ts._ntrans;
    pre->graph=nts_pre=(NODE_PRE*)calloc(n0,sizeof(NODE_PRE));
    pre->outdeg=outdeg=(int*)calloc(n0*action,sizeof(int));
    p5=pre->in_edge=in_edge=(EDGE*)malloc((m0+n0)*sizeof(EDGE));
    x=-1;
    for(auto i=0;i<n0;i++)
        for(auto j=0;j<action;j++)
            for(const auto & item : abst._ts.get_pre(i,j))
            {
                x1=item;a=j;x2=i;
                (*(outdeg+x1*action+a))++;
                if(x!=x2)
                {
                    if(x!=-1)p5++->next=0;
                    (nts_pre+x2)->in_edge=p5;x=x2;
                }
                p5++->next=1;p5->x=x1;p5->a=a;
            }
    p5->next=0;
    printf("Time used to input nts_pre = \n");
    tt.toc();
    printf("********\n\n");
}

int main() {
    TicToc tt;
    HEAD head,*p = &head;
    int dba_num[10]={6},cnt=1;
    p->nts_num = 3; p->labels_num = 2;
    
    initialization(p);
    input_labels_and_pre(p);
    optimize_pre(p);
    pre2post(p);
    for(int i=0;i<cnt;)
    {
        tt.tic();
        p->dba_num=dba_num[i++];
        post2product(p);
        buchi_and_controller(p);
        if(i!=cnt)
            free_memory(p,0);
        else
            free_memory(p,1);
        printf("Time used to solve dba%d = \n",
               dba_num[i-1]);
        tt.toc();
        printf("********\n\n");
    }
    
    
    printf("Total time used = %.4f\n",(double)clock()/CLOCKS_PER_SEC);
    
    return 1;
}

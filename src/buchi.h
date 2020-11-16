/**
 *  buchi.h
 *
 *  Headers for solving a Buchi game on a product system of an NTS and a DBA.
 *  Basic data structure for graph.
 *
 *  Created by Zhibing Sun on Mar. 27, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _buchi_h
#define _buchi_h


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/**
   This data structure stores the in-going edges of a state in nts_pre.
   It uses an array to implement the functions of a list,
   which enables us to delete useless edges easily.
   For an out-going edge x1->a->x2, we are at state x2, where
   EDGE.x == x1 is the in-going edge of x2 w.r.t. action a;
   EDGE.a == a is the action;
   EDGE.next is the # of steps we need to move to get to the next meaningful edge.
   The first edge only has EDGE.next assigned, the list ends with EDGE.next == 0.
   Aside: EDGE[0].x == indeg can be assigned as extra info.
*/
typedef struct EDGE {
    long long x;
    int a, next;
} EDGE;

/**
   This data structure stores the info of a state i,
   0 <= i <= n0-1 in nts_pre.
   NODE_PRE.num_a is the # of currently available actions;
   NODE_PRE.in_edge[m0_i+1] stores all the in-going edges of i,
   where m0_i is the in-degree of the i-th node in nts_pre.
   The first EDGE only has EDGE.next assigned, the list ends with EDGE.next == 0.
   Aside: NODE_PRE.in_edge[0].x == indeg can be assigned as extra info.
*/
typedef struct NODE_PRE {
    EDGE *in_edge;
    int num_a;
} NODE_PRE;

/**
   This data structure stores the info of nts_pre.
   For nts_pre, state indexes are from 0 to n0-1;
   NTS_PRE.graph[n0] is the set of states of nts_pre;
   NTS_PRE.outdeg[n0*action] stores the # of non-determinism of each action of each state;
   NTS_PRE.in_edge[m0+n0] stores all the in-going edges.
   NTS_PRE[i].in_edge[m0_i+1] stores all the in-going edges of i,
   where m0_i is the in-degree of the i-th node in nts_pre.
   NTS_PRE.n == n0 is the # of states;
   NTS_PRE.m == m0 is the # of transitions.
*/
typedef struct NTS_PRE {
    NODE_PRE *graph;
    int *outdeg;
    EDGE *in_edge;
    long long n,m;
} NTS_PRE;

/**
   This data structure stores the info of a state i,
   0 <= i <= n0_1-1 in nts_post.
   NODE_POST.num_a is the # of currently available actions;
   NODE_POST.label is the label of the state in decimal representation;
   NODE_POST.pos is the position of the first available action in NTS_POST.a_index and in NTS_POST.a_pos.
*/
typedef struct NODE_POST {
    int num_a,label;
    long long pos;
} NODE_POST;

/**
   This data structure stores the info of nts_post.
   For nts_post, state indexes are from 0 to n0_1-1;
   NTS_POST.graph[n0_1] is the set of states of nts_post;
   NTS_POST.a_index[m_a] stores the available actions of each state in ascending order;
   NTS_POST.n == n0_1 is the # of states;
   NTS_POST.m_a == m_a is the sum of all the available actions of all the states;
   NTS_POST.m_e == m_e is the # of transitions;
   NTS_POST.a_pos[m_a+1] stores the position of the first non-determinism of each available action of each state in NTS_POST.out_edge;
   NTS_POST.a_pos[0] == 0,
   NTS_POST.a_pos[m_a] == m_e;
   NTS_POST.out_edge[m_e] stores all the out-going edges,
   (*)    the non-determinism of each action is in ascending order.
*/
typedef struct NTS_POST {
    NODE_POST *graph;
    int *a_index;
    long long n,m_a,m_e,*a_pos,*edge;
} NTS_POST;

/**
   BOOL is a data type to record TRUE VALUES.
   Each unit is 1 byte = 8 bits, so it can record more than just TRUE == 1
   and FALSE == 0, but actually from -128 to 127.
*/
typedef char BOOL;

/**
   This data structure stores the info of a dba and the labels of nts_pre in binomial representation.
   For a dba, the indexes are from 0 to n1-1;
   DBA.n == n1 is the # of states;
   DBA.m == m1 is the # of transitions;
   DBA.ini is the unique initial state;
   DBA.k is the # of atomic propositions/ labels,
   where the indexes are from 0 to k-1;
   DBA.q_prime_size == 2^k*n1;
   DBA.q_prime[2^k][n1] is the look up table for q_prime,
   where the values are from -1 to n1:
   -1 is the initialized value, 0 to n1-1 are the live states,
   n1 is the non-live state, if exists;
   DBA.labels[n0][k] stores the labels of nts_pre in binomial representation, where the values are 0 or 1;
   DBA.acc[k] marks the accepting states, where the values are 0 or 1, we can have more than 1 accepting states which are marked as 1.
*/
typedef struct DBA {
    int n,m,ini,k,q_prime_size,*q_prime;
    BOOL *labels; // not used
    BOOL *acc;
} DBA;

/**
   This data structure stores the info to reset an out-degree of an action in nts_product.
   It uses an array to implement the functions of a list,
   which enables us to delete useless nodes easily.
   For an action a of a state (x1,q) with out-degree outdeg, we are at state (x1,q), where
   OUTDEG_RESET.a == a is the action;
   OUTDEG_RESET.outdeg == outdeg is the out-degree of action a of state (x1,q);
   OUTDEG_RESET.next is the # of steps we need to move to get to the next meaningful action.
   The first OUTDEG_RESET only has OUTDEG_RESET.next assigned, the list ends with OUTDEG_RESET.next == 0.
*/
typedef struct OUTDEG_RESET {
    int a,outdeg,next;
} OUTDEG_RESET;

/**
   This data structure stores the info of a state i,
   0 <= i <= np-1 in nts_product.
   NODE_PRODUCT.num_a is the # of currently available actions;
   NODE_PRODUCT.iterate is the current iteration # of the loops when solving buchi;
   NODE_PRODUCT.in_edge[mp_i_in+1] stores all the in-going edges of i,
   where mp_i_in is the in-degree of the i-th node in nts_product.
   The first EDGE only has EDGE.next assigned, the list ends with EDGE.next == 0.
   Aside: NODE_PRE.in_edge[0].x == indeg can be assigned as extra info.
   NODE_PRODUCT.outdeg_reset[num_a+1] stores the info to reset the out-degree of all the actions of i,
   The first OUTDEG_RESET only has OUTDEG_RESET.next assigned, the list ends with OUTDEG_RESET.next == 0.
*/
typedef struct NODE_PRODUCT {
    int num_a,iterate;
    EDGE *in_edge;
    OUTDEG_RESET *outdeg_reset;
} NODE_PRODUCT;

/**
   This data structure stores the scope of a certain property in nts_product.
   It uses an array to implement the functions of a list,
   which enables us to delete useless nodes easily.
   SCOPE.x is the state;
   SCOPE.next is the # of steps we need to move to get to the next meaningful state.
   The first SCOPE only has SCOPE.next assigned, the list ends with SCOPE.next == 0.
   Aside: SCOPE[0].x == # of states remaining in SCOPE can be assigned as extra info.
*/
typedef struct SCOPE {
    long long x, next;
} SCOPE;

/**
   This data structure stores the info of nts_product.
   For nts_product, state indexes are from 0 to np-1;
   NTS_PRODUCT.graph[np] is the set of states of nts_product;
   NTS_PRODUCT.outdeg[np*action] stores the # of non-determinism of each action of each state;
   NTS_PRODUCT.in_edge[mp_e+n.(indeg>0)] stores all the in-going edges.
   NTS_PRODUCT[i].in_edge[mp_i_in+1] stores all the in-going edges of i,
   where mp_i_in is the in-degree of the i-th node in nts_product.
   NTS_PRODUCT.outdeg_reset[mp_a+n.(outdeg/num_a >0)] stores the info to reset the out-degree of all the actions.
   NTS_PROEUCT[i].outdeg_reset[num_a+1] stores the info to reset the out-degree of all the actions of i.
   NTS_PRODUCT.acc[np*1] marks the accepting states, where the values are 0 or 1.
   NTS_PRODUCT.scope[wp+1] stores the scope of the current game graph.
   Aside: NTS_PRODUCT.scope[0].x == # of states remaining in game graph can be assigned as extra info.
   NTS_PRODUCT.target[wp+1] stores the scope of the current target.
   Aside: NTS_PRODUCT.target[0].x == # of states remaining in target can be assigned as extra info.
   NTS_PRODUCT.n == np is the # of states;
   NTS_PRODUCT.m_a == mp_a is the sum of all the available actions of all the states;
   NTS_PRODUCT.m_e == mp_e is the # of transitions;
   NTS_PRODUCT.w0 is the # of base in nts_post to take product;
   NTS_PRODUCT.w0_1 is the # of winning states of W_X0 in nts;
   NTS_PRODUCT.wp is the # of winning states of W_XP in nts_product;
   NTS_PRODUCT.queue[np] is the queue to do BFS for reachability and safety to solve buchi game in product space.
*/
typedef struct NTS_PRODUCT {
    NODE_PRODUCT *graph;
    int *outdeg;
    EDGE *in_edge;
    OUTDEG_RESET *outdeg_reset;
    BOOL *acc;
    SCOPE *scope,*target;
    long long n,m_a,m_e,w0,w0_1,wp,*queue;
} NTS_PRODUCT;

/**
   This data structure stores the control of a state (x,q) in W_XP.
   x is implicitly indicated from CTRLR.graph.
   CTRL.q == q is the state of automaton corresponding to the q in (x,q);
   CTRL.u is the control to take for state (x,q), which is an action.
*/
typedef struct CTRL {
    int q, u;
} CTRL;

/**
   This data structure stores the info of ctrlr.
   For ctrlr, state indexes are from 0 to n0_2-1;
   CTRLR.graph[n0_2] == nts_ctrlr2 is the set of states that projects from W_XP to X;
   For a state x, nts_ctrlr2[x].num_a is the # of states (x,q) in W_XP that corresponds to x;
   nts_ctrlr2[x].label is the label of the state x in decimal representation;
   nts_ctrlr2[x].pos is the position of the first state (x,q) that corresponds to x in ctrl;
   CTRLR.ctrl[wp] is the set of controls for each state in W_XP, which is sorted in multi-steps, where x is sorted implicitly in ascending order, and then q is sorted explicitly in ascending order;
   ctrl.n == n0_2 is the # of states that projects from W_XP to X;
   CTRLR.w0 is the # of winning states of W_X0 in nts;
   CTRLR.wp_acc is the # of target states in W_XP in nts_product;
   CTRLR.wp is the # of winning states of W_XP in nts_product.
*/
typedef struct CTRLR {
    NODE_POST *graph;
    CTRL *ctrl;
    long long n,w0,wp_acc,wp;
} CTRLR;

/**
   This data structure is the head of all the other data structures. It helps to transmit info of all the other data structures between functions.
   HEAD.action == action is the # of actions;
   HEAD.nts_num is the version of nts;
   HEAD.labels_num is the version of labels;
   HEAD.dba_num is the version of dba;
   HEAD.encode1[n0] encodes the states from nts_pre to nts_post, states that do not have a mapping is marked as -1;
   HEAD.decode1[n0_1] decodes the states from nts_post to nts_pre;
   HEAD.encode2[n0_1*n1] encodes the states from nts_post to nts_product, states that do not have a mapping is marked as -1;
   HEAD.decode2[np] decodes the states from nts_product to nts_post;
   HEAD.encode3[n0] encodes the states from nts_pre to nts_ctrlr2, states that do not have a mapping is marked as -1;
   HEAD.nts_pre is the sub-head for nts in pres;
   HEAD.nts_post is the sub-head for nts in posts;
   HEAD.dba is the sub-head for dba;
   HEAD.nts_product is the sub-head for the product space of nts in pres;
   HEAD.ctrlr is the sub-head for the controller.
*/
typedef struct HEAD {
    int action,nts_num,labels_num,dba_num;
    long long *encode1,*decode1,*encode2,*decode2,*encode3;
    NTS_PRE nts_pre;
    NTS_POST nts_post;
    DBA dba;
    NTS_PRODUCT nts_product;
    CTRLR ctrlr;
} HEAD;
    

/**
   This compare funciton can be used to do a bsearch to find the action when converting nts_pre into nts_post, but I didn't use it in the end.
*/
int cmp1(const void *a,const void *b);

/**
   This is the compare function for qsort to sort the q's of the same x in ctrl in ascending order. This compare function is itself not a multi-step compare function since we have already sort the x in ascending order naturally by counting the right num_a and pos in CTRLR.graph[n0_2] == nts_ctrlr2.
   The same funciton can be used to do a bsearch to find the control when doing control synthesis.
   Since n1 is small, in either case the time complexity can be considered as a small const.
*/
int cmp2(const void *a,const void *b);

void initialization(HEAD *head);

/**
 * Perform a safety preprocessing to delete transitions to/from avoid states.
 */
void safety_pre(HEAD *head);

/**
 * Convert a reverse graph back to a forward graph (post graph), and
 * assign labels to nts_post by the state labels of an abstraction.
 */
void pre2post(HEAD *head, int *labels);

/**
 * Construct a product system of a post graph and a DBA.
 */
void post2product(HEAD *head);

/**
 * Solve a buchi game on the product system and output a controller.
 */
void buchi_and_controller(HEAD *head);

/**
 * Write the synthesized controller into a .txt file. 
 */
void write_controller(HEAD *head, char *buffer);
    
/**
 * Release memory.
 * flag = 1: release all memories
 * flag = 0: release the memory related to DBA. (Use 0 when multiple DBAs are using in a single run)
 */
void free_memory(HEAD *head,int flag);


#endif

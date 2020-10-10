/**
 *  carAbst_ref.h
 *  The headers for discrete DBA control synthesis.
 *
 *  This is a reference file for testing.
 *
 *  Created by Zhibing Sun on May 26, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _abst_reference_h
#define _abst_reference_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/**
 This data structure is used to store the in-going edges of a state in nts_pre.
 It uses an array to implement the functions of a list,
 which enables us to delete useless edges easily.
 For an out-going edge x1-a->x2, we are at state x2,
 edge.a == a; edge.x == x1;
 edge.next is the # of steps we need to move to get to the next meaningful edge.
 The first edge only has edge.next assigned, the list ends with edge.next == 0.
 Aside: edge[0].x == indeg can be extra info.
 */
typedef struct EDGE
{long long x;int a,next;}EDGE;

/**
 This data structure is used to store the info of a state i,
 0 <= i <= n0-1 in nts_pre.
 node_pre.outdeg[action] is used to store the # of non-determinism of each action;
 node_pre.num is the # of currently available actions;
 node_pre.in_edge[m0i+1] stores all the in-going edges of i,
 where m0i is the in-degree of the i-th node in nts.
 The first edge only has edge.next assigned, the list ends with edge.next == 0.
 Aside: node_pre.in_edge[0].x == indeg can be extra info.
 */
typedef struct NODE_PRE
{EDGE *in_edge;int num_a;}NODE_PRE;

/**
 This data structure is used to store the info of nts_pre.
 For nts_pre, state indexes are from 0 to n0-1;
 nts_pre.graph[n0] is the set of states of nts_pre;
 nts_pre.n == n0 is the # of states;
 nts_pre.m == m0 is the # of edges.
 */
typedef struct NTS_PRE
{
    NODE_PRE *graph;
    int *outdeg;
    EDGE *in_edge;
    long long n,m;
    
}NTS_PRE;

/**
 This data structure is used to store the info of a state i,
 0 <= i <= n0_1-1 in nts_post.
 node_post.num_a is the # of currently available actions;
 node_post.label is the label of the state in decimal representation;
 node_post.pos is the pointer of the first action in nts_post.a_index.
 */
typedef struct NODE_POST
{int num_a,label;long long pos;}NODE_POST;

/**
 This data structure is used to store the info of nts_post.
 For nts_post, state indexes are from 0 to n0_1-1;
 nts_post.graph[n0_1] is the set of states;
 nts_post.a_index[nts_post.m_a] stores the available actions of each state in ascending order;
 nts_post.n == n0_1 is the # of states;
 nts_post.m_a is the sum of all the available actions of all the states;
 nts_post.m_e is the # of edges;
 nts_post.a_pos[nts_post.m_a+1] stores the position of the first non-determinism
 of each available action of each state, nts_post.a_pos[0] == 0,
 nts_post.a_pos[nts_post.m_a] == nts_post.m_e;
 nts_post.edge[nts_post.m_e] stores all the out-going edges,
 (*)    the non-determinism of each action is in ascending order.
 */
typedef struct NTS_POST
{
    NODE_POST *graph;
    int *a_index;
    long long n,m_a,m_e,*a_pos,*edge;
}NTS_POST;

/**
 BOOL is a data type to record TRUE VALUES.
 Each unit is 1 byte = 8 bits, so it can record more than just TRUE == 1
 and FALSE == 0, but actually from -128 to 127.
 */
typedef char BOOL;

/**
 This data structure is used to store the info of a dba and the labels of nts.
 For a dba, the indexes are from 0 to n1-1;
 dba.n == n1 is the # of states;
 dba.m == m1 is the # of edges;
 dba.ini is the unique initial state;
 dba.k is the # of atomic propositions/ labels,
 where the indexes are from 0 to k-1;
 q_prime_size == 2^k*n1;
 dba.q_prime[2^k][n1] are the look up table for q_prime,
 where the values are from -1 to n1,
 -1 is the initialized value, 0 to n1-1 are the live states,
 n1 is the non-live state, if exist;
 dba.labels[n0][k] are the labels of nts;
 dba.acc[k] marks the accepting states, can be more than 1.
 */
typedef struct DBA
{
    int n,m,ini,k,q_prime_size,*q_prime;
    BOOL *labels,*acc;
}DBA;

typedef struct OUTDEG_RESET
{int a,outdeg,next;}OUTDEG_RESET;

typedef struct NODE_PRODUCT
{
    int num_a,iterate;
    EDGE *in_edge;
    OUTDEG_RESET *outdeg_reset;
}NODE_PRODUCT;

typedef struct SCOPE
{long long x,next;}SCOPE;

/**
 nts_product.graph[np]
 nts_product.indeg[n0_1*n1]
 nts_product.outdeg[np*action]
 nts_product.in_edge[nts_product.m_e+n.(indeg>0)]
 nts_product.outdeg_reset[nts_product.m_a+n.(outdeg/num_a>0)]
 nts_product.acc[np*1]
 */
typedef struct NTS_PRODUCT
{
    NODE_PRODUCT *graph;
    int *outdeg;
    EDGE *in_edge;
    OUTDEG_RESET *outdeg_reset;
    BOOL *acc;
    SCOPE *scope,*target;
    long long n,m_a,m_e,w0,w0_1,wp,*queue;
}NTS_PRODUCT;

typedef struct CTRL
{int q,u;}CTRL;

/**
 ctrl.n == n0_2;
 graph[n0_2];
 w0 == # of states in W_X0;
 wp == # of states in W_P;
 ctrl[wp];
 
 */
typedef struct CTRLR
{
    NODE_POST *graph;
    CTRL *ctrl;
    long long n,w0,wp_acc,wp;
}CTRLR;

typedef struct HEAD
{
    int action,nts_num,labels_num,dba_num;
    long long *encode1,*decode1,*encode2,*decode2,*encode3;
    NTS_PRE nts_pre;
    NTS_POST nts_post;
    DBA dba;
    NTS_PRODUCT nts_product;
    CTRLR ctrlr;
}HEAD;

int cmp1(const void *a,const void *b)
{return *(int*)a-*(int*)b;}

int cmp2(const void *a,const void *b)
{
    return ((CTRL*)a)->q-((CTRL*)b)->q;
    //    CTRL *p1=(CTRL*)a,*p2=(CTRL*)b;
    //    return p1->q-p2->q;
}


/**
 p: q_prime+label*n1+q1
 a: labels+*(decode+x1)*k
 */
void q_next(HEAD *head,int *p,int q,BOOL *a)
{
    switch(head->dba_num)
    {
        case 0:
        {
            //  Fa
            switch(q)
            {
                case 0:
                {
                    *p=0;  // 1
                    break;
                }
                case 1:
                {
                    if(!*a)*p=1;  // !a
                    else *p=0; // a
                }
            }
            break;
        }
        case 1:
        {
            //  GFa && GFb && GFc (Spot's)
            switch(q)
            {
                case 0:
                {
                    if(*a&&*(a+1)&&*(a+2))*p=0;
                    else if(!*(a+2))*p=1;
                    else if(!*(a+1)&&*(a+2))*p=2;
                    else *p=3;
                    break;
                }
                case 1:
                {
                    if(*a&&*(a+1)&&*(a+2))*p=0;
                    else if(!*(a+2))*p=1;
                    else if(!*(a+1)&&*(a+2))*p=2;
                    else *p=3;
                    break;
                }
                case 2:
                {
                    if(*a&&*(a+1))*p=0;
                    else if(!*(a+1))*p=2;
                    else *p=3;
                    break;
                }
                case 3:
                {
                    if(*a)*p=0;
                    else *p=3;
                }
            }
            break;
        }
        case 2:
        {
            //  F(a && F(b && F(c && (!b U a))))  (Spot's)
            switch(q)
            {
                case 0:
                {
                    if(!*a&&!*(a+1))*p=0;
                    else if(!*a&&*(a+1))*p=1;
                    else *p=4;
                    break;
                }
                case 1:
                {
                    if(!*a&&!*(a+1)&&*(a+2))*p=0;
                    else if((!*a&&*(a+1))||!*(a+2))*p=1;
                    else *p=4;
                    break;
                }
                case 2:
                {
                    if(*a&&*(a+1)&&!*(a+2))*p=1;
                    else if(!*a)*p=2;
                    else if(*a&&!*(a+1))*p=3;
                    else *p=4;
                    break;
                }
                case 3:
                {
                    if((!*a&&*(a+1))||(*(a+1)&&!*(a+2)))*p=1;
                    else if(!*(a+1))*p=3;
                    else *p=4;
                    break;
                }
                case 4:
                {
                    *p=4;
                }
            }
            break;
        }
        case 3:
        {
            //  GFa && GFb && GFc   (Yinan's)
            switch(q)
            {
                case 0:
                {
                    if(!*(a+2))*p=1;
                    else *p=2;
                    break;
                }
                case 1:
                {
                    if(!*(a+2))*p=1;
                    else *p=2;
                    break;
                }
                case 2:
                {
                    if(!*(a+1))*p=2;
                    else *p=3;
                    break;
                }
                case 3:
                {
                    if(*a)*p=0;
                    else *p=3;
                }
            }
            break;
        }
        case 4:
        {
            //  F(a && F(b && F(c && (!b U a))))    (Yinan's)
            //  0 1 2 3 4   permutation from case 4 to case 5
            //  2 3 1 0 4
            switch(q)
            {
                case 0:
                {
                    if(!*a)*p=0;
                    else *p=1;
                    break;
                }
                case 1:
                {
                    if(!*(a+1))*p=1;
                    else *p=2;
                    break;
                }
                case 2:
                {
                    if(!*(a+2))*p=2;
                    else *p=3;
                    break;
                }
                case 3:
                {
                    if(*(a+1))*p=2;
                    else if(!*a&&!*(a+1))*p=3;
                    else *p=4;
                    break;
                }
                case 4:
                {
                    *p=4;
                }
            }
            break;
        }
        case 5:
        {
            //  F(a && F(b && F(c && (!b U a))))    (Yinan's)
            switch(q)
            {
                case 0:
                {
                    if(*(a+1))*p=1;
                    else if(!*a&&!*(a+1))*p=0;
                    else *p=4;
                    break;
                }
                case 1:
                {
                    if(!*(a+2))*p=1;
                    else *p=0;
                    break;
                }
                case 2:
                {
                    if(!*a)*p=2;
                    else *p=3;
                    break;
                }
                case 3:
                {
                    if(!*(a+1))*p=3;
                    else *p=1;
                    break;
                }
                case 4:
                {
                    *p=4;
                }
            }
            break;
        }
        case 6:
        {
            //  (F(a1 & Fa3) | F(a2 & (!a4 U a3))) & (GFc & GFd)
            switch(q)
            {
                case 0:
                {
                    if((!a[0]&&!a[1]&&a[5])||(!a[0]&&!a[4]))*p=0;
                    else if(!a[0]&&!a[1]&&a[4]&&!a[5])*p=1;
                    else if((a[0]&&a[1]&&a[2]&&a[3])||(a[1]&&a[2]&&a[3]&&a[4]))*p=2;
                    else if((a[0]&&a[1]&&!a[3])||(a[1]&&!a[3]&&a[4]))*p=3;
                    else if((a[0]&&a[1]&&!a[2]&&a[3])||(a[1]&&!a[2]&&a[3]&&a[4]))*p=4;
                    else *p=5;
                    break;
                }
                case 1:
                {
                    if(!a[0]&&!a[1]&&a[5])*p=0;
                    else if(!a[0]&&!a[1]&&!a[5])*p=1;
                    else if(a[1]&&a[2]&&a[3])*p=2;
                    else if(a[1]&&!a[3])*p=3;
                    else if(a[1]&&!a[2]&&a[3])*p=4;
                    else *p=5;
                    break;
                }
                case 2:
                {
                    if(a[2]&&a[3])*p=2;
                    else if(!a[3])*p=3;
                    else *p=4;
                    break;
                }
                case 3:
                {
                    if(a[2]&&a[3])*p=2;
                    else if(!a[3])*p=3;
                    else *p=4;
                    break;
                }
                case 4:
                {
                    if(a[2])*p=2;
                    else *p=4;
                    break;
                }
                case 5:
                {
                    if(a[1]&&a[2]&&a[3])*p=2;
                    else if(a[1]&&!a[3])*p=3;
                    else if(a[1]&&!a[2]&&a[3])*p=4;
                    else *p=5;
                }
            }
            break;
        }
    }
}


void initialization(HEAD *head)
{
    NTS_PRE *pre;
    NTS_POST *post;
    DBA *dba;
    NTS_PRODUCT *product;
    CTRLR *ctrlr;
    clock_t start_t,end_t;
    
    pre=&head->nts_pre;
    post=&head->nts_post;
    dba=&head->dba;
    product=&head->nts_product;
    ctrlr=&head->ctrlr;
    
    // MARK: initialize variables to 0
    start_t=clock();
    printf("Initialize variables to 0:\n");
    /// Initialize head
    head->encode1=head->decode1=head->encode2=head->decode2=head->encode3=0;
    
    /// Initialize nts_pre
    pre->outdeg=0;pre->in_edge=0;pre->graph=0;
    
    /// Initialize nts_post
    post->a_index=0;post->a_pos=post->edge=0;post->graph=0;
    
    /// Initialize dba
    dba->q_prime=0;dba->labels=dba->acc=0;
    
    /// Initialize nts_product
    product->outdeg=0;product->in_edge=0;product->outdeg_reset=0;
    product->acc=0;product->scope=product->target=0;
    product->queue=0;product->wp=0;product->graph=0;
    
    /// Initialize ctrlr
    ctrlr->ctrl=0;ctrlr->graph=0;
    end_t=clock();
    printf("Time used to initialize variables to 0 = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    printf("********\n\n");
}

void optimize_pre(HEAD *head)
{
    NTS_PRE *pre;
    NODE_PRE *nts_pre,*p1,*p2;
    NTS_POST *post;
    int action,*outdeg;
    long long cnt;
    long long n0,n0_1,m0,m0_1,*p3,*p4;
    long long *queue,*encode1,*decode1;
    EDGE *in_edge,*p5;
    int *p7,*p8;
    clock_t start_t0,start_t,end_t;
    
    action=head->action;
    pre=&head->nts_pre;nts_pre=pre->graph;n0=pre->n;m0=pre->m;
    outdeg=pre->outdeg;in_edge=pre->in_edge;
    post=&head->nts_post;
    printf("optimize nts_pre:\n");
    printf("# of states n0: %lld\n",n0);
    printf("# of transitions m0: %lld\n",m0);
    printf("Average degree of nts_pre m0/n0 = %lld\n",m0/n0);
    p3=queue=(long long*)malloc(n0*sizeof(long long));
    
    // MARK: do safety in nts_pre
    start_t0=start_t=clock();
    printf("Do safety to guarantee no terminal states:\n");
    ///  Set up base for safety
    for(p1=nts_pre,p2=p1+n0,p3=p4=queue,n0_1=n0,p7=outdeg;p1<p2;p1++)
    {
        for(p8=p7+action;p7<p8;)
            if(*p7++)p1->num_a++;
        if(!p1->num_a)
        {if(p1->in_edge)*p4++=p1-nts_pre;else n0_1--;}
    }
    
    ///  BFS for safety
    for(m0_1=m0;p3<p4;p3++)
        for(p5=(nts_pre+*p3)->in_edge;p5->next;)
        {
            p5+=p5->next;p1=nts_pre+p5->x;
            if(p1->num_a)
            {
                p7=outdeg+p5->x*action+p5->a;
                if(*p7)
                {
                    m0_1-=*p7;*p7=0;p1->num_a--;
                    if(!p1->num_a)
                    {if(p1->in_edge)*p4++=p5->x;else n0_1--;}
                }
            }
        }
    post->n=n0_1-=p4-queue;post->m_e=m0_1;
    end_t=clock();
    printf("Time used to do safety = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    if(!n0_1)
    {
        printf("n0_1 is empty.\n");
        free(queue);
        return;
    }
    // MARK: compression mapping1
    start_t=clock();
    printf("Compression mapping1:\n");
    p3=head->encode1=encode1=queue;
    p4=head->decode1=decode1=(long long*)malloc(n0_1*sizeof(long long));
    for(p1=nts_pre,p2=p1+n0,cnt=0;p1<p2;p1++,p3++)
        if(p1->num_a){cnt+=p1->num_a;*p3=p4-decode1;*p4++=p3-encode1;}
        else *p3=-1;
    post->m_a=cnt;
    end_t=clock();
    printf("Time used to do compression mapping1 = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    printf("Total time used to optimize nts_pre = %.4f\n",
           (double)(end_t-start_t0)/CLOCKS_PER_SEC);
    printf("********\n\n");
}

void pre2post(HEAD *head)
{
    NTS_PRE *pre;
    NODE_PRE *nts_pre,*p1;
    int action,*outdeg,*p3;
    EDGE *p4;
    NTS_POST *post;
    NODE_POST *nts_post,*p2;
    int *a_index,*p5;
    long long *a_pos,*p6,*edge;
    int k,cnt1,cnt2,cnt3;
    long long x,n0,n0_1,m_a,m_e;
    long long *encode1,*decode1,*p7,*p8;
    DBA *dba;
    BOOL *labels,*p9;
    clock_t start_t0,start_t,end_t;
    
    if(!head->nts_post.n)return;
    action=head->action;
    pre=&head->nts_pre;nts_pre=pre->graph;outdeg=pre->outdeg;n0=pre->n;
    post=&head->nts_post;n0_1=post->n;m_a=post->m_a;m_e=post->m_e;
    p7=decode1=head->decode1;
    encode1=head->encode1;
    dba=&head->dba;k=dba->k;labels=dba->labels;
    printf("pre2post:\n");
    
    /// labels are inputted at the very beginning in source file
    
    printf("# of states n0_1: %lld\n",n0_1);
    printf("# of transitions m0_1: %lld\n",m_e);
    printf("Average degree of nts_post m_e/n0_1 = %lld\n",m_e/n0_1);
    printf("Average degree of action m_a/n0_1 = %lld\n",m_a/n0_1);
    printf("Average degree of non-determinism m_e/m_a = %lld\n",m_e/m_a);
    
    // MARK: Set all the parameters in NODE_POST and NTS_POST
    start_t0=start_t=clock();
    printf("Set all the parameters in node_post and nts_post:\n");
    p2=post->graph=nts_post=(NODE_POST*)malloc(n0_1*sizeof(NODE_POST));
    p5=post->a_index=a_index=(int*)malloc(m_a*sizeof(int));
    p6=post->a_pos=a_pos=(long long*)malloc((m_a+1)*sizeof(long long));
    post->edge=edge=(long long*)malloc(m_e*sizeof(long long));
    for(p8=p7+n0_1,*p6++=x=0;p7<p8;p7++,p2++)
    {
        p1=nts_pre+*p7;p9=labels+*p7*k;
        
        // MARK: - convert labels from binomial to decimal representation
        for(cnt1=k,cnt2=*p9++;--cnt1;)
        {cnt2*=2;cnt2+=*p9++;}
        
        // MARK: set label, num_a and pos in node_post
        /// Here I set the pos of each action one unit backward so that after the conversion from pre to post it will just match up.
        p2->label=cnt2;cnt1=p2->num_a=p1->num_a;p2->pos=p5-a_index;
        
        // MARK: set parameters in a_index and a_pos in nts_post
        /// Here I set the pos of each action one unit forward so that when converting from pre to post it will just cancel the miss match.
        for(p3=outdeg+*p7*action,cnt2=cnt3=0;cnt2<cnt1;p3++,cnt3++)
            if(*p3){*p5++=cnt3;*p6++=x;x+=*p3;*p3=++cnt2;}
    }
    end_t=clock();
    printf("Time used to set all the parameters in node_post and nts_post = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: - convert edges into post after encoding
    start_t=clock();
    printf("Convert edges into post after encoding:\n");
    for(p7=decode1;p7<p8;p7++)
    {
        p4=(nts_pre+*p7)->in_edge;
        if(p4)
        {
            for(;p4->next;)
            {
                p4+=p4->next;p1=nts_pre+p4->x;
                p3=outdeg+p4->x*action+p4->a;
                if(p1->num_a&&*p3)
                {
                    p2=nts_post+*(encode1+p4->x);
                    p6=a_pos+p2->pos+*p3;
                    *(edge+(*p6)++)=p7-decode1;
                }
            }
        }
    }
    end_t=clock();
    printf("Time used to convert edges into post after encoding = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: free nts_pre
    start_t=clock();
    printf("Free nts_pre:\n");
    free(outdeg);pre->outdeg=0;
    free(pre->in_edge);pre->in_edge=0;
    free(nts_pre);pre->graph=0;
    end_t=clock();
    printf("Time used to free nts_pre = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    printf("Total time used to convert pre to post: = %.4f\n",
           (double)(end_t-start_t0)/CLOCKS_PER_SEC);
    printf("********\n\n");
}

void post2product(HEAD *head)
{
    NTS_POST *post;
    NODE_POST *nts_post,*p1;
    NTS_PRODUCT *product;
    NODE_PRODUCT *nts_product,*p3,*p4;
    int *indeg,*outdeg,*q_prime,*a_index,action,n1,ini,k,q1,q2;
    long long *queue,*a_pos,*edge,*decode1,*encode2,*decode2;
    EDGE *in_edge,*p5;
    OUTDEG_RESET *outdeg_reset,*p6;
    DBA *dba;
    BOOL *acc_dba,*acc_p,*labels,*mark,*p7;
    char buffer[40];
    int *p8,*p9,*p10;
    long long x1,xp1,xp2,n0_1,m0_1,n,np,mp_a,mp_e,w0,wp,cnt1,cnt2,i,j;
    long long *p11,*p12,*p13,*p14,*p15,*p16;
    FILE *fin;
    clock_t start_t0,start_t,end_t;
    
    if(!head->nts_post.n)return;
    action=head->action;
    post=&head->nts_post;nts_post=post->graph;
    n0_1=post->n;m0_1=post->m_e;
    a_index=post->a_index;a_pos=post->a_pos;edge=post->edge;
    product=&head->nts_product;
    dba=&head->dba;k=dba->k;
    labels=dba->labels;
    decode1=head->decode1;
    printf("post2product:\n");
    
    // MARK: input dba
    start_t=clock();
    printf("Input dba:\n");
    /*
     "dbaX.txt" gives the info of the dba,
     with a format:
     1-st line: n1 m1 ini
     where n1 is the # of states, of type int,
     the indexes are from 0 to n1-1;
     m1 is the # of transitions, of type int;
     ini is the index of the unique initial state, of type int;
     2-nd line: n1 numbers of 0 or 1, 1 for the marks of accs.
     */
    
    sprintf(buffer,"dba%d.txt",head->dba_num);
    printf("Read dba from %s:\n",buffer);
    fin=fopen(buffer,"rb");
    fscanf(fin,"%d%d%d",&n1,&dba->m,&ini);
    dba->n=n1;dba->ini=ini;
    p7=dba->acc=acc_dba=(BOOL*)malloc(n1*sizeof(BOOL));
    for(i=0;i<n1;i++)
    {fscanf(fin,"%d",&q1);*p7++=q1;}
    fclose(fin);
    
    //    for(p7=acc_dba,i=0;i<n1;i++)
    //    {
    //        printf("%d\t",*p7++);
    //    }
    //    putchar('\n');
    end_t=clock();
    printf("Time used to input dba = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: take product 1-st time, only count
    start_t0=start_t=clock();
    printf("Take product 1-st time, only count:\n");
    n=n0_1*n1;
    queue=(long long*)malloc(n*sizeof(long long));
    indeg=(int*)calloc(n,sizeof(int));
    mark=(BOOL*)calloc(n,sizeof(BOOL));
    cnt1=dba->q_prime_size=(int)pow(2,k)*n1;
    printf("2^%d*%d = %lld\n",k,n1,cnt1);
    dba->q_prime=q_prime=(int*)malloc(cnt1*sizeof(int));
    
    /// Initialize q_prime table to -1
    for(p9=q_prime,p10=p9+cnt1;p9<p10;)*p9++=-1;
    
    //    p9=q_prime+(nts_post+x)->label*n1+q;
    //    if(*p9==-1)q_next(head,p9,q,labels+*(decode1+x)*k);
    //    if(*p9!=n1)
    //    {
    //          *p9 is live
    //    }
    
    // MARK: set up base to take product 1-st time
    for(p1=nts_post,x1=0,p12=queue;x1<n0_1;p1++,x1++)
    {
        p9=q_prime+p1->label*n1+ini;
        if(*p9==-1)q_next(head,p9,ini,labels+*(decode1+x1)*k);
        if(*p9!=n1)
        {
            xp1=x1*n1+ini;
            *p12++=xp1;
            *(mark+xp1)=1;
        }
    }
    product->w0=w0=p12-queue;
    printf("Base in nts_post to take product w0 = %lld\n",w0);
    
    // MARK: take product 1-st time
    for(p11=queue,mp_a=cnt1=0;p11<p12;p11++)
    {
        xp1=*p11;x1=xp1/n1;q1=xp1%n1;
        p1=nts_post+x1;
        q2=*(q_prime+p1->label*n1+q1);
        for(i=p1->num_a,j=0,p13=a_pos+p1->pos;i--;p13++)
        {
            for(p14=edge+*p13,p15=edge+*(p13+1);p14<p15;p14++)
            {
                p9=q_prime+(nts_post+*p14)->label*n1+q2;
                if(*p9==-1)q_next(head,p9,q2,labels+*(decode1+*p14)*k);
                if(*p9==n1)break;
            }
            if(p14==p15)
                for(j++,p14=edge+*p13;p14<p15;p14++)
                {
                    xp2=*p14*n1+q2;
                    (*(indeg+xp2))++;
                    p7=mark+xp2;
                    if(!*p7){*p12++=xp2;*p7=1;}
                }
        }
        if(j){mp_a+=j;cnt1++;}
    }
    free(mark);
    product->n=wp=np=p12-queue;product->m_a=mp_a;
    end_t=clock();
    printf("Time used to take product 1-st time, only count = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    if(!np)
    {
        printf("np is empty.\n");
        free(queue);free(indeg);
        return;
    }
    // MARK: compression mapping2 and set up data
    start_t=clock();
    printf("Compression mapping2 and set up data:\n");
    p11=head->encode2=encode2=queue;
    p13=head->decode2=decode2=(long long*)malloc(np*sizeof(long long));
    p3=product->graph=nts_product=
    (NODE_PRODUCT*)calloc(np,sizeof(NODE_PRODUCT));
    p6=product->outdeg_reset=outdeg_reset=
    (OUTDEG_RESET*)malloc((mp_a+cnt1)*sizeof(OUTDEG_RESET));
    printf("n = n0_1*n1 = %lld*%d = %lld\n",n0_1,n1,n);
    printf("m = m0_1*n1 = %lld*%d = %lld\n",m0_1,n1,m0_1*n1);
    printf("np = %lld\n",np);
    printf("mp_a = %lld\n",mp_a);
    printf("# of product states with outdeg > 0: %lld\n",cnt1);
    
    
    /// Calculate mp_e, size of in_dege, and copy encode2 to decode2
    for(mp_e=cnt2=0;p11<p12;)
    {
        i=*(indeg+*p11);
        if(i){mp_e+=i;cnt2++;}
        *p13++=*p11++;
    }
    product->m_e=mp_e;
    product->in_edge=in_edge=(EDGE*)malloc((mp_e+cnt2)*sizeof(EDGE));
    printf("mp_e = %lld\n",mp_e);
    printf("# of product states with indeg > 0: %lld\n",cnt2);
    
    /// Initialize encode2 to -1
    for(p11=encode2,p12=p11+n;p11<p12;)*p11++=-1;
    
    /// Set up mapping for encode2, set up pointers for in_edge. We store the in-going edges backwards.
    for(p11=decode2,i=0,cnt2=-1;i<np;p11++,p3++,i++)
    {
        *(encode2+*p11)=i;
        j=*(indeg+*p11);
        if(j){cnt2+=j+1;p3->in_edge=in_edge+cnt2;p3->in_edge->next=0;}
    }
    free(indeg);
    end_t=clock();
    printf("Time used to do compression mapping2 and set up data = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: take product 2-nd time, store in the form of pres
    start_t=clock();
    printf("Take product 2-nd time, store in the form of pres:\n");
    
    /**
     Because of our way to take product to initialize outdeg and to reset
     outdeg, we don't need to initialize outdeg or reset_outdeg to be 0.
     */
    p10=product->outdeg=outdeg=(int*)malloc(np*action*sizeof(int));
    p7=product->acc=acc_p=(BOOL*)malloc(np*sizeof(BOOL));
    p16=product->queue=queue=(long long*)malloc(np*sizeof(long long));
    
    // MARK: take product 2-nd time
    /**
     Note that we have a nice property that the in-going edges of each state in the product space are sorted in the order xp1<, a< after compression mapping 2. This means the pointers are traversing linearly, even when deleting edges in buchi. I store the in-going edges backwards to save some operations. The order is therefore xp1>, a>. This doesn't affect efficency since it's still linear.
     */
    cnt1=0;p11=decode2;p3=nts_product;
    for(;cnt1<np;cnt1++,p11++,p3++,p10+=action)
    {
        xp1=*p11;x1=xp1/n1;q1=xp1%n1;*p7++=*(acc_dba+q1);
        p1=nts_post+x1;
        q2=*(q_prime+p1->label*n1+q1);
        i=p1->num_a;p3->outdeg_reset=p6;
        p8=a_index+p1->pos;p13=a_pos+p1->pos;
        for(;i--;p8++,p13++)
        {
            for(p14=edge+*p13,p15=edge+*(p13+1);p14<p15;p14++)
            {
                p9=q_prime+(nts_post+*p14)->label*n1+q2;
                if(*p9==n1)break;
            }
            if(p14==p15)
            {
                p3->num_a++;p6++->next=1;p6->a=*p8;p14=edge+*p13;
                *(p10+*p8)=p6->outdeg=(int)(p15-p14);
                for(;p14<p15;p14++)
                {
                    xp2=*p14*n1+q2;
                    p4=nts_product+*(encode2+xp2);
                    p5=p4->in_edge--;
                    p5->x=cnt1;p5->a=*p8;
                    p4->in_edge->next=1;
                }
            }
        }
        if(p3->num_a)p6++->next=0;
        else
        {
            p3->outdeg_reset=0;
            p3->iterate--;
            if(p3->in_edge)*p16++=cnt1;else wp--;
        }
    }
    free(encode2);head->encode2=0;
    free(acc_dba);dba->acc=0;
    end_t=clock();
    printf("Time used to take product 2-nd time, "
           "store in the form of pres = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: do safety in nts_product
    start_t=clock();
    printf("Do safety to guarantee no terminal states:\n");
    for(p15=queue;p15<p16;p15++)
        for(p5=(nts_product+*p15)->in_edge;p5->next;)
        {
            p5+=p5->next;p3=nts_product+p5->x;
            if(p3->num_a)
            {
                p10=outdeg+p5->x*action+p5->a;
                if(*p10>0)
                {
                    *p10=-1;p3->num_a--;
                    if(!p3->num_a)
                    {
                        p3->iterate--;
                        if(p3->in_edge)*p16++=p5->x;else wp--;
                    }
                }
            }
        }
    product->wp=wp-=p16-queue;
    end_t=clock();
    printf("Time used to do safety = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    if(!wp)
    {
        printf("wp is empty.\n");
        return;
    }
    printf("Total time used to take product from post: = %.4f\n",
           (double)(end_t-start_t0)/CLOCKS_PER_SEC);
    printf("********\n\n");
}

void buchi_and_controller(HEAD *head)
{
    NTS_PRODUCT *product;
    NODE_PRODUCT *nts_product,*p1,*p2;
    int action,*outdeg,*p3,*p4;
    EDGE *in_edge,*p5,*p6;
    OUTDEG_RESET *outdeg_reset,*p7,*p8;
    BOOL *acc_p,*p9;
    long long *queue,np,w0,w0_1,wp,cnt1,cnt2,i;
    long long *p10,*p11,*p12;
    SCOPE *scope,*target,*p13,*p14;
    long long *decode1,*decode2,*encode3;
    CTRLR *ctrlr;
    NODE_POST *nts_ctrlr1,*nts_ctrlr2,*p15,*p16;
    CTRL *ctrl,*p17,*p18;
    long long n0,n0_1,n0_2,xp;
    int n1;
    FILE *fout;
    clock_t start_t0,start_t,end_t;
    
    if(!head->nts_product.wp)return;
    action=head->action;
    product=&head->nts_product;nts_product=product->graph;
    outdeg=product->outdeg;in_edge=product->in_edge;
    outdeg_reset=product->outdeg_reset;acc_p=product->acc;
    np=product->n;w0_1=w0=product->w0;wp=product->wp;
    p11=p12=queue=product->queue;
    printf("buchi and controller:\n");
    
    // MARK: solve buchi game in product space
    start_t0=start_t=clock();
    printf("Solve buchi game in product space:\n");
    
    // MARK: - set up scope and target for buchi
    p13=product->scope=scope=(SCOPE*)malloc((wp+1)*sizeof(SCOPE));
    for(i=0,p1=nts_product,p9=acc_p;i<np;i++,p1++,p9++)
        if(!p1->iterate)    // iterate can be -1 or 0
        {
            p1->num_a=0;p13++->next=1;p13->x=i;
            if(*p9)*p12++=i;
        }
        else if(i<w0)
        {
            w0_1--;
            if(!w0_1)break;
        }
    p13->next=0;
    cnt1=p12-queue; // # of current acc
    if(!w0_1||!cnt1)
    {
        printf("The winning set is empty, no controller exist.\n");
        return;
    }
    p13=product->target=target=(SCOPE*)malloc((cnt1+1)*sizeof(SCOPE));
    for(;p11<p12;)
    {p13++->next=1;p13->x=*p11++;}
    p13->next=0;
    
    // MARK: compute buchi
    for(i=0;;)
    {
        // MARK: - set up base for reachability
        for(p11=p12=queue,p13=p14=target,cnt2=cnt1,wp=cnt1=0;p14->next;)
        {
            p14+=p14->next;p1=nts_product+p14->x;
            if(p1->iterate==i)
            {
                p13=p14;
                if(p1->in_edge)*p12++=p14->x;else wp++;
            }
            else p13->next+=p14->next;
        }
        p13->next=0;
        // MARK: BFS for reachability
        for(i++;p11<p12;p11++)
        {
            for(p2=nts_product+*p11,p5=p6=p2->in_edge;p6->next;)
            {
                p6+=p6->next;p1=nts_product+p6->x;
                if(p1->iterate>i-2)  // iterate can be i-2, i-1, i
                {
                    p3=outdeg+p6->x*action+p6->a;
                    if(*p3>0)  // Can be either -1 or positive
                    {
                        p5=p6;(*p3)--;
                        if(!*p3)
                        {
                            p1->num_a++;
                            if(p1->iterate<i)  // iterate == i-1
                            {
                                p1->iterate++;
                                if(*(acc_p+p6->x))cnt1++;
                                else if(p1->in_edge)*p12++=p6->x;
                                else wp++;
                            }
                        }
                    }else p5->next+=p6->next;   // delete edge is linear O(m)
                }else p5->next+=p6->next;
            }
            if(p5==p2->in_edge)p2->in_edge=0;else p5->next=0;
        }
        wp+=p12-queue;
        // MARK: end buchi if all the acc can be reached
        if(cnt1==cnt2)break;
        /* start safety */
        // MARK: set up base for safety
        for(p11=p12=queue,p13=p14=target,cnt2-=cnt1;cnt2;)
        {
            p14+=p14->next;p1=nts_product+p14->x;
            if(p1->iterate!=i)
            {
                p13->next+=p14->next;cnt2--;
                if(p1->in_edge)*p12++=p14->x;else wp--;
            }
            else p13=p14;
        }
        if(!p14->next)p13->next=0;
        // MARK: BFS for safety
        for(;p11<p12;p11++)
            for(p6=(nts_product+*p11)->in_edge;p6->next;)
            {
                p6+=p6->next;p1=nts_product+p6->x;
                if(p1->iterate==i) // iterate can be i-2, i-1, i
                {
                    p3=outdeg+p6->x*action+p6->a;
                    if(!*p3)    // Can be -1, 0, or positive
                    {
                        *p3=-1;p1->num_a--;
                        if(!p1->num_a)
                        {
                            p1->iterate--;
                            if(p1->in_edge)*p12++=p6->x;else wp--;
                            if(*(acc_p+p6->x))cnt1--;
                        }
                    }
                }
            }
        wp-=p12-queue;
        // MARK: end buchi if no acc is left
        if(!cnt1)break;
        /* end safety */
        // MARK: reset scope and label
        for(p13=p14=scope;p14->next;)
        {
            p14+=p14->next;p1=nts_product+p14->x;
            if(p1->iterate==i)
            {
                p1->num_a=0;p13=p14;
                for(p3=outdeg+p14->x*action,p7=p8=p1->outdeg_reset;p8->next;)
                {
                    p8+=p8->next;p4=p3+p8->a;
                    if(*p4){*p4=-1;p7->next+=p8->next;}
                    else{*p4=p8->outdeg;p7=p8;}
                }
                p7->next=0;
            }
            else
            {
                p13->next+=p14->next;
                if(p14->x<w0)
                {
                    w0_1--;
                    if(!w0_1)break;
                }
            }
        }
        p13->next=0;
        // MARK: end buchi if no ini is left
        if(!w0_1)break;
    }
    printf("# of loops: %lld\n",i);
    free(target);product->target=0;
    if(cnt1&&w0_1)
    {
        for(p13=p14=scope;p14->next;)
        {
            p14+=p14->next;p1=nts_product+p14->x;
            if(p1->iterate==i)
            {
                p1->num_a=0;p13=p14;
                for(p3=outdeg+p14->x*action,p7=p8=p1->outdeg_reset;p8->next;)
                {
                    p8+=p8->next;p4=p3+p8->a;
                    if(*p4){*p4=-1;p7->next+=p8->next;}
                    else{*p4=p8->outdeg;p7=p8;}
                }
                p7->next=0;
            }
            else
            {
                p13->next+=p14->next;
                if(p14->x<w0)
                {
                    w0_1--;
                    if(!w0_1)break;
                }
            }
        }
        p13->next=0;
    }
    else w0_1=0;
    if(!w0_1)
    {
        printf("The winning set is empty, no controller exist.\n");
        return;
    }
    product->w0_1=w0_1;free(outdeg_reset);product->outdeg_reset=0;
    end_t=clock();
    printf("Time used to solve buchi game in product space = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: - compute controller
    start_t=clock();
    printf("Compute controller:\n");
    decode1=head->decode1;decode2=head->decode2;
    n0=head->nts_pre.n;n0_1=head->nts_post.n;
    n1=head->dba.n;
    ctrlr=&head->ctrlr;ctrlr->w0=w0_1;
    ctrlr->wp_acc=cnt1;ctrlr->wp=wp;
    
    printf("# of states in W_X0: %lld\n",w0_1);
    
    //        printf("W_X0:\n");
    //        for(p13=scope,cnt2=0;cnt2<w0_1;cnt2++)
    //        {
    //            p13+=p13->next;
    //            printf("%lld\t",*(decode1+*(decode2+p13->x)/n1));
    //        }
    //        putchar('\n');
    printf("# of states in WP_acc: %lld\n",cnt1);
    printf("# of states in W_P: %lld\n",wp);
    
    // MARK: set up data to build controller
    /// Set nts_ctrlr1
    p15=nts_ctrlr1=(NODE_POST*)malloc(n0_1*sizeof(NODE_POST));
    for(p16=head->nts_post.graph,cnt2=0;cnt2<n0_1;cnt2++)
    {p15->num_a=0;p15++->label=p16++->label;}
    
    for(p13=scope;p13->next;)
    {
        p13+=p13->next;
        (nts_ctrlr1+*(decode2+p13->x)/n1)->num_a++;
    }
    
    /// Set encode3
    p10=head->encode3=encode3=(long long*)malloc(n0*sizeof(long long));
    for(p11=head->encode1,p12=p11+n0;p11<p12;)
        *p10++=*p11++;
    
    for(n0_2=0,p10=decode1,p15=nts_ctrlr1,p16=p15+n0_1;p15<p16;p10++,p15++)
        *(encode3+*p10)=p15->num_a?n0_2++:-1;
    ctrlr->n=n0_2;
    printf("n0_2: %lld\n",n0_2);
    
    /// Set nts_ctrlr2
    ctrlr->graph=nts_ctrlr2=(NODE_POST*)malloc(n0_2*sizeof(NODE_POST));
    for(cnt2=0,p15=nts_ctrlr1,p16=nts_ctrlr2;cnt2!=wp;p15++)
        if(p15->num_a)
        {
            p16->label=p15->label;p16->pos=p15->pos=cnt2;
            cnt2+=p16->num_a=p15->num_a;p16++;
        }
    ctrlr->ctrl=ctrl=(CTRL*)malloc(wp*sizeof(CTRL));
    
    /* for(p12=queue+cnt1;!*(acc_p+*(p12-1));p12--); */
    for(p12=p10=queue;*(acc_p+*p10);p10++)
      if((nts_product+*p10)->in_edge)*p12++=*p10;
    for(i++,p11=queue;p11<p12;p11++)
    {
        for(p2=nts_product+*p11,p5=p6=p2->in_edge;p6->next;)
        {
            p6+=p6->next;p1=nts_product+p6->x;
            if(p1->iterate>i-2)  // iterate can be i-2, i-1, i
            {
                p3=outdeg+p6->x*action+p6->a;
                if(*p3>0)  // Can be either -1 or positive
                {
                    p5=p6;(*p3)--;
                    if(!*p3)
                    {
                        p1->num_a++;
                        if(p1->iterate<i)  // iterate == i-1
                        {
                            p1->iterate++;
                            xp=*(decode2+p6->x);
                            p15=nts_ctrlr1+xp/n1;
                            p17=ctrl+p15->pos++;
                            p17->q=xp%n1;p17->u=p6->a;
                            if(!*(acc_p+p6->x)&&p1->in_edge)*p12++=p6->x;
                        }
                    }
                }else p5->next+=p6->next;   // delete edge is linear O(m)
            }else p5->next+=p6->next;
        }
        if(p5==p2->in_edge)p2->in_edge=0;else p5->next=0;
    }
    free(nts_ctrlr1);
    free(outdeg);product->outdeg=0;free(in_edge);product->in_edge=0;
    free(acc_p);product->acc=0;free(queue);product->queue=0;
    free(nts_product);product->graph=0;
    
    for(p15=nts_ctrlr2,p16=p15+n0_2;p15<p16;p15++)
        qsort(ctrl+p15->pos,p15->num_a,sizeof(CTRL),cmp2);
    end_t=clock();
    printf("Time used to compute controller = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    
    // MARK: write controller
    start_t=clock();
    printf("Write controller:\n");
    char buffer[40];
    sprintf(buffer,"controller%d_%d_%d.txt",
            head->nts_num,head->labels_num,head->dba_num);
    printf("Write controller into %s:\n",buffer);
    fout=fopen(buffer,"wb");
    
    /// W_X0 size: w0_1 * (long long == 2 * int)
    printf("W_X0\tsize: %lld * (long long == 2 * int)\n",w0_1);
    fprintf(fout,"%lld\n",w0_1);
    for(p13=scope,cnt2=0;cnt2<w0_1;cnt2++)
    {
        p13+=p13->next;
        fprintf(fout,"%lld\n",*(decode1+*(decode2+p13->x)/n1));
    }
    free(scope);product->scope=0;
    free(decode2);head->decode2=0;
    
    /// encode3 size: n0 * (long long == 2 * int)
    printf("encode3\tsize: %lld * (long long == 2 * int)\n",n0);
    fprintf(fout,"%lld\n",n0);
    for(p11=encode3,p12=p11+n0;p11<p12;)
        fprintf(fout,"%lld\n",*p11++);
    free(encode3);head->encode3=0;
    
    /// ctrlr size: n0_2 * (NODE_POST == 4 * int)
    printf("ctrlr\tsize: %lld * (NODE_POST == 4 * int)\n",n0_2);
    fprintf(fout,"%lld\n",n0_2);
    for(p15=nts_ctrlr2,p16=p15+n0_2;p15<p16;p15++)
        fprintf(fout,"%d %d %lld\n",p15->num_a,p15->label,p15->pos);
    free(nts_ctrlr2);ctrlr->graph=0;
    
    /// ctrl size: wp * (CTRL == 2 * int)
    printf("ctrl\tsize: %lld * (CTRL == 2 * int)\n",wp);
    fprintf(fout,"%lld\n",wp);
    for(p17=ctrl,p18=p17+wp;p17<p18;p17++)
        fprintf(fout,"%d %d\n",p17->q,p17->u);
    free(ctrl);ctrlr->ctrl=0;
    
    /// q_prime size: 2^k*n1 * (int)
    cnt2=head->dba.q_prime_size;
    p3=head->dba.q_prime;
    printf("q_prime\tsize: %lld * (int)\n",cnt2);
    fprintf(fout,"%lld\n",cnt2);
    for(p4=p3+cnt2;p3<p4;)
        fprintf(fout,"%d\n",*p3++);
    free(head->dba.q_prime);head->dba.q_prime=0;
    fclose(fout);
    end_t=clock();
    printf("Time used to write controller = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    printf("Total time used to solve buchi game in product space "
           "and controller = %.4f\n",
           (double)(end_t-start_t0)/CLOCKS_PER_SEC);
    printf("********\n\n");
}

void free_memory(HEAD *head,int flag)
{
    NTS_PRE *pre;
    NTS_POST *post;
    DBA *dba;
    NTS_PRODUCT *product;
    CTRLR *ctrlr;
    clock_t start_t,end_t;
    
    pre=&head->nts_pre;
    post=&head->nts_post;
    dba=&head->dba;
    product=&head->nts_product;
    ctrlr=&head->ctrlr;
    
    // MARK: free memory
    start_t=clock();
    printf("Free memory:\n");
    /// Free head
    if(flag)
    {
        if(head->encode1){free(head->encode1);head->encode1=0;}
        if(head->decode1){free(head->decode1);head->decode1=0;}
    }
    if(head->encode2){free(head->encode2);head->encode2=0;}
    if(head->decode2){free(head->decode2);head->decode2=0;}
    if(head->encode3){free(head->encode3);head->encode3=0;}
    
    /// Free nts_pre
    if(flag)    //  we only need to free it for the last time
    {
        if(pre->outdeg){free(pre->outdeg);pre->outdeg=0;}
        if(pre->in_edge){free(pre->in_edge);pre->in_edge=0;}
        if(pre->graph){free(pre->graph);pre->graph=0;}
    }
    
    /// Free nts_post
    if(flag)
    {
        if(post->a_index){free(post->a_index);post->a_index=0;}
        if(post->a_pos){free(post->a_pos);post->a_pos=0;}
        if(post->edge){free(post->edge);post->edge=0;}
        if(post->graph){free(post->graph);post->graph=0;}
    }
    
    /// Free dba
    if(flag)
    {
        if(dba->labels){free(dba->labels);dba->labels=0;}
    }
    if(dba->q_prime){free(dba->q_prime);dba->q_prime=0;}
    if(dba->acc){free(dba->acc);dba->acc=0;}
    
    /// Free nts_product
    if(product->outdeg){free(product->outdeg);product->outdeg=0;}
    if(product->in_edge){free(product->in_edge);product->in_edge=0;}
    if(product->outdeg_reset)
    {free(product->outdeg_reset);product->outdeg_reset=0;}
    if(product->acc){free(product->acc);product->acc=0;}
    if(product->scope){free(product->scope);product->scope=0;}
    if(product->target){free(product->target);product->target=0;}
    if(product->queue){free(product->queue);product->queue=0;}
    if(product->graph){free(product->graph);product->graph=0;}
    product->wp=0;
    
    /// Free ctrlr
    if(ctrlr->ctrl){free(ctrlr->ctrl);ctrlr->ctrl=0;}
    if(ctrlr->graph){free(ctrlr->graph);ctrlr->graph=0;}
    
    end_t=clock();
    printf("Time used to  = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);
    printf("********\n\n");
}
//    // MARK:
//    start_t=clock();
//    printf(":\n");
//
//
//    end_t=clock();
//    printf("Time used to  = %.4f\n",
//           (double)(end_t-start_t)/CLOCKS_PER_SEC);

//int main(int argc, const char * argv[]) {
//    HEAD head,*p=&head;
//    int dba_num[10]={3,4},i,cnt=1;
//    clock_t start_t,end_t;
//    
//    p->nts_num=1,p->labels_num=1;
//    // Remember: set dba_num in dba_num[10] and # of specifications in cnt above;
//    
//    initialization(p);
//    input_pre(p);
//    pre2post(p);
//    for(i=0;i<cnt;)
//    {
//        start_t=clock();
//        p->dba_num=dba_num[i++];
//        post2product(p);
//        buchi_and_controller(p);
//        if(i!=cnt)
//            free_memory(p,0);
//        else
//            free_memory(p,1);
//        end_t=clock();
//        printf("Time used to solve dba%d = %.4f\n\n",
//               dba_num[i-1],(double)(end_t-start_t)/CLOCKS_PER_SEC);
//    }
//    
//    printf("Total time used = %.4f\n",(double)clock()/CLOCKS_PER_SEC);
//    
//    //printf("size of EDGE: %d\n",sizeof(EDGE));
//    
//    
//    return 0;
//}
//

















#endif /* ProductAutomaton_h */

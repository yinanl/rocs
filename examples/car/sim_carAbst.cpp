/**
 *  sim_carAbst.cpp
 *
 *  To simulate the vehicle example using abstraction-based control.
 *  
 *  Authors: Yinan Li, Zhibing Sun, Jun Liu
 *  Created: May 27, 2020
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include <string>
#include <math.h>

#include "src/grid.h"
#include "src/definitions.h"
#include <boost/numeric/odeint.hpp>

/* define dynamics */

struct car_dynamics {
	rocs::Rn u;
	car_dynamics (const rocs::Rn param): u (param) {}

	void operator() (rocs::Rn &x, rocs::Rn &dxdt, double t) const 
	{
	      double alpha = std::atan(std::tan(u[1])/2.0);
	      dxdt[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
	      dxdt[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
	      dxdt[2] = u[0]*std::tan(u[1]);
	}
};

/* ??? */
typedef char BOOL;

typedef struct NODE_POST
{int num_a,label;long long pos;}NODE_POST;

typedef struct CTRL
{int q,u;}CTRL;

int cmp1(const void *a,const void *b)
{return (int)(*(long long*)a-*(long long*)b);}

int cmp2(const void *a,const void *b)
{
    return ((CTRL*)a)->q-((CTRL*)b)->q;
    //    CTRL *p1=(CTRL*)a,*p2=(CTRL*)b;
    //    return p1->q-p2->q;
}


int main(int argc, char *argv[])
{
    
    /* set the state space */
    const int x_dim = 3;
    const double theta = 3.5;
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};

    /* set the control values */
    const int u_dim = 2;
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};

    /* discretization precision */
    double mu[] = {0.3, 0.3};
    std::string scale_eta ("1");
    if (argc>1) {
	    scale_eta = argv[1];
    }
    double set_eta = 0.2/pow(2,std::stod(scale_eta)-1); 
    const double eta[] = {set_eta, set_eta, set_eta}; /* set precision */

    const double tau = 0.3;

    /* ode solver */
    const double dt = 0.001; //integration step size for odeint
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;

   /* generate grid */
    rocs::grid x_grid(x_dim,eta,xlb,xub);
    x_grid.gridding();
    rocs::grid u_grid(u_dim,mu,ulb,uub);

    std::cout << "Number of discrete states: " << x_grid._nv << "\n";
    std::cout << "Number of discrete inputs: " << u_grid._nv << "\n";

   /* simulation */

    int nts_num=1,labels_num=1,dba_num=3;

    /* choose nts, labels, dba files by argument */
    if (argc > 1){
	    nts_num = atoi(argv[1]);
           if (argc >2)  dba_num = atoi(argv[2]);
    }

    int max_num_achieve_acc=5,max_num_iteration=3000000;
    long long *w_x0,w0,*encode3,n0,n0_2,wp,q_prime_size;
    long long *p1,*p2,x_index;
    NODE_POST *nts_ctrlr,*p5,*p6;
    CTRL *ctrl,*p7,*p8;
    int n1,m1,ini,*q_prime,*p9,*p10,q,i,j;
    BOOL *acc,*p11,*p12;
    char buffer[40];
    clock_t start_t0,start_t,end_t;
    FILE *fin,*fout;

    // MARK: input dba
    start_t0=start_t=clock();
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

    sprintf(buffer,"dba%d.txt",dba_num);
    std::cout << "Read dba from " << buffer << ": ";
    fin=fopen(buffer,"rb");
    fscanf(fin,"%d%d%d",&n1,&m1,&ini);
    p11=acc=(BOOL*)malloc(n1*sizeof(BOOL));
    for(i=0;i<n1;i++)
    {
        fscanf(fin,"%d",&j);
        *p11++=j;
    }
    fclose(fin);
    printf("Done.\n");

    // MARK: input controller
    printf("input controller:\n");
    sprintf(buffer,"controller%d_%d_%d.txt",nts_num,labels_num,dba_num);
    std::cout << "Read controller from " << buffer << ":\n";
    fin=fopen(buffer,"rb");

    /// W_X0 size: w0_1 * (long long == 2 * int)
    fscanf(fin,"%lld",&w0);
    printf("W_X0\tsize: %lld * (long long == 2 * int)\n",w0);
    p1=w_x0=(long long*)malloc(w0*sizeof(long long));
    for(p2=p1+w0;p1<p2;)
        fscanf(fin,"%lld",p1++);

    /// encode3 size: n0 * (long long == 2 * int)
    fscanf(fin,"%lld",&n0);
    printf("encode3\tsize: %lld * (long long == 2 * int)\n",n0);
    p1=encode3=(long long*)malloc(n0*sizeof(long long));
    for(p2=p1+n0;p1<p2;)
        fscanf(fin,"%lld",p1++);

    /// nts_ctrlr size: n0_2 * (NODE_POST == 4 * int)
    fscanf(fin,"%lld",&n0_2);
    printf("ctrlr\tsize: %lld * (NODE_POST == 4 * int)\n",n0_2);
    p5=nts_ctrlr=(NODE_POST*)malloc(n0_2*sizeof(NODE_POST));
    for(p6=p5+n0_2;p5<p6;p5++)
        fscanf(fin,"%d %d %lld",&p5->num_a,&p5->label,&p5->pos);

    /// ctrl size: wp * (CTRL == 2 * int)
    fscanf(fin,"%lld",&wp);
    printf("ctrl\tsize: %lld * (CTRL == 2 * int)\n",wp);
    p7=ctrl=(CTRL*)malloc(wp*sizeof(CTRL));
    for(p8=p7+wp;p7<p8;p7++)
        fscanf(fin,"%d %d",&p7->q,&p7->u);

    /// q_prime size: 2^k*n1 * (int)
    fscanf(fin,"%lld",&q_prime_size);
    printf("q_prime\tsize: %lld * (int)\n",q_prime_size);
    p9=q_prime=(int*)malloc(q_prime_size*sizeof(int));
    for(p10=p9+q_prime_size;p9<p10;)
        fscanf(fin,"%d",p9++);
    fclose(fin);
    end_t=clock();
    printf("Time used to input dba and controller = %.4f\n",
           (double)(end_t-start_t)/CLOCKS_PER_SEC);

    rocs::Rn x{5, 5, 0};
    rocs::Rn u{0, 0};	
    printf("dba%d info:\n",dba_num);
    printf("n1 = %d, m1 = %d, ini = %d\n",n1,m1,ini);
    printf("acc mark: ");
    for(p11=acc,p12=p11+n1;p11<p12;)
        printf("%d ",*p11++);
    putchar('\n');
    std::cout << "Input initial state x as x[0] x[1] x[2] : " << std::endl;

    for(;scanf("%lf%lf%lf",&x[0],&x[1],&x[2])!=EOF;)
    {   
        x_index=x_grid.val_to_id(x);
	std::cout << x_index << "\n";
        p1=(long long*)bsearch(&x_index,w_x0,w0,sizeof(long long),cmp1);
        if(!p1)
        {   
            std::cout << x[0] <<  " "  << x[1] << " " << x[2];
            std::cout << " is not in the winning set, change a new initial state." << std::endl;
            continue;
        }
        fout=fopen("controller_path.txt","wb");
        for(q=ini,i=j=0;i<max_num_achieve_acc&&j<max_num_iteration;)
        {   
            p5=nts_ctrlr+*(encode3+x_index);
            p7=(CTRL*)bsearch(&q,ctrl+p5->pos,p5->num_a,sizeof(CTRL),cmp2);
            if(!p7) {std::cout << "error in ctrl" << std::endl;break;}
            u_grid.id_to_val(u,p7->u);
            
            std::cout << "# " << j++ << ":\n";  
            std::cout << "x: " << x[0] <<  " "  << x[1] << " " << x[2] << "\n";
            
            fprintf(fout,"%lf %lf\n",x[0],x[1]);
            
            std::cout << "u: " << u[0] <<  " "  << u[1] << std::endl;
            std::cout << "current dba state: " << q << "\n";
	    /* integrate the dynamics */
	    boost::numeric::odeint::integrate_const(rk45, car_dynamics(u), x, 0.0, tau, dt);
	    std::cout << "current x values: " << x[0] << " " << x[1] << " "  << x[2] << "\n";  
            x_index=x_grid.val_to_id(x);
            q=*(q_prime+p5->label*n1+q);
            if(*(acc+q))std::cout << "******** achieve acc \'" << q << "\' " << ++i << " time ********" << std::endl;
        }
        fclose(fout);
        printf("dba%d info:\n",dba_num);
        printf("n1 = %d, m1 = %d, ini = %d\n",n1,m1,ini);
        printf("acc mark: ");
        for(p11=acc,p12=p11+n1;p11<p12;)
            printf("%d ",*p11++);
        putchar('\n');
        std::cout << "Input initial state x as x[0] x[1] x[2]: " << std::endl;
    }
   
    free(acc),free(w_x0),free(encode3);
    free(nts_ctrlr),free(ctrl),free(q_prime);

    /* plot simulated trajectory 
     * to do!
     * */


//    /* test grid and odeint */
//
//    std::vector <double> x{5.03, 5, 0};
//    std::vector <double> u{0, 0};
//
//    rocs::Rn y{1, 1, 0};
//
//    std::cout<< x_grid.val_to_id(x) << "\n";
//    x_grid.id_to_val(y,45517);
//    std::cout<< y[0] << " " << y[1] << " " << y[2] << "\n";
//    //std::cout<< u_grid.val_to_id(u) << "\n";
//    //u_grid.id_to_val(u,23);
//    //std::cout<< u[0] << " " << u[1] << " " << "\n";
//
//    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;
//    y = {5,5,0};
//    y = {5,5,0.0};
//    u = {2.0,2.0};
//    boost::numeric::odeint::integrate_const(rk45, car_dynamics(u), y, 0.0, tau, dt);
//    std::cout<< y[0] << " " << y[1] << " " << y[2] << "\n";


    return 0;
}

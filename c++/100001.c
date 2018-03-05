/* 100001.c */

#include <R.h>

static double parms[8];
#define pop parms[0]
#define hatch parms[1]
#define mort parms[2]
#define Pr_M parms[3]
#define beta parms[4]
#define gamma parms[5]
#define sigma parms[6]
#define h parms[7]

static double forcs[3];
#define B forcs[0]
#define M forcs[1]
#define T forcs[2]


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
    int N=8;
    odeparms(&N, parms);
}

void initforc(void (* odeforcs)(int *, double *))
{
    int N=3;
    odeforcs(&N, forcs);
}


/* function */
void modSIRS (int *neq, double *t, double *y,
              double *ydot,
              double *yout, int *ip)
{
    
    ydot[0] = -beta*(y[1] + y[4])*y[0]                 + sigma*y[2] - h*y[0]  - mort*y[0] + B*0.5*pop*hatch;
    ydot[1] =  beta*(y[1] + y[4])*y[0]   - gamma*y[1]               + h*y[0]  - mort*y[1];
    ydot[2] =                            + gamma*y[1]  - sigma*y[2]           - mort*y[2];
    
       
    ydot[3] = -beta*(y[1] + y[4])*y[3]                 + sigma*y[5] - h*y[3]  - mort*y[3]                           + M*(Pr_M*pop*(y[0]/(y[0]+y[1]+y[2]))) + T*y[5];
    ydot[4] =  beta*(y[1] + y[4])*y[3]   - gamma*y[4]               + h*y[3]  - mort*y[4] 				   + M*(Pr_M*pop*(y[1]/(y[0]+y[1]+y[2])));
    ydot[5] =				   + gamma*y[4]  - sigma*y[5] 		  - mort*y[5] 	                         + M*(Pr_M*pop*(y[2]/(y[0]+y[1]+y[2]))) - T*y[5];
    
    yout[0] = (y[0] + y[1] + y[2] + y[3] + y[4] + y[5]);
    
}

/* event function */
void event(int *n, double *t, double *y)
{
	y[0] = y[0];
	y[1] = y[1];
	y[2] = y[2];
		  
    y[3] = 0;
    y[4] = 0;
    y[5] = 0;

    
}

/* END file 100001.c */
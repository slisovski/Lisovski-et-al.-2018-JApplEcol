/* 111100.c */

#include <R.h>

static double parms[10];
#define pop parms[0]
#define hatch parms[1]
#define mort parms[2]
#define mort_j parms[3]
#define Pr_M parms[4]
#define beta parms[5]
#define gamma parms[6]
#define sigma parms[7]
#define h parms[8]
#define b parms[9]

static double forcs[2];
#define B forcs[0]
#define M forcs[1]


/* initializers */
void initmod(void (* odeparms)(int *, double *))
{
    int N=9;
    odeparms(&N, parms);
}

void initforc(void (* odeforcs)(int *, double *))
{
    int N=2;
    odeforcs(&N, forcs);
}


/* function */
void modSIRS (int *neq, double *t, double *y,
              double *ydot,
              double *yout, int *ip)
{
    
    ydot[0] = -beta*(y[1] + y[4] + y[7])*y[0]                 + sigma*y[2] - h*y[0]  - mort*y[0];
    ydot[1] =  beta*(y[1] + y[4] + y[7])*y[0]   - gamma*y[1]               + h*y[0]  - mort*y[1];
    ydot[2] =                                   + gamma*y[1]  - sigma*y[2]           - mort*y[2];
    
    ydot[3] = -beta*(y[1] + y[4] + y[7])*y[3]                 + sigma*y[5] - h*y[3]  - mort_j*y[3] + B*0.5*pop*hatch;
    ydot[4] =  beta*(y[1] + y[4] + y[7])*y[3]   - gamma*y[4]               + h*y[3]  - mort_j*y[4];
    ydot[5] =                                   + gamma*y[4]  - sigma*y[5]           - mort_j*y[5];    
       
    ydot[6] = -beta*(y[1] + y[4] + y[7])*y[6]                 + sigma*y[8] - h*y[6]  - mort*y[6]                           + M*(Pr_M*pop*((y[0]/(y[0]+y[1]+y[2]))+((y[2]/(y[0]+y[1]+y[2]))*b)));
    ydot[7] =  beta*(y[1] + y[4] + y[7])*y[6]   - gamma*y[7]               + h*y[6]  - mort*y[7]                           + M*(Pr_M*pop*(y[1]/(y[0]+y[1]+y[2])));
    ydot[8] =							        + gamma*y[7]  - sigma*y[8] 		     - mort*y[8] + M*(Pr_M*pop*((y[2]/(y[0]+y[1]+y[2]))-((y[2]/(y[0]+y[1]+y[2]))*b)));
    
}

/* event function */
void event(int *n, double *t, double *y)
{
	y[0] = y[0]+y[3];
	y[1] = y[1]+y[4];
	y[2] = y[2]+y[5];
		  
    y[3] = 0;
    y[4] = 0;
    y[5] = 0;
    
    y[6] = 0;
    y[7] = 0;
    y[8] = 0;
    
}

/* END file 111100.c */
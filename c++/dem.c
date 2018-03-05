/* dem.c */

#include <R.h>

static double parms[4];
#define mu_a parms[0]
#define mu_j parms[1]
#define pop parms[2]
#define hatch parms[3]

static double forcs[1];
#define B forcs[0]


/* initializers */
void initmod (void (* odeparms)(int *, double *))
{
    int N=4;
    odeparms(&N, parms);
}

void initforc (void (* odeforcs)(int *, double *))
{
    int N=1;
    odeforcs(&N, forcs);
}


/* function */
void popM (int *neq, double *t, double *y,
              double *ydot,
              double *yout, int *ip)
{
    
    ydot[0] = - mu_a*y[0];
    ydot[1] = - mu_j*y[1] + B*(0.5*pop)*hatch;
    
    yout[1] = y[0] + y[1];
        
}

/* event function */
void event (int *n, double *t, double *y)
{
	y[0] = y[0] + y[1];
	y[1] = 0;
    
}

/* END file dem.c */
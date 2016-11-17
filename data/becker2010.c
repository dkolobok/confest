#include <R.h>
static double parms[8];
#define kon parms[0]
#define koff parms[1]
#define kt parms[2]
#define init_EpoR parms[3]
#define kex parms[4]
#define ke parms[5]
#define kdi parms[6]
#define kde parms[7]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N = 8;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");
ydot[0] = koff * y[2] - kon * y[1] * y[0] + kt * init_EpoR - kt * y[0] + kex * y[3];
ydot[1] = koff * y[2] - kon * y[1] * y[0] + kex * y[3];
ydot[2] = kon * y[1] * y[0] - koff * y[2] - ke * y[2];
ydot[3] = ke * y[2] - kex * y[3] - kdi * y[3] - kde * y[3];
ydot[4] = kdi * y[3];
ydot[5] = kde * y[3];
  
//yout[0] = y[0];
//yout[1] = y[1];
//yout[2] = y[2];
//yout[3] = y[3];
//yout[4] = y[4];
//yout[5] = y[5];
}

/* The Jacobian matrix  */

void jac(int *neq, double *t, double *y, int *ml, int *mu,
double *pd, int *nrowpd, double *yout, int *ip)
{
pd[0] = -kon * y[1] - kt;
pd[1] = -kon * y[1];
pd[2] = kon * y[1];
pd[3] = 0.0;
pd[4] = 0.0;
pd[5] = 0.0;

pd[(*nrowpd)] = -kon * y[0];
pd[(*nrowpd) + 1] = -kon * y[0];
pd[(*nrowpd) + 2] = kon * y[0];
pd[(*nrowpd) + 3] = 0.0;
pd[(*nrowpd) + 4] = 0.0;
pd[(*nrowpd) + 5] = 0.0;

pd[2 * (*nrowpd)] = koff;
pd[2 * (*nrowpd) + 1] = koff;
pd[2 * (*nrowpd) + 2] = -ke;
pd[2 * (*nrowpd) + 3] = ke;
pd[2 * (*nrowpd) + 4] = 0.0;
pd[2 * (*nrowpd) + 5] = 0.0;

pd[3 * (*nrowpd)] = kex;
pd[3 * (*nrowpd) + 1] = kex;
pd[3 * (*nrowpd) + 2] = 0.0;
pd[3 * (*nrowpd) + 3] = -kex - kdi - kde;
pd[3 * (*nrowpd) + 4] = kdi;
pd[3 * (*nrowpd) + 5] = kde;

pd[4 * (*nrowpd)] = 0.0;
pd[4 * (*nrowpd) + 1] = 0.0;
pd[4 * (*nrowpd) + 2] = 0.0;
pd[4 * (*nrowpd) + 3] = 0.0;
pd[4 * (*nrowpd) + 4] = 0.0;
pd[4 * (*nrowpd) + 5] = 0.0;

pd[5 * (*nrowpd)] = 0.0;
pd[5 * (*nrowpd) + 1] = 0.0;
pd[5 * (*nrowpd) + 2] = 0.0;
pd[5 * (*nrowpd) + 3] = 0.0;
pd[5 * (*nrowpd) + 4] = 0.0;
pd[5 * (*nrowpd) + 5] = 0.0;
}

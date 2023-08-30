#include <R.h>

static double parms[500];

#define n parms[0]
#define alpha parms[1]
#define beta parms[2]
#define gamma parms[3]

void initmod(void (* odeparms)(int *, double *)) {
    int N = 500;
    odeparms(&N, parms);
}

void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {

    int E_index = n;
    int I_index = 2 * n;
    int R_index = 3 * n;

    for (int i = 0; i < n; i++) {
        int m_index = 4 + i;
        double SI = 0;
        for (int j = 0; j < n; j++) {
            double I = y[I_index + j];
            SI += parms[m_index + j * (int) n] * I;
        }

        // dS
        double S = y[i];
        ydot[i] = -S*beta*SI;

        // dE
        double E = y[E_index + i];
        ydot[E_index + i] = S*beta*SI - alpha * E;

        // dI
        double I = y[I_index + i];
        ydot[I_index + i] = alpha * E - gamma * I;

        // dR
        ydot[R_index + i] = gamma * I;
    }

}

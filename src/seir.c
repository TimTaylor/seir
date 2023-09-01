#include <R.h>

static double parms[500];

#define n parms[0]
#define alpha parms[1]
#define beta parms[2]
#define gamma parms[3]

// - after gamma we then have n vaccination start times; parms[4:(4 + n - 1)].
// - we then have n vaccination end times; parms[(4+n):(4 + 2*n - 1)]
// - we then have n nu values; parms[(4 + 2*n):(4 + 3*n - 1)]
// - matrix values then occur; parms[(4 + 3*n):(4 + 3*n - 1 + n*n)]
// - the number of interventions, nint; parms[4 + 3*n + n*n]
// - the intervention start times, parms[(5 + 3*n + n*n):(4 + 3*n + n*n + nint)]
// - the intervention end times, parms[(5 + 3*n + n*n + nint):(4 + 3*n + n*n + 2*nint)]
// - the unpacked intervention matrix[(5 + 3*n + n*n + 2*nint): (4 + 3*n + n*n + 2*nint + n*nint)]

void initmod(void (* odeparms)(int *, double *)) {
    int N = 500;
    odeparms(&N, parms);
}

void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {

    int E_index = n;
    int I_index = 2 * n;
    int R_index = 3 * n;

    // We need to work out the contact reduction (if any)
    int nn = n;
    int n_int = parms[4 + 3*nn + (nn*nn)];
    int which_int = -1;
    int int_index = 0;
    if (n_int > 0) {
        for (int i = 0; i < n_int; i++) {
            if (parms[5 + 3*nn + (nn*nn) + i] <= *t && *t < parms[5 + 3*nn + (nn*nn) + n_int + i]) {
                which_int = i;
                int_index = 5 + 3*nn + (nn*nn) + (2 + i)*n_int;
                break;
            }
        }
    }

    for (int i = 0; i < n; i++) {

        // account for vaccination
        int v_start_index = 4 + i;
        int v_end_index = v_start_index + n;
        int nu_index = v_end_index + n;
        double nu = 0;
        if ((parms[v_start_index] <= *t) && (*t < parms[v_end_index])) {
            nu = parms[nu_index];
        }

        // workout contact reduction
        double reduction = 0;
        if (which_int >= 0) {
            reduction = parms[int_index + i];
        }


        // account for SI contacts
        int m_index = 4 + 3*n + i;
        double SI = 0;
        for (int j = 0; j < n; j++) {
            double I = y[I_index + j];
            SI += parms[m_index + j * (int) n] * I * (1 - reduction);
        }

        // dS
        double S = y[i];
        ydot[i] = -S*beta*SI - S*nu;

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

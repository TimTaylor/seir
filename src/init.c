#include <stdlib.h> // for NULL
# include <R.h>
# include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern void initmod(void (* odeparms)(int *, double *));
extern void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int*ip);

static const R_CMethodDef CEntries[] = {
    {"initmod",     (DL_FUNC) &initmod,     1},
    {"derivs",       (DL_FUNC) &derivs,       6},
    {NULL, NULL, 0}
};

void R_init_seir(DllInfo *dll) {
    // register entry points
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);

    // the following two lines protect against accidentially finding entry points
    R_useDynamicSymbols(dll, FALSE);  // disable dynamic searching
    //R_forceSymbols(dll, TRUE);      // entry points as R objects, not as strings
}

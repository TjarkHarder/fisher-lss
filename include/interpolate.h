#ifndef INTERPOLATE_H_INCLUDED
#define INTERPOLATE_H_INCLUDED

/* Splinter */
#include "cinterface.h"


#include "common.h"

#include "misc.h"
#include "dat.h"

#include "shape.h"


/*  ----------------------------------------------------  */
/*  -------------------   Structures   -----------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Parameters for the spline interpolation.

    */

    /* Acceleration struct from GSL (not used with SPLINTER) */
    void *acc;

    /* Interpolation (spline) function */
    void *spline;

    /* Interpolation region */
    double **xBounds;
    double *yBounds;

    /* Must filter out degenerate x data points (SPLINTER does not like them :)) */
    size_t xDim;
    size_t xDimUniq;
    size_t *xIndUniq;

    /* Extrapolation function */
    double (*extrapolate)(void*, void*, void*, size_t*, size_t); // First argument: xValues, second argument: interp_t struct, third argument: additional parameters, fourth argument: indices of xValue array that need to be extrapolated, fifth argument: size of fourth argument

} interp_t;



/*  ----------------------------------------------------  */
/*  ------------------   Initialise   ------------------  */
/*  ----------------------------------------------------  */


int interpolate_ini();
int interpolate_free();



/*  ----------------------------------------------------  */
/*  ----------------   Interpolation   -----------------  */
/*  ----------------------------------------------------  */


interp_t *interpolate_interp_init(
    interp_t *(*interpInit)(dat_t*, double (*)(void*, void*, void*, size_t*, size_t)),
    dat_t *dat,
    double (*extrapolate)(void*, void*, void*, size_t*, size_t));


interp_t *interpolate_interp_init_splinter(
    dat_t *dat,
    double (*extrapolate)(void*, void*, void*, size_t*, size_t));

interp_t *interpolate_interp_init_gsl(
    dat_t *dat,
    double (*extrapolate)(void*, void*, void*, size_t*, size_t));

/*  ----------------------------------------------------  */

double interpolate_interp_eval(
    double (*interpEval)(void*, interp_t*, void*),
    void *values,
    interp_t *interp,
    void *params);


double interpolate_interp_eval_splinter(
    void *values,
    interp_t *interp,
    void *params);


double interpolate_interp_eval_gsl(
    void *values,
    interp_t *interp,
    void *params);

double interpolate_interp_eval_deriv_gsl(
    void *values,
    interp_t *interp,
    void *params);

/*  ----------------------------------------------------  */

interp_t *interpolate_interp_free(
    interp_t *(*interpFree)(interp_t*),
    interp_t *interp);


interp_t *interpolate_interp_free_splinter(
    interp_t *interp);

interp_t *interpolate_interp_free_gsl(
    interp_t *interp);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // INTERPOLATE_H_INCLUDED

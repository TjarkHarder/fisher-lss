#ifndef INTEGRAND_AVERAGE_H_INCLUDED
#define INTEGRAND_AVERAGE_H_INCLUDED

#include "common.h"

#include "misc.h"

#include "fiducials.h"

#include "kernels.h"
#include "interpolate.h"
#include "spectra.h"


/*  ----------------------------------------------------  */
/*  --------------------   Struct   --------------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Bin-Average Struct

    */

    double (*specFunc)(void*, void*);

    void *var;
    void *params;

    double *bounds;

} cintgrnd_avr_t;



/*  ----------------------------------------------------  */
/*  -----------   Poly Spectra Bin Average   -----------  */
/*  ----------------------------------------------------  */


double integrand_average_line(
    double *var,
    size_t dim,
    void *params);


double integrand_average_tri(
    double *var,
    size_t dim,
    void *params);



/*  ----------------------------------------------------  */
/*  -------------   PP Covariance Matrix   -------------  */
/*  ----------------------------------------------------  */


double integrand_average_cov_pp_gauss(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_average_cov_pp_ngauss(
    double *var,
    size_t dim,
    void *params);

double integrand_average_cov_pp_ngauss_inf(
    double *var,
    size_t dim,
    void *params);



/*  ----------------------------------------------------  */
/*  -------------   BB Covariance Matrix   -------------  */
/*  ----------------------------------------------------  */


double integrand_average_cov_bb_gauss(
    double *var,
    size_t dim,
    void *params);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double integrand_average_cov_bb_ngauss_tp(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_average_cov_bb_ngauss_bb(
    double *var,
    size_t dim,
    void *params);



/*  ----------------------------------------------------  */
/*  -------------   PB Covariance Matrix   -------------  */
/*  ----------------------------------------------------  */


double integrand_average_cov_pb_ngauss_bp(
    double *var,
    size_t dim,
    void *params);



#endif // INTEGRAND_AVERAGE_H_INCLUDED

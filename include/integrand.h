#ifndef INTEGRAND_H_INCLUDED
#define INTEGRAND_H_INCLUDED

#include "common.h"

#include "misc.h"

#include "fiducials.h"

#include "kernels.h"
#include "interpolate.h"



/*  ----------------------------------------------------  */
/*  ----------   Non-Linear Power Spectrum   -----------  */
/*  ----------------------------------------------------  */


double integrand_spec_pnl_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_pnl_p22_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_pnl_p13_1loop(
    double *var,
    size_t dim,
    void *params);


/*  ----------------------------------------------------  */
/*  ---   Power Spectrum _1loop(Analytical) Derivatives   ----  */
/*  ----------------------------------------------------  */


double integrand_spec_dpnl_a2ga_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_a2ga_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_a2ga_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_d2ga_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_d2ga_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_d2ga_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_h_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_h_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_h_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_a3gaa_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_a3gaa_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_a3gaa_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_d3gaa_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_d3gaa_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_d3gaa_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_a3gab_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_a3gab_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_a3gab_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_d3gab_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_d3gab_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_d3gab_1loop(
    double *var,
    size_t dim,
    void *params);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double integrand_spec_dpnl_b1_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_b1_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_b1_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_f_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_f_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_f_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_b2_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_b2_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_b2_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_c2ga_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_c2ga_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_c2ga_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_bgam3_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_bgam3_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_bgam3_1loop(
    double *var,
    size_t dim,
    void *params);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double integrand_spec_dpnl_k_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_k_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_k_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_mu_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_mu_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_mu_1loop(
    double *var,
    size_t dim,
    void *params);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double integrand_spec_dpnl_pk_ang_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp22_pk_ang_1loop(
    double *var,
    size_t dim,
    void *params);

double integrand_spec_dpnl_dp13_pk_ang_1loop(
    double *var,
    size_t dim,
    void *params);

/*  ----------------------------------------------------  */

double integrand_spec_dpnl_pk_vol_1loop(
    double *var,
    size_t dim,
    void *params);


#endif // INTEGRAND_H_INCLUDED

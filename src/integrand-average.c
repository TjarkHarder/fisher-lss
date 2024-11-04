/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTEGRAND-AVERAGE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "integrand-average.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     POLY SPECTRA BIN AVERAGE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_line(double *var, size_t dim, void *params)
{
    /*

        Average over line-bins

    */

    (void) dim;

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = (spec_arg_t*) avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration variables */
    double k1 = var[0];
    double mu1 = var[1];

    /* Set variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);

    kernels_qset_nu(kern, 0, 1, -1.);

    /* Result */
    double result = 2.*M_PI * k1*k1 * avrParams -> specFunc(specArg, avrParams -> params);

    return result;
}


double integrand_average_tri(double *var, size_t dim, void *params)
{
    /*

        Average over triangle-bins

    */

    (void) dim;

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration bounds */
    double *kBounds = avrParams -> bounds;
    double *muBounds = avrParams -> bounds + 2;

    /* Integration variables */
    double k1 = var[0];
    double k2 = var[1];

    double mu1 = var[2];
    double mu2 = var[3];

    double phi = var[4];

    /* Constraining variables */
    double nu12 = cos(phi) * sqrt( fabs((1. - mu1*mu1) * (1. - mu2*mu2)) ) + mu1*mu2;
    double k12 = sqrt(k1*k1 + k2*k2 + 2.*k1*k2*nu12);
    double mu12 = - (k1*mu1 + k2*mu2) / k12;

    if ((kBounds[0] > k12 || kBounds[1] < k12) || (muBounds[0] > mu12 || muBounds[1] < mu12))
        return 0.;

    /* Set variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k12);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu12);

    kernels_qset_nu(kern, 0, 1, nu12);
    kernels_qset_nu(kern, 0, 2, -(k1 + k2*nu12) / k12);
    kernels_qset_nu(kern, 1, 2, -(k2 + k1*nu12) / k12);

    /* Result */
    double result = 2.*M_PI * k1*k1 * k2*k2 * avrParams -> specFunc(specArg, avrParams -> params);

    return result;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COVARIANCE MATRIX     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------    PP Covariance Matrix    -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_cov_pp_gauss(double *var, size_t dim, void *params)
{
    /*

        Integrand of the gaussian contribution to the PP covariance matrix.

        Assumes to have variables in kern struct:

            k = {_, _}
            mu = {_, _}
            nu = {-1}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration variables */
    double k = var[0];
    double mu = var[1];

    /* Set variables */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);
//    kernels_qset_nu(kern, 0, 1, -1.); // Should have already been set

    /* Power Spectrum */
    double pnl = _specPnl_(specArg, NULL);

    /* Result */
    double result = 2.*M_PI * k*k * pnl*pnl;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_cov_pp_ngauss(double *var, size_t dim, void *params)
{
    /*

        Integrand of the non-gaussian contribution to the PP covariance matrix (bin-averaged trispectrum).

        Assumes to have variables in kern struct:

            k = {_, _, _, _}
            mu = {_, _, _, _}
            nu = {-1, _, _, _, _, -1}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration variables */
    double k1 = var[0];
    double k2 = var[1];

    double mu1 = var[2];
    double mu2 = var[3];

    double phi = var[4];

    /* Additional variables */
    double nu = sqrt( fabs((1. - mu1*mu1) * (1. - mu2*mu2)) ) * cos(phi) + mu1*mu2;


    /* Set variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_k(kern, 2, k2);
    kernels_qset_k(kern, 3, k2);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_mu(kern, 2, mu2);
    kernels_qset_mu(kern, 3, -mu2);

//    kernels_qset_nu(kern, 0, 1, -1.); // Should have already been set
    kernels_qset_nu(kern, 0, 2, nu);
    kernels_qset_nu(kern, 0, 3, -nu);
    kernels_qset_nu(kern, 1, 2, -nu);
    kernels_qset_nu(kern, 1, 3, nu);
//    kernels_qset_nu(kern, 2, 3, -1.); // Should have already been set

    /* Result */
    double result = 2.*M_PI * k1*k1 * k2*k2 * _specTtr_(specArg, NULL);

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_cov_pp_ngauss_inf(double *var, size_t dim, void *params)
{
    /*

        Integrand of the non-gaussian contribution to the PP covariance matrix (bin-averaged trispectrum) for
        infinitesimally small bins in k and mu.

        Assumes to have variables in kern struct:

            k = {k, k, k', k'}
            mu = {mu, -mu, mu', -mu'}
            nu = {-1, _, _, _, _, -1}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration variables */
    double phi = var[0];

    /* Additional variables */
    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 2);

    double nu = sqrt( fabs((1. - mu1*mu1) * (1. - mu2*mu2)) ) * cos(phi) + mu1*mu2;

    /* Set variables */
//    kernels_qset_nu(kern, 0, 1, -1.); // Should have already been set
    kernels_qset_nu(kern, 0, 2, nu);
    kernels_qset_nu(kern, 0, 3, -nu);
    kernels_qset_nu(kern, 1, 2, -nu);
    kernels_qset_nu(kern, 1, 3, nu);
//    kernels_qset_nu(kern, 2, 3, -1.); // Should have already been set

    /* Result */
    double result = _specTtr_(specArg, NULL);

    return result;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------    BB Covariance Matrix    -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_cov_bb_gauss(double *var, size_t dim, void *params)
{
    /*

        Integrand of the gaussian contribution to the BB covariance matrix.

        Assumes to have variables in kern struct:

            k = {_, _}
            mu = {_, _}
            nu = {-1}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration bounds */
    double *kBounds = avrParams -> bounds;
    double *muBounds = avrParams -> bounds + 2;

    /* Integration variables */
    double k1 = var[0];
    double k2 = var[1];

    double mu1 = var[2];
    double mu2 = var[3];

    double phi = var[4];

    /* Constraining variables */
    double nu12 = cos(phi) * sqrt( fabs((1. - mu1*mu1) * (1. - mu2*mu2)) ) + mu1*mu2;
    double k12 = sqrt(k1*k1 + k2*k2 + 2.*k1*k2*nu12);
    double mu12 = - (k1*mu1 + k2*mu2) / k12;

    if ((kBounds[0] > k12 || kBounds[1] < k12) || (muBounds[0] > mu12 || muBounds[1] < mu12))
        return 0.;

    /* Results */

    /* Pnl1 */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);

//    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl1 = _specPnl_(specArg, NULL);

    /* Pnl2 */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);

    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);

//    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl2 = _specPnl_(specArg, NULL);

    /* Pnl3 */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k12);

    kernels_qset_mu(kern, 0, mu12);
    kernels_qset_mu(kern, 1, -mu12);

//    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl3 = _specPnl_(specArg, NULL);


    /* Result */
    double result = 2.*M_PI * k1*k1 * k2*k2 * pnl1 * pnl2 * pnl3;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_cov_bb_ngauss_tp(double *var, size_t dim, void *params)
{
    /*

        Integrand of the tp part of the non-gaussian contribution to the BB covariance matrix.

        Assumes to have variables in kern struct:

            k = {_, _, _, _}
            mu = {_, _, _, _}
            nu = {_, _, _, _, _, _}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration bounds */
    double *kBounds1 = avrParams -> bounds;
    double *muBounds1 = avrParams -> bounds + 2;
    double *kBounds2 = avrParams -> bounds + 4;
    double *muBounds2 = avrParams -> bounds + 6;

    /* Integration variables */
    double k1_1 = var[0];
    double k1_2 = var[1];
    double k2_2 = var[2];

    double mu1_1 = var[3];
    double mu1_2 = var[4];
    double mu2_2 = var[5];

    double phi1_2 = var[6];
    double phi2_2 = var[7];

    /* Constraining variables */
    double nu1_12 = cos(phi1_2) * sqrt( fabs((1. - mu1_1*mu1_1) * (1. - mu1_2*mu1_2)) ) + mu1_1*mu1_2;
    double k1_12 = sqrt(k1_1*k1_1 + k1_2*k1_2 + 2.*k1_1*k1_2*nu1_12);
    double mu1_12 = - (k1_1*mu1_1 + k1_2*mu1_2) / k1_12;

    if ((kBounds1[0] > k1_12 || kBounds1[1] < k1_12) || (muBounds1[0] > mu1_12 || muBounds1[1] < mu1_12))
        return 0.;

    double nu12_12 = cos(phi2_2) * sqrt( fabs((1. - mu1_1*mu1_1) * (1. - mu2_2*mu2_2)) ) + mu1_1*mu2_2;
    double k12_12 = sqrt(k1_1*k1_1 + k2_2*k2_2 + 2.*k1_1*k2_2*nu12_12);
    double mu12_12 = - (k1_1*mu1_1 + k2_2*mu2_2) / k12_12;

    if ((kBounds2[0] > k12_12 || kBounds2[1] < k12_12) || (muBounds2[0] > mu12_12 || muBounds2[1] < mu12_12))
        return 0.;

    /* Additional variables */
    double nu12_22 = cos(phi1_2 - phi2_2) * sqrt( fabs(1. - mu1_2*mu1_2) * fabs(1. - mu2_2*mu2_2) ) + mu1_2*mu2_2;


    /* Results */

    /* Ttr */
    kernels_qset_k(kern, 0, k1_2);
    kernels_qset_k(kern, 1, k1_12);
    kernels_qset_k(kern, 2, k2_2);
    kernels_qset_k(kern, 3, k12_12);

    kernels_qset_mu(kern, 0, mu1_2);
    kernels_qset_mu(kern, 1, mu1_12);
    kernels_qset_mu(kern, 2, -mu2_2);
    kernels_qset_mu(kern, 3, -mu12_12);

    kernels_qset_nu(kern, 0, 1, -(nu1_12 * k1_1 + k1_2) / k1_12);
    kernels_qset_nu(kern, 0, 2, -nu12_22);
    kernels_qset_nu(kern, 0, 3, (nu1_12 * k1_1 + nu12_22 * k2_2) / k12_12);
    kernels_qset_nu(kern, 1, 2, (nu12_12 * k1_1 + nu12_22 * k1_2) / k1_12);
    kernels_qset_nu(kern, 1, 3, -(k1_1 * k1_1 + k1_1 * k1_2 * nu1_12 + k1_1 * k2_2 * nu12_12 + k1_2 * k2_2 * nu12_22) / (k1_12 * k12_12));
    kernels_qset_nu(kern, 2, 3, -(k1_1 * nu12_12 + k2_2) / k12_12);

    double ttr = _specTtr_(specArg, NULL);

    /* Pnl */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_1);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, -mu1_1);

    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl = _specPnl_(specArg, NULL);

    /* Result */
    double result = 2.*M_PI * k1_1*k1_1 * k1_2*k1_2 * k2_2*k2_2 * ttr * pnl;

    return result;
}


double integrand_average_cov_bb_ngauss_bb(double *var, size_t dim, void *params)
{
    /*

        Integrand of the bb part of the non-gaussian contribution to the BB covariance matrix.

        Assumes to have variables in kern struct:

            k = {_, _, _}
            mu = {_, _, _}
            nu = {_, _, _}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration bounds */
    double *kBounds1 = avrParams -> bounds;
    double *muBounds1 = avrParams -> bounds + 2;
    double *kBounds2 = avrParams -> bounds + 4;
    double *muBounds2 = avrParams -> bounds + 6;

    /* Integration variables */
    double k1_1 = var[0];
    double k1_2 = var[1];
    double k2_2 = var[2];

    double mu1_1 = var[3];
    double mu1_2 = var[4];
    double mu2_2 = var[5];

    double phi1_2 = var[6];
    double phi2_2 = var[7];

    /* Constraining variables */
    double nu1_12 = cos(phi1_2) * sqrt( fabs((1. - mu1_1*mu1_1) * (1. - mu1_2*mu1_2)) ) + mu1_1*mu1_2;
    double k1_12 = sqrt(k1_1*k1_1 + k1_2*k1_2 + 2.*k1_1*k1_2*nu1_12);
    double mu1_12 = - (k1_1*mu1_1 + k1_2*mu1_2) / k1_12;

    if ((kBounds1[0] > k1_12 || kBounds1[1] < k1_12) || (muBounds1[0] > mu1_12 || muBounds1[1] < mu1_12))
        return 0.;

    double nu12_12 = cos(phi2_2) * sqrt( fabs((1. - mu1_1*mu1_1) * (1. - mu2_2*mu2_2)) ) + mu1_1*mu2_2;
    double k12_12 = sqrt(k1_1*k1_1 + k2_2*k2_2 - 2.*k1_1*k2_2*nu12_12);
    double mu12_12 = (k1_1*mu1_1 - k2_2*mu2_2) / k12_12;

    if ((kBounds2[0] > k12_12 || kBounds2[1] < k12_12) || (muBounds2[0] > mu12_12 || muBounds2[1] < mu12_12))
        return 0.;


    /* Results */

    /* Btr1 */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_2);
    kernels_qset_k(kern, 2, k1_12);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, mu1_2);
    kernels_qset_mu(kern, 2, mu1_12);

    kernels_qset_nu(kern, 0, 1, nu1_12);
    kernels_qset_nu(kern, 0, 2, (nu1_12 * k1_2 + k1_1) / k1_12);
    kernels_qset_nu(kern, 1, 2, (nu1_12 * k1_1 + k1_2) / k1_12);

    double btr1 = _specBtr_(specArg, NULL);

    /* Btr2 */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k2_2);
    kernels_qset_k(kern, 2, k12_12);

    kernels_qset_mu(kern, 0, -mu1_1);
    kernels_qset_mu(kern, 1, mu2_2);
    kernels_qset_mu(kern, 2, mu12_12);

    kernels_qset_nu(kern, 0, 1, -nu12_12);
    kernels_qset_nu(kern, 0, 2, (nu12_12 * k2_2 - k1_1) / k12_12);
    kernels_qset_nu(kern, 1, 2, (nu12_12 * k1_1 - k2_2) / k12_12);

    double btr2 = _specBtr_(specArg, NULL);


    /* Result */
    double result = 2.*M_PI * k1_1*k1_1 * k1_2*k1_2 * k2_2*k2_2 * btr1 * btr2;

    return result;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------    PB Covariance Matrix    -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_average_cov_pb_ngauss_bp(double *var, size_t dim, void *params)
{
    /*

        Integrand of the bispectrum x power spectrum part of the non-gaussian contribution to the PB covariance matrix.

        Assumes to have variables in kern struct:

            k = {_, _, _}
            mu = {_, _, _}
            nu = {_, _, _}

    */

    /* Not used */
    (void) (dim);

    /* Parameters */
    cintgrnd_avr_t *avrParams = (cintgrnd_avr_t*) params;

    /* specArg and kern structs */
    spec_arg_t *specArg = avrParams -> var;
    kern_t *kern = specArg -> kern;

    /* Integration bounds */
    double *kBounds = avrParams -> bounds;
    double *muBounds = avrParams -> bounds + 2;

    /* Integration variables */
    double k1_1 = var[0];
    double k2_2 = var[1];

    double mu1_1 = var[2];
    double mu2_2 = var[3];

    double phi2_2 = var[4];

    /* Constraining variables */
    double nu12_12 = cos(phi2_2) * sqrt( fabs((1. - mu1_1*mu1_1) * (1. - mu2_2*mu2_2)) ) + mu1_1*mu2_2;
    double k12_12 = sqrt(k1_1*k1_1 + k2_2*k2_2 + 2.*k1_1*k2_2*nu12_12);
    double mu12_12 = -(k1_1*mu1_1 + k2_2*mu2_2) / k12_12;

    if ((kBounds[0] > k12_12 || kBounds[1] < k12_12) || (muBounds[0] > mu12_12 || muBounds[1] < mu12_12))
        return 0.;


    /* Results */

    /* Pnl */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_1);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, -mu1_1);

    kernels_qset_nu(kern, 0, 1, -1);

    double pnl = _specPnl_(specArg, NULL);

    /* Btr */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k2_2);
    kernels_qset_k(kern, 2, k12_12);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, mu2_2);
    kernels_qset_mu(kern, 2, mu12_12);

    kernels_qset_nu(kern, 0, 1, nu12_12);
    kernels_qset_nu(kern, 0, 2, -(nu12_12 * k2_2 + k1_1) / k12_12);
    kernels_qset_nu(kern, 1, 2, -(nu12_12 * k1_1 + k2_2) / k12_12);

    double btr = _specBtr_(specArg, NULL);


    /* Result */
    double result = 2.*M_PI * k1_1*k1_1 * k2_2*k2_2 * pnl * btr;

    return result;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

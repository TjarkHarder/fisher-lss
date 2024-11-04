/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTEGRAND.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "integrand.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     NON - LINEAR POWER SPECTRUM     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------    Non-Linear Power Spectrum    ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_pnl_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of Pnl(k_):

            q^2 P^(0)(q) ( (Z2(k_-q_, q_) Z2(k_-q_, q_) P^(0)(|k_-q_|)
                                - Z2(-q_, q_) Z2(-q_, q_) P^(0)(q))
                           + 3 Z3(k_, -q_, q_) Z1(k_) P^(0)(k)) ).


    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_pnl_p22_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_pnl_p13_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_pnl_p22_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of P_22(k_):

            q^2 P^(0)(q) ( Z2(k_-q_, q_) Z2(k_-q_, q_) P^(0)(|k_-q_|)
                            - Z2(-q_, q_) Z2(-q_, q_) P^(0)(q) ).
                              \____  ___/ \____  ___/
                                   \/          \/
                                b2(z)/2      b2(z)/2      ->    renormalisation

    */


    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pq = _fidPk_(&q, _fidParamsPk_);
    double pkq = _fidPk_(&kq, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);

    kernels_qset_k(kern, 0, q);
    kernels_qset_mu(kern, 0, -muq);
    kernels_qset_nu(kern, 0, 1, -1.);

    double z2Re = kernels_z2(kern);

    /* Compute the integrand */
    double result = q*q * pq * (z2*z2 * pkq - z2Re*z2Re * pq);

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_pnl_p13_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of P_13(k_):

            3 q^2 Z3(k_, -q_, q_) Z1(k_) P^(0)(q) P^(0)(k).

    */


    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double z3 = kernels_z3(kern);

    /* Integrand */
    double result = 3. * q*q * pq * z3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------    Non-Linear Power Spectrum (Analytical) Derivatives    ---------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_a2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(a^(2)_γ):

            q^2 (2 δZ2(k_-q_, q_)/δ(a^(2)_γ) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(a^(2)_γ) P^(0)(q) Z1(k_) P^(0)(k))

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_a2ga_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_a2ga_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_a2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(a^(2)_γ):

            2 q^2 δZ2(k_-q_, q_)/δ(a^(2)_γ) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pq = _fidPk_(&q, _fidParamsPk_);
    double pkq = _fidPk_(&kq, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_a2ga(kern);

    /* Compute the integrand */
    double result = 2. * q*q * pq * z2*dz2 * pkq;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_a2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(a^(2)_γ):

            3 q^2 δZ3(k_, q_, -q_)/δ(a^(2)_γ) P^(0)(q) Z1(k_) P^(0)(k)

        Note: It is assumed that h = a^(2)_γ - 1 => δh/δa^(2)_γ = 1

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_a2ga(kern);
    double dz3_ = kernels_dz3_h(kern);

    /* Integrand */
    double result = 3. * q*q * pq * (dz3 + dz3_) * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_d2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(d^(2)_γ):

            q^2 (2 δZ2(k_-q_, q_)/δ(d^(2)_γ) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(d^(2)_γ) P^(0)(q) Z1(k_) P^(0)(k))

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_d2ga_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_d2ga_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_d2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(d^(2)_γ):

            2 q^2 δZ2(k_-q_, q_)/δ(d^(2)_γ) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_d2ga(kern);


    /* Compute the integrand */
    double result = 2. * q*q * pq * z2*dz2 * pkq;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_d2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(d^(2)_γ):

            3 q^2 δZ3(k_, q_, -q_)/δ(d^(2)_γ) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_d2ga(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_h_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(h):

            3 δZ3(k_, q_, -q_)/δ(h) P^(0)(q) Z1(k_) P^(0)(k))

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_h_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_h_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_h_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(h):

            2 q^2 δZ2(k_-q_, q_)/δ(a^(3)_γa) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|) = 0.

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0;
}


double integrand_spec_dpnl_dp13_h_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(h):

            3 q^2 δZ3(k_, q_, -q_)/δ(h) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_h(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_a3gaa_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(a^(3)_γa):

            q^2 (2 δZ2(k_-q_, q_)/δ(a^(3)_γa) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(a^(3)_γa) P^(0)(q) Z1(k_) P^(0)(k))

          = 3 q^2 δZ3(k_, q_, -q_)/δ(a^(3)_γa) P^(0)(q) Z1(k_) P^(0)(k)

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_a3gaa_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_a3gaa_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_a3gaa_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(a^(3)_γa):

            2 q^2 δZ2(k_-q_, q_)/δ(a^(3)_γa) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|) = 0.

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0.;
}


double integrand_spec_dpnl_dp13_a3gaa_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(a^(3)_γa):

            3 q^2 δZ3(k_, q_, -q_)/δ(a^(3)_γa) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_a3gaa(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_d3gaa_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(d^(3)_γa):

            q^2 (2 δZ2(k_-q_, q_)/δ(d^(3)_γa) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(d^(3)_γa) P^(0)(q) Z1(k_) P^(0)(k))

          = 3 q^2 δZ3(k_, q_, -q_)/δ(d^(3)_γa) P^(0)(q) Z1(k_) P^(0)(k)

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_d3gaa_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_d3gaa_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_d3gaa_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(d^(3)_γa):

            2 q^2 δZ2(k_-q_, q_)/δ(d^(3)_γa) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|) = 0.

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0.;
}


double integrand_spec_dpnl_dp13_d3gaa_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(d^(3)_γa):

            3 q^2 δZ3(k_, q_, -q_)/δ(d^(3)_γa) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_d3gaa(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_a3gab_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(a^(3)_γb):

            q^2 (2 δZ2(k_-q_, q_)/δ(a^(3)_γb) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(a^(3)_γb) P^(0)(q) Z1(k_) P^(0)(k))

          = 3 q^2 δZ3(k_, q_, -q_)/δ(a^(3)_γb) P^(0)(q) Z1(k_) P^(0)(k)

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_a3gab_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_a3gab_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_a3gab_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(a^(3)_γb):

            2 q^2 δZ2(k_-q_, q_)/δ(a^(3)_γb) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|) = 0.

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0.;
}


double integrand_spec_dpnl_dp13_a3gab_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(a^(3)_γb):

            3 q^2 δZ3(k_, q_, -q_)/δ(a^(3)_γb) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_a3gab(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_d3gab_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(d^(3)_γb):

            q^2 (2 δZ2(k_-q_, q_)/δ(d^(3)_γb) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(d^(3)_γb) P^(0)(q) Z1(k_) P^(0)(k))

          = 3 q^2 δZ3(k_, q_, -q_)/δ(d^(3)_γb) P^(0)(q) Z1(k_) P^(0)(k)

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_d3gab_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_d3gab_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_d3gab_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(d^(3)_γb):

            2 q^2 δZ2(k_-q_, q_)/δ(d^(3)_γb) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|) = 0.

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0.;
}


double integrand_spec_dpnl_dp13_d3gab_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(d^(3)_γb):

            3 q^2 δZ3(k_, q_, -q_)/δ(d^(3)_γb) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_d3gab(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_b1_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl/δb1(k_):

            q^2 (2 δZ2/δb1(k_-q_, q_) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3/δb1(k_, q_, -q_) P^(0)(q) Z1(k_) P^(0)(k)
                + 3 Z3(k_, q_, -q_) P^(0)(q) δZ1/δb1(k_) P^(0)(k)).

        Note: No need for P22 renormalisation here since δZ2/δb1(-q_, q_) = 0

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_b1_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_b1_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_b1_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(b1):

            2 q^2 δZ2(k_-q_, q_)/δ(b1) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_b1(kern);


    /* Compute the integrand */
    double result = 2. * q*q * pq * z2*dz2 * pkq;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_b1_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(b1):

            3 q^2 P^(0)(q) (δZ3(k_, q_, -q_)/δ(b1) Z1(k_) + Z3(k_, q_, -q_) δZ1/δ(b1)(k_)) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz1 = kernels_dz1_b1(kern);

    double z3 = kernels_z3(kern);
    double dz3 = kernels_dz3_b1(kern);

    /* Integrand */
    double result = 3. * q*q * pq * (dz3 * z1 + z3 * dz1) * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_f_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl/δf(k_):

            q^2 (2 δZ2/δf(k_-q_, q_) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3/δf(k_, q_, -q_) P^(0)(q) Z1(k_) P^(0)(k)
                + 3 Z3(k_, q_, -q_) P^(0)(q) δZ1/δf(k_) P^(0)(k)).

        Note: No need for P22 renormalisation here since δZ2/δf(-q_, q_) = 0

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_f_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_f_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_f_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(f):

            2 q^2 δZ2(k_-q_, q_)/δ(f) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_f(kern);


    /* Compute the integrand */
    double result = 2. * q*q * pq * z2*dz2 * pkq;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_f_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(f):

            3 q^2 P^(0)(q) (δZ3(k_, q_, -q_)/δ(f) Z1(k_) + Z3(k_, q_, -q_) δZ1/δ(f)(k_)) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz1 = kernels_dz1_f(kern);

    double z3 = kernels_z3(kern);
    double dz3 = kernels_dz3_f(kern);

    /* Integrand */
    double result = 3. * q*q * pq * (dz3 * z1 + z3 * dz1) * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_b2_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl/δb2(k_):

            2 q^2 P^(0)(q)(δZ2/δb2(k_-q_, q_) Z2(k_-q_, q_) P^(0)(|k_-q_|)
                            - δZ2/δb2(-q_, q_) Z2(-q_, q_) P^(0)(q))
                              \______  ______/ \____  ___/
                                     \/             \/
                                     1/2           b2/2     ->   renormalisation

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_b2_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_b2_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_b2_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22/δb2(k_):

            2 q^2 P^(0)(q)(δZ2/δb2(k_-q_, q_) Z2(k_-q_, q_) P^(0)(|k_-q_|)
                            - δZ2/δb2(-q_, q_) Z2(-q_, q_) P^(0)(q))
                              \______  ______/ \____  ___/
                                     \/             \/
                                     1/2           b2/2     ->   renormalisation

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_b2(kern);

    kernels_qset_k(kern, 0, q);
    kernels_qset_mu(kern, 0, -muq);
    kernels_qset_nu(kern, 0, 1, -1.);

    double z2Re = kernels_z2(kern);
    double dz2Re = kernels_dz2_b2(kern);


    /* Compute the integrand */
    double result = 2. * q*q * pq * (z2*dz2 * pkq - z2Re*dz2Re * pq);

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_b2_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(b2):

            3 q^2 P^(0)(q) δZ3(k_, q_, -q_)/δ(b2) Z1(k_) P^(0)(k) = 0

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0.;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_c2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(c^(2)_γ):

            q^2 (2 δZ2(k_-q_, q_)/δ(c^(2)_γ) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(c^(2)_γ) P^(0)(q) Z1(k_) P^(0)(k))

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_c2ga_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_c2ga_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_c2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(c^(2)_γ):

            2 q^2 δZ2(k_-q_, q_)/δ(c^(2)_γ) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_c2ga(kern);


    /* Compute the integrand */
    double result = 2. * q*q * pq * z2*dz2 * pkq;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_c2ga_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(c^(2)_γ):

            3 q^2 δZ3(k_, q_, -q_)/δ(c^(2)_γ) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_c2ga(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_bgam3_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δ(b_Γ3):

            q^2 (2 δZ2(k_-q_, q_)/δ(b_Γ3) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|)
                + 3 δZ3(k_, q_, -q_)/δ(b_Γ3) P^(0)(q) Z1(k_) P^(0)(k))
          = 3 q^2 δZ3(k_, q_, -q_)/δ(b_Γ3) P^(0)(q) Z1(k_) P^(0)(k)

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_bgam3_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_bgam3_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_bgam3_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δ(b_Γ3):

            2 q^2 δZ2(k_-q_, q_)/δ(b_Γ3) Z2(k_-q_, q_) P^(0)(q) P^(0)(|k_-q_|) = 0

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;
    (void) var;
    (void) params;

    return 0.;
}


double integrand_spec_dpnl_dp13_bgam3_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δ(b_Γ3):

            3 q^2 δZ3(k_, q_, -q_)/δ(b_Γ3) P^(0)(q) Z1(k_) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz3 = kernels_dz3_bgam3(kern);

    /* Integrand */
    double result = 3. * q*q * pq * dz3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_k_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δk:

            q^2 P^(0)(q) (2 Z2(k_-q_, q_) δZ2/δk(k_-q_, q_) P^(0)(|k_-q_|)
                            + (Z2(k_-q_, q_))^2 δP^(0)/δk(|k_-q_|)
                            + 3 δZ3/δk(k_, q_, -q_) Z1(k_) P^(0)(k)
                            + 3 Z3(k_, q_, -q_) Z1(k_) δP^(0)/δk(k))

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_k_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_k_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_k_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δk:

            q^2 P^(0)(q) (2 Z2(k_-q_, q_) δZ2/δk(k_-q_, q_) P^(0)(|k_-q_|)
                            + (Z2(k_-q_, q_))^2 δP^(0)/δk(|k_-q_|))

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    double dkq = (k - q*nu) / kq;
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    double dmukq = (mu - mukq * dkq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    double dnukq = (nu - nukq * dkq) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double dpkq = _fidDPk_(&kq, _fidParamsDPk_) * dkq;
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_k(kern) * dkq + kernels_dz2_nu(kern) * dnukq + kernels_dz2_mu(kern) * dmukq;


    /* Compute the integrand */
    double result = q*q * pq * z2 * (2. * dz2 * pkq + z2 * dpkq);

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_k_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δk:

            3 q^2 P^(0)(q) Z1(k_) (δZ3(k_, q_, -q_)/δk P^(0)(k) + Z3(k_, q_, -q_) δP^(0)(k)/δk)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double dpk = _fidDPk_(&k, _fidParamsDPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double z3 = kernels_z3(kern);
    double dz3 = kernels_dz3_k(kern);

    /* Integrand */
    double result = 3. * q*q * pq * (dz3 * pk + z3 * dpk) * z1;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_mu_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δPnl(k_)/δmu:

            q^2 P^(0)(q) (2 Z2(k_-q_, q_) δZ2/δmu(k_-q_, q_) P^(0)(|k_-q_|)
                            + 3 δZ3/δmu(k_, q_, -q_) Z1(k_) P^(0)(k)
                            + 3 Z3(k_, q_, -q_) δZ1(k_)/δmu P^(0)(k))

    */


    /* Not used */
    (void) dim;

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_mu_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_mu_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_dp22_mu_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP22(k_)/δmu:

            2 q^2 P^(0)(q) Z2(k_-q_, q_) δZ2/δmu(k_-q_, q_) P^(0)(|k_-q_|)

        Note: No need for P22 renormalisation here

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    /* |q_| */
    kernels_qset_k(kern, 1, q);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    double dmuq = -mu * sqrt( (1. - nu*nu) / (1. - mu*mu) ) * cphi + nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    double dmukq = (k - q*dmuq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);
    double dz2 = kernels_dz2_mu(kern) * dmukq;

      {
        /*

            Derivative w.r.t to muq -> (q_, k_-q_)

        */

        kernels_qset_k(kern, 0, q);
        kernels_qset_k(kern, 1, kq);

        kernels_qset_mu(kern, 0, muq);
        kernels_qset_mu(kern, 1, mukq);

        dz2 += kernels_dz2_mu(kern) * dmuq;
      }

    /* Compute the integrand */
    double result = 2. * q*q * pq * z2*dz2 * pkq;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_mu_1loop(double *var, size_t dim, void *params)
{
    /*

        Integrand of δP13(k_)/δmu:

            3 q^2 P^(0)(q) (δZ3(k_, q_, -q_)/δmu Z1(k_) + Z3(k_, q_, -q_) δZ1(k_)/δk) P^(0)(k)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);

    /* dmuq / dmu */
    double dmuq = -mu * sqrt( (1. - nu*nu) / (1. - mu*mu) ) * cphi + nu;


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double dz1 = kernels_dz1_mu(kern);

    double z3 = kernels_z3(kern);
    double dz3 = kernels_dz3_mu(kern);

      {
        /*

            Derivative w.r.t to muq -> (q_, k_, -q_)

        */

        kernels_qset_k(kern, 0, q);
        kernels_qset_k(kern, 1, k);
        kernels_qset_k(kern, 2, q);

        kernels_qset_mu(kern, 0, muq);
        kernels_qset_mu(kern, 1, mu);
        kernels_qset_mu(kern, 2, -muq);

        kernels_qset_nu(kern, 0, 1, nu);
        kernels_qset_nu(kern, 0, 2, -1.);
        kernels_qset_nu(kern, 1, 2, -nu);

        dz3 += kernels_dz3_mu(kern) * dmuq;
      }

      {
        /*

            Derivative w.r.t to -muq -> (-q_, k_, q_)

            Note: Must treat -muq as independent variable and cannot simply do dz3 += 2 * (...)!

        */

        kernels_qset_k(kern, 0, q);
        kernels_qset_k(kern, 1, k);
        kernels_qset_k(kern, 2, q);

        kernels_qset_mu(kern, 0, -muq);
        kernels_qset_mu(kern, 1, mu);
        kernels_qset_mu(kern, 2, muq);

        kernels_qset_nu(kern, 0, 1, -nu);
        kernels_qset_nu(kern, 0, 2, -1.);
        kernels_qset_nu(kern, 1, 2, nu);

        dz3 -= kernels_dz3_mu(kern) * dmuq;
    }

    /* Integrand */
    double result = 3. * q*q * pq * (dz3 * z1 + z3 * dz1) * pk;

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    /* Reset b2 */
    kern -> bias -> b2 = b2;

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_pk_ang_1loop(double *var, size_t dim, void *params)
{
    /*

        (Angular) Integrand of δPnl(k_)/δP^(0)(k'):

             2 ( (Z2(k_-q_, q_))^2 P^(0)(|k_-q_|) - (Z2(-q_, q_))^2 P^(0)(k) )
              + 3 Z3(k_, q_, -q_) Z1(k_) P^(0)(k) ).

        where q = k'.

    */


    /* Not used */
    (void) (dim);

    /* P22 integrand */
    double integrandP22 = integrand_spec_dpnl_dp22_pk_ang_1loop(var, dim, params);

    /* P13 integrand */
    double integrandP13 = integrand_spec_dpnl_dp13_pk_ang_1loop(var, dim, params);

    /* Combine the terms */
    double result = integrandP22 + integrandP13;

    return result;
}


double integrand_spec_dpnl_dp22_pk_ang_1loop(double *var, size_t dim, void *params)
{
    /*

        (Angular) Integrand of δP22(k_)/δP^(0)(k'):

             2 ( (Z2(k_-q_, q_))^2 P^(0)(|k_-q_|) - (Z2(-q_, q_))^2 P^(0)(q) )

        where q = k'. The factor of two is due to the symmetry q_ <-> k_ - q_.

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    double q = kernels_qget_k(kern, 1);

    /* Integration variables */
    double nu = var[0];
    double cphi = cos(var[1]);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* |k_ - q_| */
    double kq = sqrt( q*q + k*k - 2.*k*q*nu );
    kernels_qset_k(kern, 0, kq);

    /* (k_ - q_).s_ / |k_ - q_| */
    double mukq = (k*mu - q*muq) / kq;
    kernels_qset_mu(kern, 0, mukq);

    /* (k_ - q_).q_ / (|k_ - q_| q) */
    double nukq = (k*nu - q) / kq;
    kernels_qset_nu(kern, 0, 1, nukq);


    /* Contributions */

    /* Power spectrum */
    double pkq = _fidPk_(&kq, _fidParamsPk_);
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Compute kernels */
    double z2 = kernels_z2(kern);

    kernels_qset_k(kern, 0, q);
    kernels_qset_mu(kern, 0, -muq);
    kernels_qset_nu(kern, 0, 1, -1.);

    double z2Re = kernels_z2(kern);


    /* Compute the integrand */
    double result = 2. * (z2*z2 * pkq - z2Re*z2Re * pq);

    /* Reset kern */
    kernels_qset_k(kern, 0, k);
    kernels_qset_mu(kern, 0, mu);

    return result;
}


double integrand_spec_dpnl_dp13_pk_ang_1loop(double *var, size_t dim, void *params)
{
    /*

        (Angular) Integrand of δP13(k_)/δP^(0)(k'):

              3 Z3(k_, q_, -q_) Z1(k_) P^(0)(k)

        where q = k'.

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    double q = kernels_qget_k(kern, 1);

    /* Integration variables */
    double nu = var[0];
    double cphi = cos(var[1]);

    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pk = _fidPk_(&k, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double z3 = kernels_z3(kern);

    /* Integrand */
    double result = 3. * z3 * z1 * pk;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double integrand_spec_dpnl_pk_vol_1loop(double *var, size_t dim, void *params)
{
    /*

        (Volume) Integrand of δPnl(k_)/δP^(0)(k'):

              3 q^2 Z1(k_) Z3(k_, q_, -q_) P^(0)(q)

    */

    /* Not used */
    (void) dim;

    /* Kern parameters */
    kern_t *kern = (kern_t*) params;

    /* Renormalise b2 -> 0 */
    double b2 = kern -> bias -> b2;
    kern -> bias -> b2 = 0.;


    /* Use correct variables */

    /* Scales and angles at which loop integral is performed */
//    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    /* Integration variables */
    double q = var[0];
    double nu = var[1];
    double cphi = cos(var[2]);

    kernels_qset_k(kern, 1, q);
    kernels_qset_k(kern, 2, q);

    kernels_qset_nu(kern, 0, 1, nu);
    kernels_qset_nu(kern, 0, 2, -nu);
    kernels_qset_nu(kern, 1, 2, -1.);

    /* q_.s_ / q */
    double muq = sqrt( (1. - mu*mu) * (1. - nu*nu) ) * cphi + mu*nu;
    kernels_qset_mu(kern, 1, muq);

    /* -q_.s_ / q */
    kernels_qset_mu(kern, 2, -muq);


    /* Contributions */

    /* Power spectrum */
    double pq = _fidPk_(&q, _fidParamsPk_);

    /* Kernels */
    double z1 = kernels_z1(kern);
    double z3 = kernels_z3(kern);

    /* Integrand */
    double result = 3. * q*q * pq * z3 * z1;

    /* Reset b2 */
    kern -> bias -> b2 = b2;


    return result;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

#ifndef KERNELS_H_INCLUDED
#define KERNELS_H_INCLUDED

#include "common.h"

#include "misc.h"

#include "fiducials.h"



/*  ----------------------------------------------------  */
/*  ----------------   Kernel Structure   --------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Kernel parameters

     */

    /* Spectrum order (2 for power spectrum, 3 for bispectrum, ...) */
    size_t specOrder;


    /* Kernel order TODO: Needed? */
    size_t kernOrder;

    /* Max order which dictates the amount of memory used for k, mu and nu */
    size_t maxOrder;

    /* Variables */
    double z;

    double *k;
    double *nu;
    double *mu;

    /* Linear growth factor */
    double growth;

    /* Fiducials */
    fid_lcdm_t *lcdm;

    fid_btst_t *btst;

    fid_bias_t *bias;
    fid_rsd_t *rsd;

    fid_ctr_t *ctr;

    fid_surv_t *surv;


    /* Internal kern_t structs */
    void *kernWork;
    bool computeWork;

} kern_t;



/*  ----------------------------------------------------  */
/*  ------------------   Kern Struct   -----------------  */
/*  ----------------------------------------------------  */


kern_t *kernels_new(
    size_t specOrder,
    size_t loopOrder);

kern_t *kernels_new_order(
    size_t specOrder,
    size_t kernOrder);

kern_t *kernels_free(
    kern_t *kernVar);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int kernels_set_z(
    kern_t *kern,
    double z);


int kernels_set_k(
    kern_t *kern,
    size_t index,
    double k);

int kernels_qset_k(
    kern_t *kern,
    size_t index,
    double k);


int kernels_set_mu(
    kern_t *kern,
    size_t index,
    double mu);

int kernels_qset_mu(
    kern_t *kern,
    size_t index,
    double mu);


int kernels_set_nu(
    kern_t *kern,
    size_t index1,
    size_t index2,
    double nu);

int kernels_qset_nu(
    kern_t *kern,
    size_t index1,
    size_t index2,
    double nu);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double kernels_get_z(
    kern_t *kern);


double kernels_get_k(
    kern_t *kern,
    size_t index);

double kernels_qget_k(
    kern_t *kern,
    size_t index);


double kernels_get_mu(
    kern_t *kern,
    size_t index);

double kernels_qget_mu(
    kern_t *kern,
    size_t index);


double kernels_get_nu(
    kern_t *kern,
    size_t index1,
    size_t index2);

double kernels_qget_nu(
    kern_t *kern,
    size_t index1,
    size_t index2);


size_t kernels_get_nu_index(
    size_t order,
    size_t index1,
    size_t index2);



/*  ----------------------------------------------------  */
/*  --------------   Raw bootstrap Kernels   -----------  */
/*  ----------------------------------------------------  */


double kernels_alpha(
    kern_t *kernVar);

double kernels_dalpha_k(
    kern_t *kernVar);

double kernels_dalpha_nu(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_beta(
    kern_t *kernVar);

double kernels_dbeta_k(
    kern_t *kernVar);

double kernels_dbeta_nu(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_gamma(
    kern_t *kernVar);

double kernels_gamma_nu(
    kern_t *kernVar);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int kernels_btst_h2(
    kern_t *kernVar,
    double *h2Kernels);

/*  ----------------------------------------------------  */

int kernels_btst_dh2_k(
    kern_t *kernVar,
    double *dh2Kernels);

int kernels_btst_dh2_nu(
    kern_t *kernVar,
    double *dh2Kernels);

/*  ----------------------------------------------------  */

int kernels_btst_dh2_a2ga(
    kern_t *kernVar,
    double *dh2Kernels);

int kernels_btst_dh2_d2ga(
    kern_t *kernVar,
    double *dh2Kernels);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int kernels_btst_h3(
    kern_t *kernVar,
    double *h3Kernels);

/*  ----------------------------------------------------  */

int kernels_btst_dh3_k(
    kern_t *kernVar,
    double *dh3Kernels);

int kernels_btst_dh3_nu(
    kern_t *kernVar,
    double *dh3Kernels);

/*  ----------------------------------------------------  */

int kernels_btst_dh3_a2ga(
    kern_t *kernVar,
    double *dh3Kernels);

int kernels_btst_dh3_d2ga(
    kern_t *kernVar,
    double *dh3Kernels);

int kernels_btst_dh3_h(
    kern_t *kernVar,
    double *dh3Kernels);


int kernels_btst_dh3_a3gaa(
    kern_t *kernVar,
    double *dh3Kernels);

int kernels_btst_dh3_a3gab(
    kern_t *kernVar,
    double *dh3Kernels);

int kernels_btst_dh3_d3gaa(
    kern_t *kernVar,
    double *dh3Kernels);

int kernels_btst_dh3_d3gab(
    kern_t *kernVar,
    double *dh3Kernels);



/*  ----------------------------------------------------  */
/*  ----------------   Redshift Kernels   --------------  */
/*  ----------------------------------------------------  */


double kernels_smooth(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dsmooth_sigv(
    kern_t *kernVar);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double kernels_z1(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz1_mu(
    kern_t *kernVar);

double kernels_dz1_b1(
    kern_t *kernVar);

double kernels_dz1_f(
    kern_t *kernVar);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double kernels_z2(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz2_k(
    kern_t *kernVar);

double kernels_dz2_nu(
    kern_t *kernVar);

double kernels_dz2_mu(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz2_a2ga(
    kern_t *kernVar);

double kernels_dz2_d2ga(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz2_b1(
    kern_t *kernVar);

double kernels_dz2_b2(
    kern_t *kernVar);

double kernels_dz2_c2ga(
    kern_t *kernVar);

double kernels_dz2_f(
    kern_t *kernVar);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double kernels_z3(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz3_k(
    kern_t *kernVar);

double kernels_dz3_nu(
    kern_t *kernVar);

double kernels_dz3_mu(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz3_a2ga(
    kern_t *kernVar);

double kernels_dz3_d2ga(
    kern_t *kernVar);

double kernels_dz3_h(
    kern_t *kernVar);


double kernels_dz3_a3gaa(
    kern_t *kernVar);

double kernels_dz3_a3gab(
    kern_t *kernVar);

double kernels_dz3_d3gaa(
    kern_t *kernVar);

double kernels_dz3_d3gab(
    kern_t *kernVar);

/*  ----------------------------------------------------  */

double kernels_dz3_b1(
    kern_t *kernVar);

double kernels_dz3_b2(
    kern_t *kernVar);

double kernels_dz3_c2ga(
    kern_t *kernVar);

double kernels_dz3_bgam3(
    kern_t *kernVar);

double kernels_dz3_f(
    kern_t *kernVar);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // KERNELS_H_INCLUDED

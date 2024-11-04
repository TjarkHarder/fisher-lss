#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include "common.h"

#include "misc.h"



/*  ----------------------------------------------------  */
/*  --------------   Integrate Struct   ----------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Parameters for the CQUAD routine by GSL

    */

    double relErr; // desired relative error
    double absErr; // desired absolute error

    size_t limit; // maximum number of subintervals

} intgrt_cquad_t;


typedef struct
{
    /*

        Parameters for the vegas routine by GSL

    */

    size_t mainCalls;
    size_t warmUpCalls;
    double updateCalls;

    int maxLoop;

    double chisqErr;
    double relErr;

    int verbose;

} intgrt_vegas_t;


typedef struct
{
    /*

        Parameters for the Divonne routine by CUBA

    */

    int nComp; // Number of components of integrand
    int nVec; // Maximum number of points given to the integrand routine

    double epsRel; // Relative accuracy goal
    double epsAbs; // Absolute accuracy goal

    int seed; // Seed of pseudo random generator

    int minEval; // Minimum number of integrand evaluations
    int maxEval; // (Approximate) Maximum number of integrand evaluations

    int key1; // Determines sampling in the partitioning phase
    int key2; // Determines sampling in the final integration phase
    int key3; // Sets strategy of the refinement phase

    int maxPass; // Thoroughness of the partitioning phase

    double border; // Width of the border of the integration region in which points are not sampled directly

    double maxChisq; // Maximum χ^2 value each subregion is allowed to have
    double minDeviation; // Determines if a region is further examined if it failed the χ^2 test (given as fraction of requested error of integral)

//    char *stateFile; // Filename for storing internal state

//    int *spin; // Spinning cores

    int verbose; // Verbose level of output

} intgrt_divonne_t;


typedef struct
{
    /*

        Parameters for the integral function

    */

    size_t dim; // Dimension of integral

    double *lowerBounds; // Lower bounds
    double *upperBounds; // Upper bounds

    intgrt_cquad_t *cquad; // CQUAD struct
    intgrt_vegas_t *vegas; // Vegas struct
    intgrt_divonne_t *divonne; // Divonne struct

    const char *routine;

    double (*integrand)(double*, size_t, void*); // Integrand
    void *params; // Parameters of the integrand

} intgrt_t;




/*  ----------------------------------------------------  */
/*  ---------------   External Variables   -------------  */
/*  ----------------------------------------------------  */


/**  Integration Routines  **/

extern const char *_idIntgrtCQUAD_;
extern const char *_idIntgrtVegas_;
extern const char *_idIntgrtDivonne_;



/*  ----------------------------------------------------  */
/*  --------------   Integration Routines   ------------  */
/*  ----------------------------------------------------  */


const char *integrate_find_routine(
    const char *routine);



/*  ----------------------------------------------------  */
/*  ------------------   CQUAD Struct   -----------------  */
/*  ----------------------------------------------------  */


intgrt_cquad_t *integrate_cquad_new(void);

intgrt_cquad_t *integrate_cquad_free(
    intgrt_cquad_t *cquad);

/*  ----------------------------------------------------  */

intgrt_cquad_t *integrate_cquad_cp(
    intgrt_cquad_t *cquad);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_cquad_set_abserr(
    intgrt_cquad_t *cquad,
    double abserr);

int integrate_cquad_set_relerr(
    intgrt_cquad_t *cquad,
    double relErr);

int integrate_cquad_set_limit(
    intgrt_cquad_t *cquad,
    size_t limit);



/*  ----------------------------------------------------  */
/*  ------------------   Vegas Struct   ----------------  */
/*  ----------------------------------------------------  */


intgrt_vegas_t *integrate_vegas_new(void);

intgrt_vegas_t *integrate_vegas_free(
    intgrt_vegas_t *vegas);

/*  ----------------------------------------------------  */

intgrt_vegas_t *integrate_vegas_cp(
    intgrt_vegas_t *vegas);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_vegas_set_wucalls(
    intgrt_vegas_t *vegas,
    size_t warmUpCalls);

int integrate_vegas_set_mcalls(
    intgrt_vegas_t *vegas,
    size_t mainCalls);

int integrate_vegas_set_ucalls(
    intgrt_vegas_t *vegas,
    double updateCalls);

int integrate_vegas_set_mloop(
    intgrt_vegas_t *vegas,
    int maxLoop);

int integrate_vegas_set_chisqerr(
    intgrt_vegas_t *vegas,
    double chisqErr);

int integrate_vegas_set_relerr(
    intgrt_vegas_t *vegas,
    double relErr);

int integrate_vegas_set_verbose(
    intgrt_vegas_t *vegas,
    int verbose);



/*  ----------------------------------------------------  */
/*  -----------------   Divonne Struct   ---------------  */
/*  ----------------------------------------------------  */


intgrt_divonne_t *integrate_divonne_new(void);

intgrt_divonne_t *integrate_divonne_free(
    intgrt_divonne_t *divonne);

/*  ----------------------------------------------------  */

intgrt_divonne_t *integrate_divonne_cp(
    intgrt_divonne_t *divonne);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_divonne_set_ncomp(
    intgrt_divonne_t *divonne,
    int nComp);

int integrate_divonne_set_nvec(
    intgrt_divonne_t *divonne,
    int nVec);


int integrate_divonne_set_epsrel(
    intgrt_divonne_t *divonne,
    double epsRel);

int integrate_divonne_set_epsabs(
    intgrt_divonne_t *divonne,
    double epsAbs);


int integrate_divonne_set_seed(
    intgrt_divonne_t *divonne,
    int seed);


int integrate_divonne_set_mineval(
    intgrt_divonne_t *divonne,
    int minEval);

int integrate_divonne_set_maxeval(
    intgrt_divonne_t *divonne,
    int maxEval);


int integrate_divonne_set_key1(
    intgrt_divonne_t *divonne,
    int key1);

int integrate_divonne_set_key2(
    intgrt_divonne_t *divonne,
    int key2);

int integrate_divonne_set_key3(
    intgrt_divonne_t *divonne,
    int key3);


int integrate_divonne_set_maxpass(
    intgrt_divonne_t *divonne,
    int maxPass);


int integrate_divonne_set_border(
    intgrt_divonne_t *divonne,
    double border);


int integrate_divonne_set_maxchisq(
    intgrt_divonne_t *divonne,
    double maxChisq);

int integrate_divonne_set_mindeviation(
    intgrt_divonne_t *divonne,
    double minDeviation);


int integrate_divonne_set_verbose(
    intgrt_divonne_t *divonne,
    int verbose);



/*  ----------------------------------------------------  */
/*  ----------------   Integrate Struct   --------------  */
/*  ----------------------------------------------------  */


intgrt_t *integrate_new(void);

intgrt_t *integrate_free(
    intgrt_t *intgrt);

/*  ----------------------------------------------------  */

intgrt_t *integrate_cp(
    intgrt_t *intgrt);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_set_bounds(
    intgrt_t *intgrt,
    size_t dim,
    double *ubounds,
    double *lbounds);


int integrate_set_dim(
    intgrt_t *intgrt,
    size_t dim);

int integrate_set_bounds_upper(
    intgrt_t *intgrt,
    double *ubounds);

int integrate_set_bounds_lower(
    intgrt_t *intgrt,
    double *lbounds);


int integrate_set_cquad(
    intgrt_t *intgrt,
    intgrt_cquad_t *cquad);

int integrate_set_vegas(
    intgrt_t *intgrt,
    intgrt_vegas_t *vegas);

int integrate_set_divonne(
    intgrt_t *intgrt,
    intgrt_divonne_t *divonne);


int integrate_set_routine(
    intgrt_t *intgrt,
    const char *routine);


int integrate_set_integrand(
    intgrt_t *intgrt,
    double (*integrand)(double*, size_t, void*));

int integrate_set_params(
    intgrt_t *intgrt,
    void *params);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


intgrt_cquad_t *integrate_get_cquad(
    intgrt_t *intgrt);

intgrt_vegas_t *integrate_get_vegas(
    intgrt_t *intgrt);

intgrt_divonne_t *integrate_get_divonne(
    intgrt_t *intgrt);




/*  ----------------------------------------------------  */
/*  --------------   Integrate Function   --------------  */
/*  ----------------------------------------------------  */


int integrate(
    double integrand(double*, size_t, void*),
    intgrt_t *intgrt,
    double *result);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_cquad(
    double integrand(double*, size_t, void*),
    intgrt_t *intgrt,
    double *result);

/*  ----------------------------------------------------  */

double integrate_cquad_integrand(
    double x,
    void *params);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_plain(
    double integrand(double*, size_t, void*),
    intgrt_t *intgrt,
    double *result);


int integrate_vegas(
    double integrand(double*, size_t, void*),
    intgrt_t *intgrt,
    double *result);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int integrate_divonne(
    double integrand(double*, size_t, void*),
    intgrt_t *intgrt,
    double *result);

/*  ----------------------------------------------------  */

int integrate_divonne_integrand(
    const int *nDim,
    const cubareal xx[],
    const int *nComp,
    cubareal ff[],
    void *usrData);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // INTEGRATE_H_INCLUDED

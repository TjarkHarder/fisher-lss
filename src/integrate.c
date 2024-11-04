/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTEGRATE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "integrate.h"


/**  Integration Routines  **/

const char *integrate_find_routine(const char *routine)
{
    /*

        Find the integration routine

    */

    /* CQUAD */
    if (misc_sinci(routine, _idIntgrtCQUAD_))
        return _idIntgrtCQUAD_;

    /* Vegas */
    if (misc_sinci(routine, _idIntgrtVegas_))
        return _idIntgrtVegas_;

    /* Divonne */
    if (misc_sinci(routine, _idIntgrtDivonne_))
        return _idIntgrtDivonne_;

    /* No routine found */
    printf("Did not find any '%s' integration routine.\n", routine);
    exit(1);

    return NULL;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTEGRATE STRUCTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   CQUAD Struct   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


intgrt_cquad_t *integrate_cquad_new(void)
{
    /*

        Create a new intgrt_cquad_t struct

    */

    intgrt_cquad_t *cquad = malloc(sizeof(intgrt_cquad_t));

    cquad -> absErr = 0.;
    cquad -> relErr = 0.;

    cquad -> limit = 0;

    return cquad;
}


intgrt_cquad_t *integrate_cquad_free(intgrt_cquad_t *cquad)
{
    /*

        Free cquad

    */

    /* Free cquad (no need to check if cquad is NULL) */
    free(cquad);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


intgrt_cquad_t *integrate_cquad_cp(intgrt_cquad_t *cquad)
{
    /*

        Copy cquad

    */

    /* Check for NULL */
    if (cquad == NULL)
        return NULL;

    intgrt_cquad_t *cquadCp = malloc(sizeof(intgrt_cquad_t));

    cquadCp -> absErr = cquad -> absErr;
    cquadCp -> relErr = cquad -> relErr;

    cquadCp -> limit = cquad -> limit;

    return cquadCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_cquad_set_abserr(intgrt_cquad_t *cquad, double absErr)
{
    /*

        Set the absolute error for cquad

    */

    cquad -> absErr = absErr;

    return 0;
}


int integrate_cquad_set_relerr(intgrt_cquad_t *cquad, double relErr)
{
    /*

        Set the relative error for cquad

    */

    cquad -> relErr = relErr;

    return 0;
}


int integrate_cquad_set_limit(intgrt_cquad_t *cquad, size_t limit)
{
    /*

        Set the maximum number of subintervals for cquad

    */

    cquad -> limit = limit;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Vegas Struct   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


intgrt_vegas_t *integrate_vegas_new(void)
{
    /*

        Create a new intgrt_vegas_t struct

    */

    intgrt_vegas_t *vegas = malloc(sizeof(intgrt_vegas_t));

    vegas -> warmUpCalls = 0;
    vegas -> mainCalls = 0;
    vegas -> updateCalls = 1.;

    vegas -> maxLoop = 0;

    vegas -> chisqErr = 0.;
    vegas -> relErr = 0.;

    vegas -> verbose = 0;

    return vegas;
}


intgrt_vegas_t *integrate_vegas_free(intgrt_vegas_t *vegas)
{
    /*

        Free vegas

    */

    /* Free vegas (no need to check if vegas is NULL) */
    free(vegas);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


intgrt_vegas_t *integrate_vegas_cp(intgrt_vegas_t *vegas)
{
    /*

        Copy vegas

    */

    /* Check for NULL */
    if (vegas == NULL)
        return NULL;

    intgrt_vegas_t *vegasCp = malloc(sizeof(intgrt_vegas_t));

    vegasCp -> warmUpCalls = vegas -> warmUpCalls;
    vegasCp -> mainCalls = vegas -> mainCalls;
    vegasCp -> updateCalls = vegas -> updateCalls;

    vegasCp -> maxLoop = vegas -> maxLoop;

    vegasCp -> chisqErr = vegas -> chisqErr;
    vegasCp -> relErr = vegas -> relErr;

    vegasCp -> verbose = vegas -> verbose;

    return vegasCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_vegas_set_wucalls(intgrt_vegas_t *vegas, size_t warmUpCalls)
{
    /*

        Set the warm up calls for vegas

    */

    vegas -> warmUpCalls = warmUpCalls;

    return 0;
}


int integrate_vegas_set_mcalls(intgrt_vegas_t *vegas, size_t mainCalls)
{
    /*

        Set the main calls for vegas

    */

    vegas -> mainCalls = mainCalls;

    return 0;
}


int integrate_vegas_set_ucalls(intgrt_vegas_t *vegas, double updateCalls)
{
    /*

        Set the factor by which the main calls are increased for vegas

    */

    vegas -> updateCalls = updateCalls;

    return 0;
}


int integrate_vegas_set_mloop(intgrt_vegas_t *vegas, int maxLoop)
{
    /*

        Set the max loop for vegas

    */

    vegas -> maxLoop = maxLoop;

    return 0;
}


int integrate_vegas_set_chisqerr(intgrt_vegas_t *vegas, double chisqErr)
{
    /*

        Set the chi squared error for vegas

    */

    vegas -> chisqErr = chisqErr;

    return 0;
}


int integrate_vegas_set_relerr(intgrt_vegas_t *vegas, double relErr)
{
    /*

        Set the min relative error for vegas

    */

    vegas -> relErr = relErr;

    return 0;
}


int integrate_vegas_set_verbose(intgrt_vegas_t *vegas, int verbose)
{
    /*

        Set the verbose level of output for vegas (ie 0 for no output and 1 if integral converged)

    */

    vegas -> verbose = verbose;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------------   Divonne Struct   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


intgrt_divonne_t *integrate_divonne_new(void)
{
    /*

        Create a new intgrt_divonne_t struct

    */

    intgrt_divonne_t *divonne = malloc(sizeof(intgrt_divonne_t));

    divonne -> nComp = 1;
    divonne -> nVec = 1;

    divonne -> epsRel = 1.e-3;
    divonne -> epsAbs = 1.e-12;

    divonne -> seed = 0;

    divonne -> minEval = 0;
    divonne -> maxEval = 50000;

    divonne -> key1 = 47;
    divonne -> key2 = 1;
    divonne -> key3 = 1;

    divonne -> maxPass = 5;

    divonne -> border = 0.;

    divonne -> maxChisq = 10.;
    divonne -> minDeviation = 0.25;

    divonne -> verbose = 0;

    return divonne;
}


intgrt_divonne_t *integrate_divonne_free(intgrt_divonne_t *divonne)
{
    /*

        Free divonne

    */

    /* Free divonne */
    free(divonne);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


intgrt_divonne_t *integrate_divonne_cp(intgrt_divonne_t *divonne)
{
    /*

        Copy divonne

    */

    /* Check for NULL */
    if (divonne == NULL)
        return NULL;

    intgrt_divonne_t *divonneCp = malloc(sizeof(intgrt_divonne_t));

    divonneCp -> nComp = divonne -> nComp;
    divonneCp -> nVec = divonne -> nVec;

    divonneCp -> epsRel = divonne -> epsRel;
    divonneCp -> epsAbs = divonne -> epsAbs;

    divonneCp -> seed = divonne -> seed;

    divonneCp -> minEval = divonne -> minEval;
    divonneCp -> maxEval = divonne -> maxEval;

    divonneCp -> key1 = divonne -> key1;
    divonneCp -> key2 = divonne -> key2;
    divonneCp -> key3 = divonne -> key3;

    divonneCp -> maxPass = divonne -> maxPass;

    divonneCp -> border = divonne -> border;

    divonneCp -> maxChisq = divonne -> maxChisq;
    divonneCp -> minDeviation = divonne -> minDeviation;

    divonneCp -> verbose = divonne -> verbose;

    return divonneCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_divonne_set_ncomp(intgrt_divonne_t *divonne, int nComp)
{
    /*

        Set ncomp in divonne

    */

    divonne -> nComp = nComp;

    return 0;
}


int integrate_divonne_set_nvec(intgrt_divonne_t *divonne, int nVec)
{
    /*

        Set nVec in divonne

    */

    divonne -> nVec = nVec;

    return 0;
}


int integrate_divonne_set_epsrel(intgrt_divonne_t *divonne, double epsRel)
{
    /*

        Set epsRel in divonne

    */

    divonne -> epsRel = epsRel;

    return 0;
}


int integrate_divonne_set_epsabs(intgrt_divonne_t *divonne, double epsAbs)
{
    /*

        Set epsAbs in divonne

    */

    divonne -> epsAbs = epsAbs;

    return 0;
}


int integrate_divonne_set_seed(intgrt_divonne_t *divonne, int seed)
{
    /*

        Set seed in divonne

    */

    divonne -> seed = seed;

    return 0;
}


int integrate_divonne_set_mineval(intgrt_divonne_t *divonne, int minEval)
{
    /*

        Set minEval in divonne

    */

    divonne -> minEval = minEval;

    return 0;
}


int integrate_divonne_set_maxeval(intgrt_divonne_t *divonne, int maxEval)
{
    /*

        Set maxEval in divonne

    */

    divonne -> maxEval = maxEval;

    return 0;
}


int integrate_divonne_set_key1(intgrt_divonne_t *divonne, int key1)
{
    /*

        Set key1 in divonne

    */

    divonne -> key1 = key1;

    return 0;
}


int integrate_divonne_set_key2(intgrt_divonne_t *divonne, int key2)
{
    /*

        Set key2 in divonne

    */

    divonne -> key2 = key2;

    return 0;
}


int integrate_divonne_set_key3(intgrt_divonne_t *divonne, int key3)
{
    /*

        Set key3 in divonne

    */

    divonne -> key3 = key3;

    return 0;
}


int integrate_divonne_set_maxpass(intgrt_divonne_t *divonne, int maxPass)
{
    /*

        Set maxPass in divonne

    */

    divonne -> maxPass = maxPass;

    return 0;
}


int integrate_divonne_set_border(intgrt_divonne_t *divonne, double border)
{
    /*

        Set border in divonne

    */

    divonne -> border = border;

    return 0;
}


int integrate_divonne_set_maxchisq(intgrt_divonne_t *divonne, double maxChisq)
{
    /*

        Set maxChisq in divonne

    */

    divonne -> maxChisq = maxChisq;

    return 0;
}


int integrate_divonne_set_mindeviation(intgrt_divonne_t *divonne, double minDeviation)
{
    /*

        Set minDeviation in divonne

    */

    divonne -> minDeviation = minDeviation;

    return 0;
}


int integrate_divonne_set_verbose(intgrt_divonne_t *divonne, int verbose)
{
    /*

        Set verbose in divonne

    */

    divonne -> verbose = verbose;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Integrate Struct   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


intgrt_t *integrate_new(void)
{
    /*

        Create a new intgrt_t struct

    */

    intgrt_t *intgrt = malloc(sizeof(intgrt_t));

    intgrt -> dim = 0;

    intgrt -> upperBounds = NULL;
    intgrt -> lowerBounds = NULL;

    intgrt -> cquad = NULL;
    intgrt -> vegas = NULL;
    intgrt -> divonne = NULL;

    intgrt -> routine = _idIntgrtVegas_;

    intgrt -> params = NULL;

    return intgrt;
}


intgrt_t *integrate_free(intgrt_t *intgrt)
{
    /*

        Free intgrt

    */

    /* Check for NULL */
    if (intgrt == NULL)
        return NULL;

    /* Free bounds */
    free(intgrt -> upperBounds);
    free(intgrt -> lowerBounds);

    /* Free cquad */
    intgrt -> cquad = integrate_cquad_free(intgrt -> cquad);

    /* Free vegas */
    intgrt -> vegas = integrate_vegas_free(intgrt -> vegas);

    /* Free divonne */
    intgrt -> divonne = integrate_divonne_free(intgrt -> divonne);

    /* Free intgrt itself */
    free(intgrt);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


intgrt_t *integrate_cp(intgrt_t *intgrt)
{
    /*

        Copy intgrt

    */

    /* Check for NULL */
    if (intgrt == NULL)
        return NULL;

    intgrt_t *intgrtCp = integrate_new();

    integrate_set_dim(intgrtCp, intgrt -> dim);
    integrate_set_bounds_upper(intgrtCp, intgrt -> upperBounds);
    integrate_set_bounds_lower(intgrtCp, intgrt -> lowerBounds);

    integrate_set_cquad(intgrtCp, intgrt -> cquad);
    integrate_set_vegas(intgrtCp, intgrt -> vegas);
    integrate_set_divonne(intgrtCp, intgrt -> divonne);

    integrate_set_routine(intgrtCp, intgrt -> routine);

    integrate_set_params(intgrtCp, intgrt -> params);

    return intgrtCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_set_bounds(intgrt_t *intgrt, size_t dim, double *ubounds, double *lbounds)
{
    /*

        Set the integration bounds for intgrt

    */

    intgrt -> dim = dim;

    intgrt -> upperBounds = realloc(intgrt -> upperBounds, sizeof(double) * dim);
    intgrt -> lowerBounds = realloc(intgrt -> lowerBounds, sizeof(double) * dim);

    for (size_t i = 0; i < dim; i++)
      {
        intgrt -> upperBounds[i] = ubounds[i];
        intgrt -> lowerBounds[i] = lbounds[i];
      }

    return 0;
}


int integrate_set_dim(intgrt_t *intgrt, size_t dim)
{
    /*

        Set the integration dimension for intgrt

    */

    intgrt -> dim = dim;

    intgrt -> upperBounds = realloc(intgrt -> upperBounds, sizeof(double) * dim);
    intgrt -> lowerBounds = realloc(intgrt -> lowerBounds, sizeof(double) * dim);

    return 0;
}


int integrate_set_bounds_upper(intgrt_t *intgrt, double *ubounds)
{
    /*

        Set the integration upper bounds for intgrt

    */

    if (intgrt -> dim == 0)
      {
        printf("Cannot set the integration upper bounds if the dimension is not yet specified.\n");
        exit(1);

        return 1;
      }

    for (size_t i = 0; i < intgrt -> dim; i++)
      {
        intgrt -> upperBounds[i] = ubounds[i];
      }

    return 0;
}


int integrate_set_bounds_lower(intgrt_t *intgrt, double *lbounds)
{
    /*

        Set the integration lower bounds for intgrt

    */

    if (intgrt -> dim == 0)
      {
        printf("Cannot set the integration lower bounds if the dimension is not yet specified.\n");
        exit(1);

        return 1;
      }

    for (size_t i = 0; i < intgrt -> dim; i++)
      {
        intgrt -> lowerBounds[i] = lbounds[i];
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int integrate_set_routine(intgrt_t *intgrt, const char *routine)
{
    /*

        Set the integration routine

    */

    intgrt -> routine = integrate_find_routine(routine);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int integrate_set_cquad(intgrt_t *intgrt, intgrt_cquad_t *cquad)
{
    /*

        Set the cquad parameters for intgrt

    */

    intgrt -> cquad = integrate_cquad_free(intgrt -> cquad);
    intgrt -> cquad = integrate_cquad_cp(cquad);

    return 0;
}


int integrate_set_vegas(intgrt_t *intgrt, intgrt_vegas_t *vegas)
{
    /*

        Set the vegas parameters for intgrt

    */

    intgrt -> vegas = integrate_vegas_free(intgrt -> vegas);
    intgrt -> vegas = integrate_vegas_cp(vegas);

    return 0;
}


int integrate_set_divonne(intgrt_t *intgrt, intgrt_divonne_t *divonne)
{
    /*

        Set the divonne parameters for intgrt

    */

    intgrt -> divonne = integrate_divonne_free(intgrt -> divonne);
    intgrt -> divonne = integrate_divonne_cp(divonne);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int integrate_set_integrand(intgrt_t *intgrt, double (*integrand)(double*, size_t, void*))
{
    /*

        Set the integrand for intgrt

    */

    intgrt -> integrand = integrand;

    return 0;
}


int integrate_set_params(intgrt_t *intgrt, void *params)
{
    /*

        Set the additional parameters of the integrand for intgrt

    */

    intgrt -> params = params;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


intgrt_cquad_t *integrate_get_cquad(intgrt_t *intgrt)
{
    /*

        Get the cquad parameters from intgrt

    */

    return intgrt -> cquad;
}


intgrt_vegas_t *integrate_get_vegas(intgrt_t *intgrt)
{
    /*

        Get the vegas parameters from intgrt

    */

    return intgrt -> vegas;
}


intgrt_divonne_t *integrate_get_divonne(intgrt_t *intgrt)
{
    /*

        Get the divonne parameters from intgrt

    */

    return intgrt -> divonne;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTEGRATE FUNCTION     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------   General Integrate Function   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate(double integrand(double*, size_t, void*), intgrt_t *intgrt, double *result)
{
    /*

        Integrate a function over a region (provided in intgrt) using the integration routine provided in intgrt -> routine.

    */

    /* CQUAD */
    if (misc_sinci((const char*) intgrt -> routine, _idIntgrtCQUAD_))
      {
        integrate_cquad(integrand, intgrt, result);

        return 0;
      }

    /* Vegas */
    if (misc_sinci((const char*) intgrt -> routine, _idIntgrtVegas_))
      {
        integrate_vegas(integrand, intgrt, result);

        return 0;
      }

    /* Divonne */
    if (misc_sinci((const char*) intgrt -> routine, _idIntgrtDivonne_))
      {
        integrate_divonne(integrand, intgrt, result);

        return 0;
      }

    /* Did not find routine */
    printf("Did not find any '%s' integration routine.\n", intgrt -> routine);
    exit(1);

    return 1;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   GSL's CQUAD   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_cquad(double integrand(double*, size_t, void*), intgrt_t *intgrt, double *result)
{
    /*

        Integrate a function over a 1d region (provided in intgrt) using the CQUAD routine from GSL.

    */

    /* Reset result */
    result[0] = 0.;
    result[1] = 0.;

    /* Set integrand */
    integrate_set_integrand(intgrt, integrand);

    /* intgrt_cquad_t struct */
    intgrt_cquad_t *cquad = intgrt -> cquad;

    /* Allocate working memory for CQUAD */
    gsl_integration_cquad_workspace *work = gsl_integration_cquad_workspace_alloc(cquad -> limit);

    /* Integrand to be passed to CQUAD */
    gsl_function integrandCQUAD;
    integrandCQUAD.function = &integrate_cquad_integrand;
    integrandCQUAD.params = intgrt;

    /* Integrate */
    gsl_integration_cquad(&integrandCQUAD, *(intgrt -> lowerBounds), *(intgrt -> upperBounds), cquad -> absErr, cquad -> relErr, work, &result[0], &result[1], NULL);

    /* Free allocated memory */
    gsl_integration_cquad_workspace_free(work);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double integrate_cquad_integrand(double x, void *params)
{
    /*

        Integrand for integrate_divonne

    */

    intgrt_t *intgrt = (intgrt_t*) params;

    return intgrt -> integrand(&x, intgrt -> dim, intgrt -> params);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   GSL's VEGAS   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_vegas(double integrand(double*, size_t, void*), intgrt_t *intgrt, double *result)
{
    /*

        Integrate a function over a region (provided in intgrt) using the Monte-Carlo integration algorithm VEGAS from GSL.

    */

    /* Reset result */
    result[0] = 0.;
    result[1] = 0.;
    result[2] = 0.;

    /* intgrt_vegas_t struct */
    intgrt_vegas_t *vegas = intgrt -> vegas;

    /* Setup VEGAS */

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function monteIntegral = {integrand, intgrt -> dim, intgrt -> params};

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(intgrt -> dim);


    /* VEGAS warm up */

    /* Integrate */
    gsl_monte_vegas_integrate(&monteIntegral, intgrt -> lowerBounds, intgrt -> upperBounds, intgrt -> dim, vegas -> warmUpCalls, r, s, &result[0], &result[1]);

    /* χ^2 */
    result[2] = gsl_monte_vegas_chisq(s);


    /* VEGAS main stage */

    int loop = 0;
    size_t mainCalls = vegas -> mainCalls;
    double relErr = 0.;


    do
      {
        /* If the algorithm does not converge break the loop after vegas -> maxLoop iterations */
        if (loop == vegas -> maxLoop)
          {
            if (vegas -> verbose)
                printf("Integral did not converge after %d iterations.\n", vegas -> maxLoop);

            break;
          }

        loop++;

        /* Integrate */
        gsl_monte_vegas_integrate(&monteIntegral, intgrt -> lowerBounds, intgrt -> upperBounds, intgrt -> dim, mainCalls, r, s, &result[0], &result[1]);

        /* χ^2 */
        result[2] = gsl_monte_vegas_chisq(s);

        /* Check for relative error */
        if (vegas -> relErr != 0.)
          {
            /* Relative error */
            relErr = (result[0] != 0.) ? result[1] / fabs(result[0]) : 0.;

            /* Relative error out of bounds -> increase number of function evaluations and continue loop */
            if (relErr > vegas -> relErr)
              {
                mainCalls = (size_t) (((double) mainCalls) * vegas -> updateCalls);

                continue;
              }

            /* Relative error in bounds and chisq error in bounds -> break loop */
            else if (fabs(result[2] - 1.) < vegas -> chisqErr)
              {
                if (vegas -> verbose)
                    printf("Integral converged after %d iterations.\n", loop);

                break;
              }
          }

        /* Chisq error in bounds -> break loop */
        if (fabs(result[2] - 1.) < vegas -> chisqErr)
          {
            if (vegas -> verbose)
                printf("Integral converged after %d iterations.\n", loop);

            break;
          }
      }

    while (true);


    /* Free allocated memory */

    gsl_rng_free(r);
    gsl_monte_vegas_free(s);

    return 0;
}


int integrate_plain(double integrand(double*, size_t, void*), intgrt_t *intgrt, double *result)
{
    /*

        Integrate a function over a region (provided in "intgrt") using the plain Monte-Carlo integration from GSL.

    */

    /* Reset result */
    result[0] = 0.;
    result[1] = 0.;
    result[2] = 0.;

    /* Setup */
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function monteIntegral = {integrand, intgrt -> dim, intgrt -> params};

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_plain_state *s = gsl_monte_plain_alloc(intgrt -> dim);

    /* Integrate */
    gsl_monte_plain_integrate(&monteIntegral, intgrt -> lowerBounds, intgrt -> upperBounds, intgrt -> dim, intgrt -> vegas -> mainCalls, r, s, &result[0], &result[1]);

    /* Free allocated memory */
    gsl_rng_free(r);
    gsl_monte_plain_free(s);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   CUBA's Divonne   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_divonne(double integrand(double*, size_t, void*), intgrt_t *intgrt, double *result)
{
    /*

        Integrate a function over a region (provided in "intgrt") with the Divonne integration from CUBA

    */

    int nRegions, nEval, fail;
    cubareal integral, error, prob;

    /* Divonne struct */
    intgrt_divonne_t *divonne = intgrt -> divonne;

    /* Set integrand */
    integrate_set_integrand(intgrt, integrand);

    /* Integrate */
    Divonne((int) intgrt -> dim, divonne -> nComp, integrate_divonne_integrand, intgrt, divonne -> nVec,
            divonne -> epsRel, divonne -> epsAbs,
            divonne -> verbose,
            divonne -> seed,
            divonne -> minEval, divonne -> maxEval,
            divonne -> key1, divonne -> key2, divonne -> key3,
            divonne -> maxPass,
            divonne -> border,
            divonne -> maxChisq, divonne -> minDeviation,
            0, 0, NULL, 0, NULL, NULL, NULL, // Never used
            &nRegions, &nEval, &fail, &integral, &error, &prob); // Output

    /* Result */
    result[0] = (double) integral;
    result[1] = (double) error;
    result[2] = (double) prob;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int integrate_divonne_integrand(const int *nDim, const cubareal xx[], const int *nComp, cubareal ff[], void *usrData)
{
    /*

        Integrand for integrate_divonne

    */

    (void) nComp;

    /* Integrate struct */
    intgrt_t *intgrt = (intgrt_t*) usrData;

    /* x Data */
    double *xData = malloc(sizeof(double) * ((size_t) *nDim));

    /* Jacobian to obtain unicube as integration region */
    double jac = 1.;

    /* Transform to unicube */
    for (size_t i = 0; i < (size_t) *nDim; i++)
      {
        xData[i] = intgrt -> lowerBounds[i] + xx[i] * (intgrt -> upperBounds[i] - intgrt -> lowerBounds[i]);
        jac *= intgrt -> upperBounds[i] - intgrt -> lowerBounds[i];
      }

    /* Integrand */
    ff[0] = intgrt -> integrand(xData, (size_t) *nDim, intgrt -> params) * jac;

    /* Free memory */
    free(xData);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     AVERAGE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "average.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE / SETUP / FREE VARIABLES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Volume Average  **/

/* Line */
static int (*_avrLineVolFunc)(spec_arg_t*, double*) = NULL;

static const size_t _avrLineVolOrder = 2;

static intgrt_t *_avrLineVolIntgrt = NULL;

static const size_t _avrLineVolIntgrtDim = 2;
static const char *_avrLineVolIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrLineVolIntgrtVegas = {1000, 100, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrLineVolIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};

/* Triangle */
static int (*_avrTriVolFunc)(spec_arg_t*, double*) = NULL;

static const size_t _avrTriVolOrder = 3;

static intgrt_t *_avrTriVolIntgrt = NULL;

static const size_t _avrTriVolIntgrtDim = 5;
static const char *_avrTriVolIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrTriVolIntgrtVegas = {10000, 1000, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrTriVolIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};


/**  Shape Average  **/

/* Line */
static const size_t _avrLineOrder = 2;

static intgrt_t *_avrLineIntgrt = NULL;

static const size_t _avrLineIntgrtDim = 2;
static const char *_avrLineIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrLineIntgrtVegas = {1000, 100, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrLineIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};

/* Triangle */
static const size_t _avrTriOrder = 3;

static intgrt_t *_avrTriIntgrt = NULL;

static const size_t _avrTriIntgrtDim = 5;
static const char *_avrTriIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrTriIntgrtVegas = {10000, 1000, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrTriIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};


/**  PP Covarianece Matrix Average  **/

/* Gaussian */
static intgrt_t *_avrCovPPGaussIntgrt = NULL;

static const size_t _avrCovPPGaussIntgrtDim = 2;
static const char *_avrCovPPGaussIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrCovPPGaussIntgrtVegas = {1000, 100, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrCovPPGaussIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};


/* Non-Gaussian (Infinitesimal limit) */
static intgrt_t *_avrCovPPNGaussInfIntgrt = NULL;

static const size_t _avrCovPPNGaussInfIntgrtDim = 1;
static const char *_avrCovPPNGaussInfIntgrtRoutine = "cquad";
static const intgrt_cquad_t _avrCovPPNGaussInfIntgrtCQUAD = {0., 1e-4, 1000};
static const intgrt_vegas_t _avrCovPPNGaussInfIntgrtVegas = {1000, 100, 1., 10, 0.5, 0., 0};

/* Non-Gaussian */
static intgrt_t *_avrCovPPNGaussIntgrt = NULL;

static const size_t _avrCovPPNGaussIntgrtDim = 5;
static const char *_avrCovPPNGaussIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrCovPPNGaussIntgrtVegas = {1000, 100, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrCovPPNGaussIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};


/**  BB Covarianece Matrix Average  **/

/* Gaussian */
static intgrt_t *_avrCovBBGaussIntgrt = NULL;

static const size_t _avrCovBBGaussIntgrtDim = 5;
static const char *_avrCovBBGaussIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrCovBBGaussIntgrtVegas = {100000, 10000, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _avrCovBBGaussIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};

/* Non-Gaussian */
static intgrt_t *_avrCovBBNGaussIntgrt = NULL;

static const size_t _avrCovBBNGaussIntgrtDim = 8;
static const char *_avrCovBBNGaussIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrCovBBNGaussIntgrtVegas = {100000, 10000, 1., 10, 0.5, 0., 1};
static const intgrt_divonne_t _avrCovBBNGaussIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};


/**  PB Covarianece Matrix Average  **/

/* Gaussian (does not exist) */

/* Non-Gaussian */
static intgrt_t *_avrCovPBNGaussIntgrt = NULL;

static const size_t _avrCovPBNGaussIntgrtDim = 5;
static const char *_avrCovPBNGaussIntgrtRoutine = "vegas";
static const intgrt_vegas_t _avrCovPBNGaussIntgrtVegas = {10000, 1000, 1., 10, 0.5, 0., 1};
static const intgrt_divonne_t _avrCovPBNGaussIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Setters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_set_shape_vol_func(const char *idShape, const char *idFunc)
{
    /*

        Set the label for the shape volume function

    */

    /* Line */
    if (!strcmp(idShape, _idAvrLineVol_))
      {
        if (!strcmp(idFunc, _idAvrLineVolFuncNum_))
          {
            _avrLineVolFunc = avr_shape_line_vol;
            return 0;
          }

        if (!strcmp(idFunc, _idAvrLineVolFuncAna_))
          {
            _avrLineVolFunc = avr_shape_line_vol_ana;
            return 0;
          }

        printf("Line-Bin-Volume function with id '%s' does not exist.\n", idFunc);
        exit(1);

        return 1;
      }

    /* Triangle */
    if (!strcmp(idShape, _idAvrTriVol_))
      {
        if (!strcmp(idFunc, _idAvrTriVolFuncNum_))
          {
            _avrTriVolFunc = avr_shape_tri_vol;
            return 0;
          }

        if (!strcmp(idFunc, _idAvrTriVolFuncAna_))
          {
            _avrTriVolFunc = avr_shape_tri_vol_ana;
            return 0;
          }

        printf("Triangle-Bin-Volume function with id '%s' does not exist.\n", idFunc);
        exit(1);

        return 1;
      }

    printf("Bin-Volume function for shape with id '%s' does not exist.\n", idShape);
    exit(1);

    return 1;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Getters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


intgrt_t *avr_get_integrate(const char *id)
{
    /*

        Get the intgrt_vegas_t struct for a given id

    */

    /* Line Volume */
    if (!strcmp(id, _idAvrLineVol_))
        return _avrLineVolIntgrt;

    /* Triangle Volume */
    if (!strcmp(id, _idAvrTriVol_))
        return _avrTriVolIntgrt;


    /* Line Average */
    if (!strcmp(id, _idAvrLine_))
        return _avrLineIntgrt;

    /* Triangle Average */
    if (!strcmp(id, _idAvrTri_))
        return _avrTriIntgrt;


    /* Gaussian Cov PP Average */
    if (!strcmp(id, _idAvrCovPPGauss_))
        return _avrCovPPGaussIntgrt;

    /* Non-Gaussian Cov PP Average (Infinitesimal Limit) */
    if (!strcmp(id, _idAvrCovPPNGaussInf_))
        return _avrCovPPNGaussInfIntgrt;

    /* Non-Gaussian Cov PP Average */
    if (!strcmp(id, _idAvrCovPPNGauss_))
        return _avrCovPPNGaussIntgrt;


    /* Gaussian Cov BB Average */
    if (!strcmp(id, _idAvrCovBBGauss_))
        return _avrCovBBGaussIntgrt;

    /* Non-Gaussian Cov BB Average */
    if (!strcmp(id, _idAvrCovBBNGauss_))
        return _avrCovBBNGaussIntgrt;


    /* Non-Gaussian Cov PB Average */
    if (!strcmp(id, _idAvrCovPBNGauss_))
        return _avrCovPBNGaussIntgrt;


    printf("Cannot get the 'intgrt_t' struct for unknown id '%s'\n.", id);
    exit(1);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Initialise Local Variables   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_ini(void)
{
    /*


        Initialise local variables

    */

    /**  Shape Volume  **/

    /* Line */
    avr_set_shape_vol_func(_idAvrLineVol_, _idAvrLineVolFuncNum_);

    _avrLineVolIntgrt = integrate_new();

    integrate_set_dim(_avrLineVolIntgrt, (size_t) _avrLineVolIntgrtDim);

    integrate_set_routine(_avrLineVolIntgrt, (const char*) _avrLineVolIntgrtRoutine);

    integrate_set_vegas(_avrLineVolIntgrt, (intgrt_vegas_t*) &_avrLineVolIntgrtVegas);
    integrate_set_divonne(_avrLineVolIntgrt, (intgrt_divonne_t*) &_avrLineVolIntgrtDivonne);

    /* Triangle */
    avr_set_shape_vol_func(_idAvrTriVol_, _idAvrTriVolFuncNum_);

    _avrTriVolIntgrt = integrate_new();

    integrate_set_dim(_avrTriVolIntgrt, (size_t) _avrTriVolIntgrtDim);

    integrate_set_routine(_avrTriVolIntgrt, (const char*) _avrTriVolIntgrtRoutine);

    integrate_set_vegas(_avrTriVolIntgrt, (intgrt_vegas_t*) &_avrTriVolIntgrtVegas);
    integrate_set_divonne(_avrTriVolIntgrt, (intgrt_divonne_t*) &_avrTriVolIntgrtDivonne);


    /**  Shape Average  **/

    /* Line */
    _avrLineIntgrt = integrate_new();

    integrate_set_dim(_avrLineIntgrt, (size_t) _avrLineIntgrtDim);

    integrate_set_routine(_avrLineIntgrt, (const char*) _avrLineIntgrtRoutine);

    integrate_set_vegas(_avrLineIntgrt, (intgrt_vegas_t*) &_avrLineIntgrtVegas);
    integrate_set_divonne(_avrLineIntgrt, (intgrt_divonne_t*) &_avrLineIntgrtDivonne);

    /* Triangle */
    _avrTriIntgrt = integrate_new();

    integrate_set_dim(_avrTriIntgrt, (size_t) _avrTriIntgrtDim);

    integrate_set_routine(_avrTriIntgrt, (const char*) _avrTriIntgrtRoutine);

    integrate_set_vegas(_avrTriIntgrt, (intgrt_vegas_t*) &_avrTriIntgrtVegas);
    integrate_set_divonne(_avrTriIntgrt, (intgrt_divonne_t*) &_avrTriIntgrtDivonne);


    /**  PP Covariance Matrix Average  **/

    /* Gaussian */
    _avrCovPPGaussIntgrt = integrate_new();

    integrate_set_dim(_avrCovPPGaussIntgrt, (size_t) _avrCovPPGaussIntgrtDim);

    integrate_set_routine(_avrCovPPGaussIntgrt, (const char*) _avrCovPPGaussIntgrtRoutine);

    integrate_set_vegas(_avrCovPPGaussIntgrt, (intgrt_vegas_t*) &_avrCovPPGaussIntgrtVegas);
    integrate_set_divonne(_avrCovPPGaussIntgrt, (intgrt_divonne_t*) &_avrCovPPGaussIntgrtDivonne);

    /* Non-Gaussian (Infinitesimal limit) */
    _avrCovPPNGaussInfIntgrt = integrate_new();

    integrate_set_dim(_avrCovPPNGaussInfIntgrt, (size_t) _avrCovPPNGaussInfIntgrtDim);

    integrate_set_routine(_avrCovPPNGaussInfIntgrt, (const char*) _avrCovPPNGaussInfIntgrtRoutine);

    integrate_set_cquad(_avrCovPPNGaussInfIntgrt, (intgrt_cquad_t*) &_avrCovPPNGaussInfIntgrtCQUAD);
    integrate_set_vegas(_avrCovPPNGaussInfIntgrt, (intgrt_vegas_t*) &_avrCovPPNGaussInfIntgrtVegas);

    /* Non-Gaussian */
    _avrCovPPNGaussIntgrt = integrate_new();

    integrate_set_dim(_avrCovPPNGaussIntgrt, (size_t) _avrCovPPNGaussIntgrtDim);

    integrate_set_routine(_avrCovPPNGaussIntgrt, (const char*) _avrCovPPNGaussIntgrtRoutine);

    integrate_set_vegas(_avrCovPPNGaussIntgrt, (intgrt_vegas_t*) &_avrCovPPNGaussIntgrtVegas);
    integrate_set_divonne(_avrCovPPNGaussIntgrt, (intgrt_divonne_t*) &_avrCovPPNGaussIntgrtDivonne);


    /**  BB Covariance Matrix Average  **/

    /* Gaussian */
    _avrCovBBGaussIntgrt = integrate_new();

    integrate_set_dim(_avrCovBBGaussIntgrt, (size_t) _avrCovBBGaussIntgrtDim);

    integrate_set_routine(_avrCovBBGaussIntgrt, (const char*) _avrCovBBGaussIntgrtRoutine);

    integrate_set_vegas(_avrCovBBGaussIntgrt, (intgrt_vegas_t*) &_avrCovBBGaussIntgrtVegas);
    integrate_set_divonne(_avrCovBBGaussIntgrt, (intgrt_divonne_t*) &_avrCovBBGaussIntgrtDivonne);

    /* Non-Gaussian */
    _avrCovBBNGaussIntgrt = integrate_new();

    integrate_set_dim(_avrCovBBNGaussIntgrt, (size_t) _avrCovBBNGaussIntgrtDim);

    integrate_set_routine(_avrCovBBNGaussIntgrt, (const char*) _avrCovBBNGaussIntgrtRoutine);

    integrate_set_vegas(_avrCovBBNGaussIntgrt, (intgrt_vegas_t*) &_avrCovBBNGaussIntgrtVegas);
    integrate_set_divonne(_avrCovBBNGaussIntgrt, (intgrt_divonne_t*) &_avrCovBBNGaussIntgrtDivonne);


    /**  PB Covariance Matrix Average  **/

    /* Gaussian (does not exist) */

    /* Non-Gaussian */
    _avrCovPBNGaussIntgrt = integrate_new();

    integrate_set_dim(_avrCovPBNGaussIntgrt, (size_t) _avrCovPBNGaussIntgrtDim);

    integrate_set_routine(_avrCovPBNGaussIntgrt, (const char*) _avrCovPBNGaussIntgrtRoutine);

    integrate_set_vegas(_avrCovPBNGaussIntgrt, (intgrt_vegas_t*) &_avrCovPBNGaussIntgrtVegas);
    integrate_set_divonne(_avrCovPBNGaussIntgrt, (intgrt_divonne_t*) &_avrCovPBNGaussIntgrtDivonne);


    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Free Local Variables   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_free(void)
{
    /*

        Free local variables

    */


    /**  Shape Volume  **/

    /* Line */
    _avrLineVolIntgrt = integrate_free(_avrLineVolIntgrt);

    /* Triangle */
    _avrTriVolIntgrt = integrate_free(_avrTriVolIntgrt);


    /**  Shape Average  **/

    /* Line */
    _avrLineIntgrt = integrate_free(_avrLineIntgrt);

    /* Triangle */
    _avrTriIntgrt = integrate_free(_avrTriIntgrt);


    /**  PP Covariance Matrix Average  **/

    /* Gaussian */
    _avrCovPPGaussIntgrt = integrate_free(_avrCovPPGaussIntgrt);

    /* Non-Gaussian (Infinitesimal limit) */
    _avrCovPPNGaussInfIntgrt = integrate_free(_avrCovPPNGaussInfIntgrt);

    /* Non-Gaussian */
    _avrCovPPNGaussIntgrt = integrate_free(_avrCovPPNGaussIntgrt);


    /**  BB Covariance Matrix Average  **/

    /* Gaussian */
    _avrCovBBGaussIntgrt = integrate_free(_avrCovBBGaussIntgrt);

    /* Non-Gaussian */
    _avrCovBBNGaussIntgrt = integrate_free(_avrCovBBNGaussIntgrt);


    /**  BB Covariance Matrix Average  **/

    /* Gaussian (does not exist) */

    /* Non-Gaussian */
    _avrCovPBNGaussIntgrt = integrate_free(_avrCovPBNGaussIntgrt);


    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     AVERAGE SHAPES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Direct Functions   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double avr_shape_direct(int (*avr)(double (*)(void*, void*), spec_arg_t*, double*), double (*specFunc)(void*, void*), spec_arg_t *specArg)
{
    /*

        Return the result of avr directly

    */

    /* Get the result */
    double result[3];
    avr(specFunc, specArg, result);

    return result[0];
}


double avr_shape_vol_direct(int (*vol)(spec_arg_t*, double*), spec_arg_t *specArg)
{
    /*

        Return the result of vol directly

    */

    /* Get the result */
    double result[3];
    vol(specArg, result);

    return result[0];
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Line-Bin-Average   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_line(double (*specFunc)(void*, void*), spec_arg_t *specArg, double *result)
{
    /*

        Calculate the line-bin-average over a given function

    */

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    double dk = specArg -> dk;
    double dmu = specArg -> dmu;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrLineIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrLineIntgrt -> dim);

    upperBounds[0] = k + dk / 2.;
    upperBounds[1] = mu + dmu / 2.;

    lowerBounds[0] = k - dk / 2.;
    lowerBounds[1] = mu - dmu / 2.;

    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL; // TODO: Remove this?

    double bounds[4] = {lowerBounds[0], upperBounds[0], lowerBounds[1], upperBounds[1]};
    avrParams -> bounds = bounds;

    avrParams -> specFunc = specFunc;

    /* Integrate struct */
    intgrt_t *intgrt = (specFunc == _oneFunc_) ? integrate_cp(_avrLineVolIntgrt) : integrate_cp(_avrLineIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);

    /* Integrate */
    integrate(integrand_average_line, intgrt, result);

    result[0] /= pow(2. * M_PI, 3.) * specArg -> binVolume;
    result[1] /= pow(2. * M_PI, 3.) * specArg -> binVolume;


    /* Reset variables */
    kernels_qset_k(kern, 0, k);
    kernels_qset_k(kern, 1, k);

    kernels_qset_mu(kern, 0, mu);
    kernels_qset_mu(kern, 1, -mu);

    kernels_qset_nu(kern, 0, 1, -1.);

    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_line_vol(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the line-bin-volume

    */

    /* Get the bin volume */
    avr_shape_line(_oneFunc_, specArg, result);

    /* Must correct for (possibly incorrect) normalisation */
    result[0] *= specArg -> binVolume;
    result[1] *= specArg -> binVolume;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_line_vol_ana(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the line-bin-volume analytically

    */

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_qget_k(kern, 0);
//    double mu = kernels_qget_mu(kern, 0);

    double dk = specArg -> dk;
    double dmu = specArg -> dmu;

    /* Result */
    result[0] = 1. / pow(2. * M_PI, 2.) * (k*k * dk + k * dk*dk / 12.) * dmu;
    result[1] = 0.;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Triangle-Bin-Average   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_tri(double (*specFunc)(void*, void*), spec_arg_t *specArg, double *result)
{
    /*

        Calculate the triangle-bin-average over a given function

    */

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);
    double nu13 = kernels_qget_nu(kern, 0, 2);
    double nu23 = kernels_qget_nu(kern, 1, 2);

    double dk = specArg -> dk;
    double dmu = specArg -> dmu;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrTriIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrTriIntgrt -> dim);

    upperBounds[0] = k1 + dk / 2.;
    upperBounds[1] = k2 + dk / 2.;
    upperBounds[2] = mu1 + dmu / 2.;
    upperBounds[3] = mu2 + dmu / 2.;
    upperBounds[4] = 2. * M_PI;

    lowerBounds[0] = k1 - dk / 2.;
    lowerBounds[1] = k2 - dk / 2.;
    lowerBounds[2] = mu1 - dmu / 2.;
    lowerBounds[3] = mu2 - dmu / 2.;
    lowerBounds[4] = 0.;

    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL; // TODO: Remove this?

    double bounds[4] = {k3 - dk / 2., k3 + dk / 2., mu3 - dmu / 2., mu3 + dmu / 2.};
    avrParams -> bounds = bounds;

    avrParams -> specFunc = specFunc;

    /* Integrate struct */
    intgrt_t *intgrt = (specFunc == _oneFunc_) ? integrate_cp(_avrTriVolIntgrt) : integrate_cp(_avrTriIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);


    /* Integrate */
    integrate(integrand_average_tri, intgrt, result);

    result[0] /= pow(2. * M_PI, 6.) * specArg -> binVolume;
    result[1] /= pow(2. * M_PI, 6.) * specArg -> binVolume;


    /* Reset variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k3);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu3);

    kernels_qset_nu(kern, 0, 1, nu12);
    kernels_qset_nu(kern, 0, 2, nu13);
    kernels_qset_nu(kern, 1, 2, nu23);

    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_tri_vol(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the triangle-bin-volume

    */

    /* Get the bin volume */
    avr_shape_tri(_oneFunc_, specArg, result);

    /* Must correct for (possibly incorrect) normalisation */
    result[0] *= specArg -> binVolume;
    result[1] *= specArg -> binVolume;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_tri_vol_ana(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the triangle-bin-volume analytically (YP (2018))

    */

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);

    double nu12 = kernels_qget_nu(kern, 0, 1);

    double dk = specArg -> dk;
    double dmu = specArg -> dmu;

    /* Result */
    result[0] = 2. / pow(2. * M_PI, 5.) * k1*k2*k3 * dk*dk*dk * dmu*dmu / sqrt(1. - nu12*nu12 - mu1*mu1 - mu2*mu2 + 2.*mu1*mu2*nu12);
    result[1] = 0.;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   Get Average and Volume Functions   --------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*avr_shape_get_func(size_t order))(double(*)(void*, void*), spec_arg_t*, double*)
{
    /*

        Get the average function of a shape of given order

    */

    /* Line */
    if (order == _avrLineOrder)
      {
        return avr_shape_line;
      }

    /* Triangle */
    if (order == _avrTriOrder)
      {
        return avr_shape_tri;
      }

    printf("A shape-bin-average function of order '%ld' is not supported.\n", order);
    exit(1);

    return NULL;
}


int (*avr_shape_vol_get_func(size_t order))(spec_arg_t*, double*)
{
    /*

        Get the volume function of a shape of given order

    */

    /* Line */
    if (order == _avrLineVolOrder)
      {
        return _avrLineVolFunc;
      }

    /* Triangle */
    if (order == _avrTriVolOrder)
      {
        return _avrTriVolFunc;
      }

    printf("A shape-bin-volume function of order '%ld' is not supported.\n", order);
    exit(1);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------   Average over Infinitesimally Small Bins   -----------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_inf(double (*specFunc)(void*, void*), spec_arg_t *specArg, double *result)
{
    /*

        Bin-average for infinitesimally small bins simply returns the function itself

    */

    result[0] = specFunc(specArg, NULL);
    result[1] = 0.;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------   Shape Bin-Volume for Sample of Shapes   ------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_shape_sample_vol(sample_shape_t *sampleShape)
{
    /*

        Calculate the shape-bin-volumes for a sample of shapes

    */

    /* Average function */
    int (*avrFunc)(spec_arg_t*, double*) = avr_shape_vol_get_func(sampleShape -> dim);

    #pragma omp parallel
      {

        /* New spec_arg_t struct */
        spec_arg_t *specArg = spec_arg_new(NULL);

        /* New kern_t struct */
        kern_t *kern = kernels_new(sampleShape -> dim, 0);
        spec_arg_set_kern(specArg, kern);

        /* Set step sizes */
        spec_arg_set_dk(specArg, sampleShape -> sampleRawLength -> sampleArg -> step);
        spec_arg_set_dmu(specArg, sampleShape -> sampleRawOrientation -> sampleArg -> step);

        /* Calculate the bin-volumes */
        #pragma omp for
        for (size_t n = 0; n < sampleShape -> size; n++)
          {
            /* Get the n'th shape */
            shape_t *shape = sampleShape -> arrayShape[n];

            /* Set the variables */
            for (size_t i = 0; i < shape -> dim; i++)
              {
                kernels_qset_k(kern, i, shape_get_vertex_length(shape, i));
                kernels_qset_mu(kern, i, shape_get_vertex_orientation(shape, i));

                for (size_t j = i + 1; j < shape -> dim; j++)
                    kernels_qset_nu(kern, i, j, shape_get_vertex_angle(shape, i, j));
              }

            /* Calculate the bin-volume */
            double vol[3];
            avrFunc(specArg, vol);

            /* Insert the bin-volume */
            shape_set_volume(shape, vol[0]);
            shape_set_volume_err(shape, vol[1]);
          }

        /* Free memory */
        specArg = spec_arg_free(specArg);

      }

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     AVERAGE COVARIANCE MATRICES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Direct Function   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double avr_cov_direct(int (*avr)(cov_arg_t*, double*), cov_arg_t *covArg)
{
    /*

        Return the result of avr directly

    */

    /* Get the result */
    double result[3];
    avr(covArg, result);

    return result[0];
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   PP Covariance Matrix Average   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_cov_pp_gauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged gaussian pp covariance matrix element

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_qget_k(kern, 0);
    double mu = kernels_qget_mu(kern, 0);

    double dk1 = covArg -> dk1;
    double dmu1 = covArg -> dmu1;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovPPGaussIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovPPGaussIntgrt -> dim);

    upperBounds[0] = k + dk1 / 2.;
    upperBounds[1] = mu + dmu1 / 2.;

    lowerBounds[0] = k - dk1 / 2.;
    lowerBounds[1] = mu - dmu1 / 2.;

    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;
    avrParams -> bounds = NULL;
    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovPPGaussIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);

    /* Integrate */
    integrate(integrand_average_cov_pp_gauss, intgrt, result);

    result[0] /= pow(2. * M_PI, 3.) * shape_get_volume(covArg -> shape1);
    result[1] /= pow(2. * M_PI, 3.) * shape_get_volume(covArg -> shape1);


    /* Reset variables */
    kernels_qset_k(kern, 0, k);
    kernels_qset_k(kern, 1, k);

    kernels_qset_mu(kern, 0, mu);
    kernels_qset_mu(kern, 1, -mu);

    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


int avr_cov_pp_gauss_inf(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged gaussian pp covariance matrix element for infinitesimally small bins in k and mu

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* Get the result */
    result[0] = pow(_specPnl_(specArg, NULL), 2.);
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_cov_pp_ngauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged non-gaussian pp covariance matrix element

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 2);

    double dk1 = covArg -> dk1;
    double dk2 = covArg -> dk2;
    double dmu1 = covArg -> dmu1;
    double dmu2 = covArg -> dmu2;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovPPNGaussIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovPPNGaussIntgrt -> dim);

    upperBounds[0] = k1 + dk1 / 2.;
    upperBounds[1] = k2 + dk2 / 2.;
    upperBounds[2] = mu1 + dmu1 / 2.;
    upperBounds[3] = mu2 + dmu2 / 2.;
    upperBounds[4] = 2. * M_PI;

    lowerBounds[0] = k1 - dk1 / 2.;
    lowerBounds[1] = k2 - dk2 / 2.;
    lowerBounds[2] = mu1 - dmu1 / 2.;
    lowerBounds[3] = mu2 - dmu2 / 2.;
    lowerBounds[4] = 0.;

    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;
    avrParams -> bounds = NULL;
    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovPPNGaussIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);


    /* Integrate */
    integrate(integrand_average_cov_pp_ngauss, intgrt, result);

    result[0] /= pow(2. * M_PI, 6.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);
    result[1] /= pow(2. * M_PI, 6.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);


    /* Reset variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_k(kern, 2, k2);
    kernels_qset_k(kern, 3, k2);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_mu(kern, 2, mu2);
    kernels_qset_mu(kern, 3, -mu2);

    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


int avr_cov_pp_ngauss_inf(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged non-gaussian pp covariance matrix element for infinitesimally small bins in k and mu

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovPPNGaussInfIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovPPNGaussInfIntgrt -> dim);

    upperBounds[0] = 2. * M_PI;
    lowerBounds[0] = 0.;

    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;
    avrParams -> bounds = NULL;
    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovPPNGaussInfIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);

    /* Integrate */
    integrate(integrand_average_cov_pp_ngauss_inf, intgrt, result);

    result[0] /= 2. * M_PI;
    result[1] /= 2. * M_PI;


    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   BB Covariance Matrix Average   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_cov_bb_gauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged gaussian bb covariance matrix element

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double dk1 = covArg -> dk1;
    double dmu1 = covArg -> dmu1;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovBBGaussIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovBBGaussIntgrt -> dim);

    upperBounds[0] = k1 + dk1 / 2.;
    upperBounds[1] = k2 + dk1 / 2.;
    upperBounds[2] = mu1 + dmu1 / 2.;
    upperBounds[3] = mu2 + dmu1 / 2.;
    upperBounds[4] = 2. * M_PI;

    lowerBounds[0] = k1 - dk1 / 2.;
    lowerBounds[1] = k2 - dk1 / 2.;
    lowerBounds[2] = mu1 - dmu1 / 2.;
    lowerBounds[3] = mu2 - dmu1 / 2.;
    lowerBounds[4] = 0.;


    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;

    double bounds[4] = {k3 - dk1 / 2., k3 + dk1 / 2., mu3 - dmu1 / 2., mu3 + dmu1 / 2.};
    avrParams -> bounds = bounds;

    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovBBGaussIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);


    /* Integrate */
    integrate(integrand_average_cov_bb_gauss, intgrt, result);

    result[0] /= pow(2. * M_PI, 6.) * shape_get_volume(covArg -> shape1);
    result[1] /= pow(2. * M_PI, 6.) * shape_get_volume(covArg -> shape1);


    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


int avr_cov_bb_gauss_inf(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged gaussian bb covariance matrix element for infinitesimally small bins in k and mu

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);


    /* Results */

    /* Pnl1 */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);

    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl1 = _specPnl_(specArg, NULL);

    /* Pnl2 */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);

    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);

    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl2 = _specPnl_(specArg, NULL);

    /* Pnl3 */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);

    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);

    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl3 = _specPnl_(specArg, NULL);


    /* Result */

    result[0] = pnl1 * pnl2 * pnl3;
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_cov_bb_ngauss_tp(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged trispectrum x power spectrum part of the non-gaussian bb covariance matrix element

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1_1 = kernels_qget_k(kern, 0);
    double k1_2 = kernels_qget_k(kern, 1);
    double k1_3 = kernels_qget_k(kern, 2);

    double mu1_1 = kernels_qget_mu(kern, 0);
    double mu1_2 = kernels_qget_mu(kern, 1);
    double mu1_3 = kernels_qget_mu(kern, 2);

//    double k2_1 = kernels_qget_k(kern, 3); // ~ k1_1
    double k2_2 = kernels_qget_k(kern, 4);
    double k2_3 = kernels_qget_k(kern, 5);

//    double mu2_1 = kernels_qget_mu(kern, 3); // ~ mu1_1
    double mu2_2 = kernels_qget_mu(kern, 4);
    double mu2_3 = kernels_qget_mu(kern, 5);


    double dk1 = covArg -> dk1;
    double dk2 = covArg -> dk2;
    double dmu1 = covArg -> dmu1;
    double dmu2 = covArg -> dmu2;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovBBNGaussIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovBBNGaussIntgrt -> dim);

    upperBounds[0] = k1_1 + dk1 / 2.;
    upperBounds[1] = k1_2 + dk1 / 2.;
    upperBounds[2] = k2_2 + dk2 / 2.;
    upperBounds[3] = mu1_1 + dmu1 / 2.;
    upperBounds[4] = mu1_2 + dmu1 / 2.;
    upperBounds[5] = mu2_2 + dmu2 / 2.;
    upperBounds[6] = 2. * M_PI;
    upperBounds[7] = 2. * M_PI;

    lowerBounds[0] = k1_1 - dk1 / 2.;
    lowerBounds[1] = k1_2 - dk1 / 2.;
    lowerBounds[2] = k2_2 - dk2 / 2.;
    lowerBounds[3] = mu1_1 - dmu1 / 2.;
    lowerBounds[4] = mu1_2 - dmu1 / 2.;
    lowerBounds[5] = mu2_2 - dmu2 / 2.;
    lowerBounds[6] = 0.;
    lowerBounds[7] = 0.;


    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;

    double bounds[8] = {k1_3 - dk1 / 2., k1_3 + dk1 / 2., mu1_3 - dmu1 / 2., mu1_3 + dmu1 / 2.,
                        k2_3 - dk2 / 2., k2_3 + dk2 / 2., mu2_3 - dmu2 / 2., mu2_3 + dmu2 / 2.};
    avrParams -> bounds = bounds;

    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovBBNGaussIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);


    /* Integrate */
    integrate(integrand_average_cov_bb_ngauss_tp, intgrt, result);

    result[0] /= pow(2. * M_PI, 9.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);
    result[1] /= pow(2. * M_PI, 9.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);


    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


int avr_cov_bb_ngauss_tp_inf(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged trispectrum x power spectrum part of the non-gaussian bb covariance matrix element for infinitesimally small bins in k and mu

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1_1 = kernels_qget_k(kern, 0);
    double k1_2 = kernels_qget_k(kern, 1);
    double k1_3 = kernels_qget_k(kern, 2);

    double mu1_1 = kernels_qget_mu(kern, 0);
    double mu1_2 = kernels_qget_mu(kern, 1);
    double mu1_3 = kernels_qget_mu(kern, 2);

    double nu1_12 = kernels_qget_nu(kern, 0, 1);
    double nu1_13 = kernels_qget_nu(kern, 0, 2);
    double nu1_23 = kernels_qget_nu(kern, 1, 2);

//    double k2_1 = kernels_qget_k(kern, 3); // ~ k1_1
    double k2_2 = kernels_qget_k(kern, 4);
    double k2_3 = kernels_qget_k(kern, 5);

//    double mu2_1 = kernels_qget_mu(kern, 3); // ~ mu1_1
    double mu2_2 = kernels_qget_mu(kern, 4);
    double mu2_3 = kernels_qget_mu(kern, 5);

    double nu2_12 = kernels_qget_nu(kern, 3, 4);
    double nu2_13 = kernels_qget_nu(kern, 3, 5);
    double nu2_23 = kernels_qget_nu(kern, 4, 5);


    /* Additional variables */
    double cosPhi1_2 = (nu1_12 - mu1_1 * mu1_2) / (sqrt(1. - mu1_1*mu1_1) * sqrt(1. - mu1_2*mu1_2));
    double cosPhi1_3 = (nu1_13 - mu1_1 * mu1_3) / (sqrt(1. - mu1_1*mu1_1) * sqrt(1. - mu1_3*mu1_3));

    double sinPhi1_2 = (fabs(fabs(cosPhi1_2) - 1.) < __ABSTOL__) ? 0. : sqrt(1. - cosPhi1_2*cosPhi1_2);
    double sinPhi1_3 = (fabs(fabs(cosPhi1_3) - 1.) < __ABSTOL__) ? 0. : -sqrt(1. - cosPhi1_3*cosPhi1_3);


    double cosPhi2_2 = (nu2_12 - mu1_1 * mu2_2) / (sqrt(1. - mu1_1*mu1_1) * sqrt(1. - mu2_2*mu2_2));
    double cosPhi2_3 = (nu2_13 - mu1_1 * mu2_3) / (sqrt(1. - mu1_1*mu1_1) * sqrt(1. - mu2_3*mu2_3));

    double sinPhi2_2 = (fabs(fabs(cosPhi2_2) - 1.) < __ABSTOL__) ? 0. : sqrt(1. - cosPhi2_2*cosPhi2_2);
    double sinPhi2_3 = (fabs(fabs(cosPhi2_3) - 1.) < __ABSTOL__) ? 0. :  -sqrt(1. - cosPhi2_3*cosPhi2_3);


    /* Results */

    /* Ttr */
    kernels_qset_k(kern, 0, k1_2);
    kernels_qset_k(kern, 1, k1_3);
    kernels_qset_k(kern, 2, k2_2);
    kernels_qset_k(kern, 3, k2_3);

    kernels_qset_mu(kern, 0, mu1_2);
    kernels_qset_mu(kern, 1, mu1_3);
    kernels_qset_mu(kern, 2, -mu2_2);
    kernels_qset_mu(kern, 3, -mu2_3);

    kernels_qset_nu(kern, 0, 1, nu1_23);
    kernels_qset_nu(kern, 2, 3, nu2_23);

    double ttr = 0.;

    for (size_t i = 0; i < 4; i++)
      {
        double nu12_22 = (cosPhi1_2*cosPhi2_2 + sinPhi1_2*sinPhi2_2) * sqrt( fabs((1. - mu1_2*mu1_2) * (1. - mu2_2*mu2_2)) ) + mu1_2*mu2_2;
        double nu12_23 = (cosPhi1_2*cosPhi2_3 + sinPhi1_2*sinPhi2_3) * sqrt( fabs((1. - mu1_2*mu1_2) * (1. - mu2_3*mu2_3)) ) + mu1_2*mu2_3;
        double nu12_32 = (cosPhi1_3*cosPhi2_2 + sinPhi1_3*sinPhi2_2) * sqrt( fabs((1. - mu1_3*mu1_3) * (1. - mu2_2*mu2_2)) ) + mu1_3*mu2_2;
        double nu12_33 = (cosPhi1_3*cosPhi2_3 + sinPhi1_3*sinPhi2_3) * sqrt( fabs((1. - mu1_3*mu1_3) * (1. - mu2_3*mu2_3)) ) + mu1_3*mu2_3;


        nu12_22 = (fabs(nu12_22) - 1. < 0.) ? nu12_22
                                            : (nu12_22 > 0.) ? 1.
                                                             : -1.;
        nu12_23 = (fabs(nu12_23) - 1. < 0.) ? nu12_23
                                            : (nu12_23 > 0.) ? 1.
                                                             : -1.;
        nu12_32 = (fabs(nu12_32) - 1. < 0.) ? nu12_32
                                            : (nu12_32 > 0.) ? 1.
                                                             : -1.;
        nu12_33 = (fabs(nu12_33) - 1. < 0.) ? nu12_33
                                            : (nu12_33 > 0.) ? 1.
                                                             : -1.;

        kernels_qset_nu(kern, 0, 2, -nu12_22);
        kernels_qset_nu(kern, 0, 3, -nu12_23);
        kernels_qset_nu(kern, 1, 2, -nu12_32);
        kernels_qset_nu(kern, 1, 3, -nu12_33);

        ttr += _specTtr_(specArg, NULL) / 4.;

        sinPhi1_2 *= (i % 2 == 0) ? -1. : 1.;
        sinPhi1_3 *= (i % 2 == 0) ? -1. : 1.;
        sinPhi2_2 *= (i % 2 == 1) ? -1. : 1.;
        sinPhi2_3 *= (i % 2 == 1) ? -1. : 1.;
      }

    /* Pnl */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_1);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, -mu1_1);

    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl = _specPnl_(specArg, NULL);

    /* Line-Bin-Volume */
    specArg -> dk = covArg -> dk2;
    specArg -> dmu = covArg -> dmu2;

    double lineVol = avr_shape_vol_direct(avr_shape_line_vol_ana, specArg);


    /* Get the result */
    result[0] = pnl * ttr / lineVol;
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int avr_cov_bb_ngauss_bb(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged bispectrum x bispectrum part of the non-gaussian bb covariance matrix element

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1_1 = kernels_qget_k(kern, 0);
    double k1_2 = kernels_qget_k(kern, 1);
    double k1_3 = kernels_qget_k(kern, 2);

    double mu1_1 = kernels_qget_mu(kern, 0);
    double mu1_2 = kernels_qget_mu(kern, 1);
    double mu1_3 = kernels_qget_mu(kern, 2);


//    double k2_1 = kernels_qget_k(kern, 3); // ~ k1_1
    double k2_2 = kernels_qget_k(kern, 4);
    double k2_3 = kernels_qget_k(kern, 5);

//    double mu2_1 = kernels_qget_mu(kern, 3); // ~ -mu1_1
    double mu2_2 = kernels_qget_mu(kern, 4);
    double mu2_3 = kernels_qget_mu(kern, 5);

    double dk1 = covArg -> dk1;
    double dk2 = covArg -> dk2;
    double dmu1 = covArg -> dmu1;
    double dmu2 = covArg -> dmu2;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovBBNGaussIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovBBNGaussIntgrt -> dim);

    upperBounds[0] = k1_1 + dk1 / 2.;
    upperBounds[1] = k1_2 + dk1 / 2.;
    upperBounds[2] = k2_2 + dk2 / 2.;
    upperBounds[3] = mu1_1 + dmu1 / 2.;
    upperBounds[4] = mu1_2 + dmu1 / 2.;
    upperBounds[5] = mu2_2 + dmu2 / 2.;
    upperBounds[6] = 2. * M_PI;
    upperBounds[7] = 2. * M_PI;

    lowerBounds[0] = k1_1 - dk1 / 2.;
    lowerBounds[1] = k1_2 - dk1 / 2.;
    lowerBounds[2] = k2_2 - dk2 / 2.;
    lowerBounds[3] = mu1_1 - dmu1 / 2.;
    lowerBounds[4] = mu1_2 - dmu1 / 2.;
    lowerBounds[5] = mu2_2 - dmu2 / 2.;
    lowerBounds[6] = 0.;
    lowerBounds[7] = 0.;


    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;

    double bounds[8] = {k1_3 - dk1 / 2., k1_3 + dk1 / 2., mu1_3 - dmu1 / 2., mu1_3 + dmu1 / 2.,
                        k2_3 - dk2 / 2., k2_3 + dk2 / 2., mu2_3 - dmu2 / 2., mu2_3 + dmu2 / 2.};
    avrParams -> bounds = bounds;

    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovBBNGaussIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);


    /* Integrate */
    integrate(integrand_average_cov_bb_ngauss_bb, intgrt, result);

    result[0] /= pow(2. * M_PI, 9.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);
    result[1] /= pow(2. * M_PI, 9.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);


    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


int avr_cov_bb_ngauss_bb_inf(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged bispectrum x bispectrum part of the non-gaussian bb covariance matrix element for infinitesimally small bins in k and mu

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1_1 = kernels_qget_k(kern, 0);
    double k1_2 = kernels_qget_k(kern, 1);
    double k1_3 = kernels_qget_k(kern, 2);

    double mu1_1 = kernels_qget_mu(kern, 0);
    double mu1_2 = kernels_qget_mu(kern, 1);
    double mu1_3 = kernels_qget_mu(kern, 2);

    double nu1_12 = kernels_qget_nu(kern, 0, 1);
    double nu1_13 = kernels_qget_nu(kern, 0, 2);
    double nu1_23 = kernels_qget_nu(kern, 1, 2);


    double k2_1 = kernels_qget_k(kern, 3);
    double k2_2 = kernels_qget_k(kern, 4);
    double k2_3 = kernels_qget_k(kern, 5);

    double mu2_1 = kernels_qget_mu(kern, 3);
    double mu2_2 = kernels_qget_mu(kern, 4);
    double mu2_3 = kernels_qget_mu(kern, 5);

    double nu2_12 = kernels_qget_nu(kern, 3, 4);
    double nu2_13 = kernels_qget_nu(kern, 3, 5);
    double nu2_23 = kernels_qget_nu(kern, 4, 5);


    /* Results */

    /* Btr1 */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_2);
    kernels_qset_k(kern, 2, k1_3);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, mu1_2);
    kernels_qset_mu(kern, 2, mu1_3);

    kernels_qset_nu(kern, 0, 1, nu1_12);
    kernels_qset_nu(kern, 0, 2, nu1_13);
    kernels_qset_nu(kern, 1, 2, nu1_23);

    double btr1 = _specBtr_(specArg, NULL);

    /* Btr1 */
    kernels_qset_k(kern, 0, k2_1);
    kernels_qset_k(kern, 1, k2_2);
    kernels_qset_k(kern, 2, k2_3);

    kernels_qset_mu(kern, 0, mu2_1);
    kernels_qset_mu(kern, 1, mu2_2);
    kernels_qset_mu(kern, 2, mu2_3);

    kernels_qset_nu(kern, 0, 1, nu2_12);
    kernels_qset_nu(kern, 0, 2, nu2_13);
    kernels_qset_nu(kern, 1, 2, nu2_23);

    double btr2 = _specBtr_(specArg, NULL);

    /* Line-Bin-Volume */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_1);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, -mu1_1);

    kernels_qset_nu(kern, 0, 1, -1.);

    specArg -> dk = covArg -> dk2;
    specArg -> dmu = covArg -> dmu2;

    double lineVol = avr_shape_vol_direct(avr_shape_line_vol_ana, specArg);


    /* Get the result */
    result[0] = btr1 * btr2 / lineVol;
    result[1] = 0.;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   BB Covariance Matrix Average   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int avr_cov_pb_ngauss_bp(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged bispectrum x power specturm part of the non-gaussian pb covariance matrix element

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1_1 = kernels_qget_k(kern, 0);
//    double k1_2 = kernels_qget_k(kern, 1);

    double mu1_1 = kernels_qget_mu(kern, 0);
//    double mu1_2 = kernels_qget_mu(kern, 1);


//    double k2_1 = kernels_qget_k(kern, 2); // ~ k1_1
    double k2_2 = kernels_qget_k(kern, 3);
    double k2_3 = kernels_qget_k(kern, 4);

//    double mu2_1 = kernels_qget_mu(kern, 2); // ~ -mu1_1
    double mu2_2 = kernels_qget_mu(kern, 3);
    double mu2_3 = kernels_qget_mu(kern, 4);

    double dk1 = covArg -> dk1;
    double dk2 = covArg -> dk2;
    double dmu1 = covArg -> dmu1;
    double dmu2 = covArg -> dmu2;

    /* Integration bounds */
    double *upperBounds = malloc(sizeof(double) * _avrCovPBNGaussIntgrt -> dim);
    double *lowerBounds = malloc(sizeof(double) * _avrCovPBNGaussIntgrt -> dim);

    upperBounds[0] = k1_1 + dk1 / 2.;
    upperBounds[1] = k2_2 + dk2 / 2.;
    upperBounds[2] = mu1_1 + dmu1 / 2.;
    upperBounds[3] = mu2_2 + dmu2 / 2.;
    upperBounds[4] = 2. * M_PI;

    lowerBounds[0] = k1_1 - dk1 / 2.;
    lowerBounds[1] = k2_2 - dk2 / 2.;
    lowerBounds[2] = mu1_1 - dmu1 / 2.;
    lowerBounds[3] = mu2_2 - dmu2 / 2.;
    lowerBounds[4] = 0.;


    /* Integrand parameters */
    cintgrnd_avr_t *avrParams = malloc(sizeof(cintgrnd_avr_t));

    avrParams -> var = specArg;
    avrParams -> params = NULL;

    double bounds[4] = {k2_3 - dk2 / 2., k2_3 + dk2 / 2., mu2_3 - dmu2 / 2., mu2_3 + dmu2 / 2.};
    avrParams -> bounds = bounds;

    avrParams -> specFunc = NULL;

    /* Integrate struct */
    intgrt_t *intgrt = integrate_cp(_avrCovPBNGaussIntgrt);

    integrate_set_bounds_upper(intgrt, upperBounds);
    integrate_set_bounds_lower(intgrt, lowerBounds);
    integrate_set_params(intgrt, avrParams);


    /* Integrate */
    integrate(integrand_average_cov_pb_ngauss_bp, intgrt, result);

    result[0] /= pow(2. * M_PI, 6.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);
    result[1] /= pow(2. * M_PI, 6.) * shape_get_volume(covArg -> shape1) * shape_get_volume(covArg -> shape2);


    /* Free memory */
    free(upperBounds);
    free(lowerBounds);

    free(avrParams);

    intgrt = integrate_free(intgrt);

    return 0;
}


int avr_cov_pb_ngauss_bp_inf(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the bin-averaged bispectrum x power spectrum part of the non-gaussian pb covariance matrix element for infinitesimally small bins in k and mu

    */

    /* specArg */
    spec_arg_t *specArg = covArg -> specArg;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1_1 = kernels_qget_k(kern, 0);
    double k1_2 = kernels_qget_k(kern, 1);

    double mu1_1 = kernels_qget_mu(kern, 0);
    double mu1_2 = kernels_qget_mu(kern, 1);

    double nu1_12 = kernels_qget_nu(kern, 0, 1);


    double k2_1 = kernels_qget_k(kern, 2);
    double k2_2 = kernels_qget_k(kern, 3);
    double k2_3 = kernels_qget_k(kern, 4);

    double mu2_1 = kernels_qget_mu(kern, 2);
    double mu2_2 = kernels_qget_mu(kern, 3);
    double mu2_3 = kernels_qget_mu(kern, 4);

    double nu2_12 = kernels_qget_nu(kern, 2, 3);
    double nu2_13 = kernels_qget_nu(kern, 2, 4);
    double nu2_23 = kernels_qget_nu(kern, 3, 4);


    /* Results */

    /* Pnl */
    kernels_qset_k(kern, 0, k1_1);
    kernels_qset_k(kern, 1, k1_2);

    kernels_qset_mu(kern, 0, mu1_1);
    kernels_qset_mu(kern, 1, mu1_2);

    kernels_qset_nu(kern, 0, 1, nu1_12);

    double pnl = _specPnl_(specArg, NULL);

    /* Btr */
    kernels_qset_k(kern, 0, k2_1);
    kernels_qset_k(kern, 1, k2_2);
    kernels_qset_k(kern, 2, k2_3);

    kernels_qset_mu(kern, 0, mu2_1);
    kernels_qset_mu(kern, 1, mu2_2);
    kernels_qset_mu(kern, 2, mu2_3);

    kernels_qset_nu(kern, 0, 1, nu2_12);
    kernels_qset_nu(kern, 0, 2, nu2_13);
    kernels_qset_nu(kern, 1, 2, nu2_23);

    double btr = _specBtr_(specArg, NULL);

    /* Line-Bin-Volume */
    kernels_qset_k(kern, 0, k2_1);
    kernels_qset_k(kern, 1, k2_1);

    kernels_qset_mu(kern, 0, mu2_1);
    kernels_qset_mu(kern, 1, -mu2_1);

    kernels_qset_nu(kern, 0, 1, -1.);

    specArg -> dk = covArg -> dk2; // dk comes from the triangle bins
    specArg -> dmu = covArg -> dmu2; // dmu comes from the triangle bins

    double lineVol = avr_shape_vol_direct(avr_shape_line_vol_ana, specArg);


    /* Get the result */
    result[0] = pnl * btr / lineVol;
    result[1] = 0.;

    return 0;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

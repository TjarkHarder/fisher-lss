/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     FLSS.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "flss.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PRE-INITIALISE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------    Pre-Initialise Functions    ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double _zeroFunc_()
{
    /*

        Returns (double) zero, regardless of input

    */

    return 0.;
}


double _oneFunc_()
{
    /*

        Returns (double) one, regardless of input

    */

    return 1.;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------    Pre-Initialise Variables    ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  True + False  **/

const bool _true_ = true;
const bool _false_ = false;


/**  Random Number Generator  **/

gsl_rng *_randomNumberGen_ = NULL;


/**  Miscellaneous Identifiers  **/

const char *_idIntError_ = "[int. error]";

const char *_idTree_ = "tree";
const char *_idSn_ = "sn";
const char *_idCtr_ = "ctr";

const char *_idGauss_ = "gaussian";
const char *_idNGauss_ = "non-gaussian";

const char *_idFull_ = "full";

const char *_idInf_ = "inf";


/**  Integration Routines  **/

const char *_idIntgrtCQUAD_ = "cquad";
const char *_idIntgrtVegas_ = "vegas";
const char *_idIntgrtDivonne_ = "divonne";


/**  Shape Identifiers  **/

const char *_idShapeLine_ = "line";
const char *_idShapeTri_ = "tri";


/**  Variable Identifiers  **/

const char *_idVarZ_ = "z";
const char *_idVarK_ = "k";
const char *_idVarMu_ = "mu";


/**  Fisher Parameter Identifiers  **/

const char *_idParamsA2Ga_ = "a2ga";
const char *_idParamsD2Ga_ = "d2ga";
const char *_idParamsA3GaA_ = "a3gaa";
const char *_idParamsA3GaB_ = "a3gab";
const char *_idParamsD3GaA_ = "d3gaa";
const char *_idParamsD3GaB_ = "d3gab";

const char *_idParamsB1_ = "b1";
const char *_idParamsB2_ = "b2";
const char *_idParamsC2Ga_ = "c2ga";
const char *_idParamsBGam3_ = "bgam3";

const char *_idParamsF_ = "f";
const char *_idParamsSigv_ = "sigv";

const char *_idParamsD_ = "d";
const char *_idParamsH_ = "h";

const char *_idParamsC0_ = "c0";
const char *_idParamsC2_ = "c2";
const char *_idParamsC4_ = "c4";

const char *_idParamsPsn_ = "psn";
const char *_idParamsBsn1_ = "bsn1";
const char *_idParamsBsn2_ = "bsn2";

const char *_idParamsPk_ = "pk";


/**  Spectra Identifiers  **/

const char *_idSpecPnl_ = "pnl";
const char *_idSpecDPnl_ = "dpnl";

const char *_idSpecBtr_ = "btr";
const char *_idSpecDBtr_ = "dbtr";

const char *_idSpecTtr_ = "ttr";


/**  Covariance Matrix Identifiers  **/

const char *_idCovPP_ = "covpp";
const char *_idCovPB_ = "covpb";
const char *_idCovBB_ = "covbb";


/**  Fisher Matrix Identifiers  **/

const char *_idFishPP_ = "fishpp";
const char *_idFishPB_ = "fishpb";
const char *_idFishBB_ = "fishbb";


/**  Average Identifiers  **/

const char *_idAvrLine_ = "avr-line";
const char *_idAvrTri_ = "avr-tri";

const char *_idAvrLineVol_ = "avr-line-vol";
const char *_idAvrTriVol_ = "avr-tri-vol";

const char *_idAvrLineVolFuncNum_ = "avr-line-vol-num";
const char *_idAvrLineVolFuncAna_ = "avr-line-vol-ana";

const char *_idAvrTriVolFuncNum_ = "avr-tri-vol-num";
const char *_idAvrTriVolFuncAna_ = "avr-tri-vol-ana";

const char *_idAvrCovPPGauss_ = "avr-covpp-gaussian";
const char *_idAvrCovPPNGaussInf_ = "avr-covpp-non-gaussian-inf";
const char *_idAvrCovPPNGauss_ = "avr-covpp-non-gaussian";

const char *_idAvrCovBBGauss_ = "avr-covbb-gaussian";
const char *_idAvrCovBBNGauss_ = "avr-covbb-non-gaussian";

const char *_idAvrCovPBNGauss_ = "avr-covpb-non-gaussian";


/**  Fisher Parameter Multiplicities  **/

/* mult = 0 : no multiplicity, ie parameter dθ(x)/dx = 0 is a constant wrt the variable x          */
/* mult != 0 : multiplicity, ie parameter dθ(x)/dx != 0 is not a constant wrt the variable x       */
/*     - mult = 1 : one-loop power spectrum depends on integral of parameter θ over the variable x */
/*     - mult = -1 : one-loop power spectrum only depends on parameter θ at a single point x       */

const int _multParamsA2Ga_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsD2Ga_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsA3GaA_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsA3GaB_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsD3GaA_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsD3GaB_[__ID_VAR_SIZE__] = {-1, 0, 0};

const int _multParamsB1_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsB2_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsC2Ga_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsBGam3_[__ID_VAR_SIZE__] = {-1, 0, 0};

const int _multParamsF_[__ID_VAR_SIZE__] = {1, 0, 0};
const int _multParamsSigv_[__ID_VAR_SIZE__] = {-1, 0, 0};

const int _multParamsD_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsH_[__ID_VAR_SIZE__] = {-1, 0, 0};

const int _multParamsC0_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsC2_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsC4_[__ID_VAR_SIZE__] = {-1, 0, 0};

const int _multParamsPsn_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsBsn1_[__ID_VAR_SIZE__] = {-1, 0, 0};
const int _multParamsBsn2_[__ID_VAR_SIZE__] = {-1, 0, 0};

const int _multParamsPk_[__ID_VAR_SIZE__] = {0, 1, 0};


/**  Flags  **/

bool _fidInclBias_ = true;
bool _fidInclRSD_ = true;
bool _fidInclCtr_ = true;
bool _fidInclSn_ = true;

bool _avrShapeLine_ = false;
bool _avrShapeTri_ = false;


/**  Sample  **/

/* Redshift (Raw) Sample */
sample_raw_t *_sampleRedshift_ = NULL;

/* Line (Shape) Sample */
sample_shape_t *_sampleShapeLine_ = NULL;

/* Triangle (Shape) Sample */
sample_shape_t *_sampleShapeTri_ = NULL;


/**  Fiducial Functions  **/

/* LCDM */

fid_lcdm_t *(*_fidLCDM_)(void*, void*) = NULL;
void *_fidParamsLCDM_ = NULL;

/*  ----------------------------------------------------  */

double (*_fidLCDMxOmegaM0_)(void*, void*) = _zeroFunc_;
void *_fidParamsLCDMxOmegaM0_ = NULL;

double (*_fidLCDMxNs_)(void*, void*) = _zeroFunc_;
void *_fidParamsLCDMxNs_ = NULL;

double (*_fidLCDMxGrowthIndex_)(void*, void*) = _zeroFunc_;
void *_fidParamsLCDMxGrowthIndex_ = NULL;


/* Bootstrap */

fid_btst_t *(*_fidBTST_)(void*, void*) = NULL;
void *_fidParamsBTST_ = NULL;

/*  ----------------------------------------------------  */

double (*_fidBTSTxA2Ga_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxA2Ga_ = NULL;

double (*_fidBTSTxD2Ga_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxD2Ga_ = NULL;


double (*_fidBTSTxA3GaA_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxA3GaA_ = NULL;

double (*_fidBTSTxA3GaB_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxA3GaB_ = NULL;

double (*_fidBTSTxD3GaA_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxD3GaA_ = NULL;

double (*_fidBTSTxD3GaB_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxD3GaB_ = NULL;


double (*_fidBTSTxH_)(void*, void*) = _zeroFunc_;
void *_fidParamsBTSTxH_ = NULL;


/* Bias */

fid_bias_t *(*_fidBias_)(void*, void*) = NULL;
void *_fidParamsBias_ = NULL;

/*  ----------------------------------------------------  */

double (*_fidBiasxB1_)(void*, void*) = _zeroFunc_;
void *_fidParamsBiasxB1_ = NULL;

double (*_fidBiasxB2_)(void*, void*) = _zeroFunc_;
void *_fidParamsBiasxB2_ = NULL;


double (*_fidBiasxBG2_)(void*, void*) = _zeroFunc_;
void *_fidParamsBiasxBG2_ = NULL;

double (*_fidBiasxC2Ga_)(void*, void*) = _zeroFunc_;
void *_fidParamsBiasxC2Ga_ = NULL;


double (*_fidBiasxBGam3_)(void*, void*) = _zeroFunc_;
void *_fidParamsBiasxBGam3_ = NULL;


/* RSD */

fid_rsd_t *(*_fidRSD_)(void*, void*) = NULL;
void *_fidParamsRSD_ = NULL;

/*  ----------------------------------------------------  */

double (*_fidRSDxF_)(void*, void*) = _zeroFunc_;
void *_fidParamsRSDxF_ = NULL;


double (*_fidRSDxSigv_)(void*, void*) = _zeroFunc_;
void *_fidParamsRSDxSigv_ = NULL;

double (*_fidRSDxSigs_)(void*, void*) = _zeroFunc_;
void *_fidParamsRSDxSigs_ = NULL;


/* Ctr */

fid_ctr_t *(*_fidCtr_)(void*, void*) = NULL;
void *_fidParamsCtr_ = _zeroFunc_;

/*  ----------------------------------------------------  */

double (*_fidCtrxC0_)(void*, void*) = _zeroFunc_;
void *_fidParamsCtrxC0_ = NULL;

double (*_fidCtrxC2_)(void*, void*) = _zeroFunc_;
void *_fidParamsCtrxC2_ = NULL;

double (*_fidCtrxC4_)(void*, void*) = _zeroFunc_;
void *_fidParamsCtrxC4_ = NULL;


/* Survey */

fid_surv_t *(*_fidSurv_)(void*, void*) = NULL;
void *_fidParamsSurv_ = NULL;

/*  ----------------------------------------------------  */

double (*_fidSurvxV_)(void*, void*) = _zeroFunc_;
void *_fidParamsSurvxV_ = NULL;

double (*_fidSurvxN_)(void*, void*) = _zeroFunc_;
void *_fidParamsSurvxN_ = NULL;

double (*_fidSurvxSn_)(void*, void*) = _zeroFunc_;
void *_fidParamsSurvxSn_ = NULL;


/* Growth */

double (*_fidGrowth_)(void*, void*) = _zeroFunc_;
void *_fidParamsGrowth_ = NULL;


/* Linear Power Spectrum */

double (*_fidPk_)(void*, void*) = _zeroFunc_;
double (*_fidExtrapPk_)(void*, void*, void*, size_t*, size_t) = _zeroFunc_;
void *_fidParamsPk_ = NULL;

double (*_fidDPk_)(void*, void*) = _zeroFunc_;
double (*_fidExtrapDPk_)(void*, void*, void*, size_t*, size_t) = _zeroFunc_;
void *_fidParamsDPk_ = NULL;


/**  Spectra Functions  **/

/* Power Spectrum */

double (*_specPtr_)(void*, void*) = _zeroFunc_;
double (*_specPnl_)(void*, void*) = _zeroFunc_;


double (*_specDPnlA2Ga_)(void*, void*) = _zeroFunc_;
double (*_specDPnlD2Ga_)(void*, void*) = _zeroFunc_;
double (*_specDPnlA3GaA_)(void*, void*) = _zeroFunc_;
double (*_specDPnlA3GaB_)(void*, void*) = _zeroFunc_;
double (*_specDPnlD3GaA_)(void*, void*) = _zeroFunc_;
double (*_specDPnlD3GaB_)(void*, void*) = _zeroFunc_;

double (*_specDPnlB1_)(void*, void*) = _zeroFunc_;
double (*_specDPnlB2_)(void*, void*) = _zeroFunc_;
double (*_specDPnlC2Ga_)(void*, void*) = _zeroFunc_;
double (*_specDPnlBGam3_)(void*, void*) = _zeroFunc_;

double (*_specDPnlF_)(void*, void*) = _zeroFunc_;
double (*_specDPnlSigv_)(void*, void*) = _zeroFunc_;

double (*_specDPnlD_)(void*, void*) = _zeroFunc_;
double (*_specDPnlH_)(void*, void*) = _zeroFunc_;

double (*_specDPnlC0_)(void*, void*) = _zeroFunc_;
double (*_specDPnlC2_)(void*, void*) = _zeroFunc_;
double (*_specDPnlC4_)(void*, void*) = _zeroFunc_;

double (*_specDPnlPsn_)(void*, void*) = _zeroFunc_;
double (*_specDPnlBsn1_)(void*, void*) = _zeroFunc_;
double (*_specDPnlBsn2_)(void*, void*) = _zeroFunc_;

double (*_specDPnlPk_)(void*, void*) = _zeroFunc_;

double (*(*_specDPnl_)(const char*))(void*, void*) = NULL;


/* Bispectrum */

double (*_specBtr_)(void*, void*) = _zeroFunc_;


double (*_specDBtrA2Ga_)(void*, void*) = _zeroFunc_;
double (*_specDBtrD2Ga_)(void*, void*) = _zeroFunc_;
double (*_specDBtrA3GaA_)(void*, void*) = _zeroFunc_;
double (*_specDBtrA3GaB_)(void*, void*) = _zeroFunc_;
double (*_specDBtrD3GaA_)(void*, void*) = _zeroFunc_;
double (*_specDBtrD3GaB_)(void*, void*) = _zeroFunc_;

double (*_specDBtrB1_)(void*, void*) = _zeroFunc_;
double (*_specDBtrB2_)(void*, void*) = _zeroFunc_;
double (*_specDBtrC2Ga_)(void*, void*) = _zeroFunc_;
double (*_specDBtrBGam3_)(void*, void*) = _zeroFunc_;

double (*_specDBtrF_)(void*, void*) = _zeroFunc_;
double (*_specDBtrSigv_)(void*, void*) = _zeroFunc_;

double (*_specDBtrD_)(void*, void*) = _zeroFunc_;
double (*_specDBtrH_)(void*, void*) = _zeroFunc_;

double (*_specDBtrC0_)(void*, void*) = _zeroFunc_;
double (*_specDBtrC2_)(void*, void*) = _zeroFunc_;
double (*_specDBtrC4_)(void*, void*) = _zeroFunc_;

double (*_specDBtrPsn_)(void*, void*) = _zeroFunc_;
double (*_specDBtrBsn1_)(void*, void*) = _zeroFunc_;
double (*_specDBtrBsn2_)(void*, void*) = _zeroFunc_;

double (*_specDBtrPk_)(void*, void*) = _zeroFunc_;

double (*(*_specDBtr_)(const char*))(void*, void*) = NULL;



/* Trispectrum */

double (*_specTtr_)(void*, void*) = _zeroFunc_;


/* Poly Spectrum */

double (*(*_specPoly_)(const char*, const char*))(void*, void*) = NULL;


/**  Covariance Matrix Functions  **/

/* PP Covariance Matrix */

mat_t *(*_covPPMat_)(void*, void*) = NULL;

double (*_covPP_)(void*, void*) = _zeroFunc_;


/* BB Covariance Matrix */

mat_t *(*_covBBMat_)(void*, void*) = NULL;

double (*_covBB_)(void*, void*) = _zeroFunc_;


/* PB Covariance Matrix */

mat_t *(*_covPBMat_)(void*, void*) = NULL;

double (*_covPB_)(void*, void*) = _zeroFunc_;


/* Poly Spectrum Covariance Matrix */

mat_t *(*(*_covPolyMat_)(const char*))(void*, void*) = NULL;

double (*(*_covPoly_)(const char*, const char*))(void*, void*) = NULL;






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LOCAL PARAMETERS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/* Local parameters */

static unsigned int _threadsMaxNum;





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE MODULES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_threads(void);
static int _ini_random(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int flss_ini(void)
{
    /*

        Initialise all modules.

    */


    /* Threads */
    _ini_threads();

    /* Random Numbers */
    _ini_random();


    /* Additional ini functions */

    spec_ini();
    cov_ini();
    fish_ini();

    fiducials_ini();

    avr_ini();

    prnt_ini();


    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------    Ini Functions    -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_threads(void)
{
    /*

        Threads

    */

    _threadsMaxNum = (unsigned int) omp_get_max_threads(); /* Store the maximum number of threads */

    /* Oopen MP */
    omp_set_num_threads(1);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_random(void)
{
    /*

        Random number generator

    */

    _randomNumberGen_ = gsl_rng_alloc(gsl_rng_taus);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SETTERS AND GETTERS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------    Getters    -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *flss_get_dir(void)
{
    /*

        Get the directory of fisher-lss

    */

    return __FISHERLSSDIR__;
}


bool flss_get_avr_shape_flag(const char *id)
{
    /*

        Get the average flag for a given id

    */

    /* Line */
    if (!strcmp(id, _idAvrLine_) || !strcmp(id, _idShapeLine_) || !strcmp(id, _idSpecPnl_) || !strcmp(id, _idSpecDPnl_) || !strcmp(id, _idCovPP_) || !strcmp(id, _idFishPP_))
      {
        return _avrShapeLine_;
      }

    /* Triangle */
    if (!strcmp(id, _idAvrTri_) || !strcmp(id, _idShapeTri_) || !strcmp(id, _idSpecBtr_) || !strcmp(id, _idSpecDBtr_) || !strcmp(id, _idCovBB_) || !strcmp(id, _idFishBB_))
      {
        return _avrShapeTri_;
      }

    /* Wrong id */
    printf("Cannot get the average flag for a shape of unknown id '%s'.\n", id);
    exit(1);

    return false;
}


sample_shape_t *flss_get_sample_shape(const char *id)
{
    /*

        Get the sample for the shape that enters the spectrum of given id

    */

    /* Power Spectrum */
    if (!strcmp(id, _idShapeLine_) || !strcmp(id, _idSpecPnl_) || !strcmp(id, _idSpecDPnl_) || !strcmp(id, _idCovPP_) || !strcmp(id, _idFishPP_))
        return _sampleShapeLine_;

    /* Bispectrum */
    if (!strcmp(id, _idShapeTri_) || !strcmp(id, _idSpecBtr_) || !strcmp(id, _idSpecDBtr_) || !strcmp(id, _idCovBB_) || !strcmp(id, _idFishBB_))
        return _sampleShapeTri_;

    /* Wrong id */
    printf("Cannot get the sample of the shape related to an unknown spectrum's id %s.\n", id);
    exit(1);

    return NULL;
}


sample_raw_t *flss_get_sample_redshift(void)
{
    /*

        Get the sample for the redshift

    */

    return _sampleRedshift_;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------    Setters    -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int flss_set_threads(unsigned int numberOfThreads)
{
    /*

        Set the number of threads

    */

    if (numberOfThreads > _threadsMaxNum)
        numberOfThreads = _threadsMaxNum;

    omp_set_num_threads((int) numberOfThreads);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int flss_set_avr_shape_flag(const char *id, bool avr)
{
    /*

        Set the average flag for a given id

    */

    /* Line */
    if (!strcmp(id, _idAvrLine_) || !strcmp(id, _idShapeLine_) || !strcmp(id, _idSpecPnl_) || !strcmp(id, _idSpecDPnl_) || !strcmp(id, _idCovPP_) || !strcmp(id, _idFishPP_))
      {
        _avrShapeLine_ = avr;

        return 0;
      }

    /* Triangle */
    if (!strcmp(id, _idAvrTri_) || !strcmp(id, _idShapeTri_) || !strcmp(id, _idSpecBtr_) || !strcmp(id, _idSpecDBtr_) || !strcmp(id, _idCovBB_) || !strcmp(id, _idFishBB_))
      {
        _avrShapeTri_ = avr;

        return 0;
      }

    /* Wrong id */
    printf("Cannot set the average flag for a shape of unknown id '%s'.\n", id);
    exit(1);

    return 1;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int flss_set_sample_shape(const char *id, sample_shape_t *sampleShape)
{
    /*

        Set the sample for the shape related to the spectrum with given id.

    */

    /* Line */
    if (!strcmp(id, _idShapeLine_) || !strcmp(id, _idSpecPnl_) || !strcmp(id, _idSpecDPnl_) || !strcmp(id, _idCovPP_) || !strcmp(id, _idFishPP_))
      {
        _sampleShapeLine_ = sample_shape_free(_sampleShapeLine_);
        _sampleShapeLine_ = sample_shape_cp(sampleShape);

        return 0;
      }

    /* Triangle */
    if (!strcmp(id, _idShapeTri_) || !strcmp(id, _idSpecBtr_) || !strcmp(id, _idSpecDBtr_) || !strcmp(id, _idCovBB_) || !strcmp(id, _idFishBB_))
      {
        _sampleShapeTri_ = sample_shape_free(_sampleShapeTri_);
        _sampleShapeTri_ = sample_shape_cp(sampleShape);

        return 0;
      }

    /* Wrong id */
    printf("Cannot set the sample for the shape related to an unknown spectrum's id %s.\n", id);
    exit(1);

    return 1;
}


int flss_set_sample_redshift(sample_raw_t *sampleRedshift)
{
    /*

        Set the sample for the redshift

    */

    _sampleRedshift_ = sample_raw_free(_sampleRedshift_);
    _sampleRedshift_ = sample_raw_cp(sampleRedshift);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int flss_set_random_seed(size_t seed)
{
    /*

        Set the seed of the random number generator

    */

    gsl_rng_set(_randomNumberGen_, seed);

    return 0;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     FREE MODULES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_random(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int flss_free(void)
{
    /*

        Free variables from all modules.

    */

    /* Random Numbers */
    _free_random();

    /* Additional free functions */
    spec_free();
    cov_free();
    fish_free();

    fiducials_free();

    avr_free();

    prnt_free();

    /* Free global variables */
    _sampleRedshift_ = sample_raw_free(_sampleRedshift_);
    _sampleShapeLine_ = sample_shape_free(_sampleShapeLine_);
    _sampleShapeTri_ = sample_shape_free(_sampleShapeTri_);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------    Free Functions    ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_random(void)
{
    /*

        Random number generator

    */

    gsl_rng_free(_randomNumberGen_);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

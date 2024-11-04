/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SPECTRA.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "spectra.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE / SETUP / FREE VARIABLES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Local Variables   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Zero function for spectra  **/

static int _spec_zero(spec_arg_t *specArg, double *result)
{
    /*

        Fills result with zeros (always assuming that result has size 2 !!!)

    */

    (void) specArg;

    result[0] = 0.;
    result[1] = 0.;

    return 0;
}


/**  Output Directory and File Extension  **/

static const char *_outDir = "/output/spec/";
static const char *_outExt = ".dat";


/**  Information about Spectra  **/

typedef struct
{
    /*

        Information about the spectra

    */

    const char *id;

    size_t specOrder;

    size_t loopOrderMax;
    size_t loopOrder; // Can change this
    intgrt_t **loopIntgrt; // Can change this

    size_t partsSize;

    char **partsLabels;
    bool **partsExist;
    bool **partsHasErr;

    dat_t **partsDerivVals;
    int **partsDerivMult;
    bool *partsDerivLog;

} _spec_info_t;


/*  ------------------------------------------------------------------------------------------------------  */


static _spec_info_t *_pnlInfo = NULL;
static _spec_info_t *_dpnlInfo = NULL;

static _spec_info_t *_btrInfo = NULL;
static _spec_info_t *_dbtrInfo = NULL;

static _spec_info_t *_ttrInfo = NULL;


static _spec_info_t *_spec_info_get_struct(const char *id)
{
    /*

        Get the info struct with given id

    */

    /* Non-Linear Power Spectrum */
    if (!strcmp(id, _idSpecPnl_))
        return _pnlInfo;

    /* Non-Linear Power Spectrum Derivatives */
    if (!strcmp(id, _idSpecDPnl_))
        return _dpnlInfo;

    /* Tree-Level Bispectrum */
    if (!strcmp(id, _idSpecBtr_))
        return _btrInfo;

    /* Tree-Level Bispectrum Derivatives */
    if (!strcmp(id, _idSpecDBtr_))
        return _dbtrInfo;

    /* Tree-Level Trispectrum */
    if (!strcmp(id, _idSpecTtr_))
        return _ttrInfo;

    printf("Spectra with ID ’%s’ does not exist.\n", id);
    exit(1);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


static _spec_info_t *_spec_info_free(_spec_info_t *info)
{
    /*

        Free a _spec_info_t struct

    */

    if (info == NULL)
        return 0;

    /* Loop integration variables */
    for (size_t i = 0; i < info -> loopOrderMax; i++)
      {
        info -> loopIntgrt[i] = integrate_free(info -> loopIntgrt[i]);
      }

    free(info -> loopIntgrt);

    /* Parts variables */
    free(info -> partsLabels);
    free(info -> partsExist);
    free(info -> partsHasErr);

    /* Additional derivative variables */
    if (info -> partsDerivVals != NULL)
      {
        for (size_t i = 0; i < info -> partsSize; i++)
            info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);
      }

    free(info -> partsDerivMult);
    free(info -> partsDerivVals);
    free(info -> partsDerivLog);

    free(info);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool spec_info_isderiv(const char *id)
{
    /*

        Check if a spectrum with given ID is a "derivative"

    */

    _spec_info_t *info = _spec_info_get_struct(id);

    return info -> partsDerivLog != NULL; // Only need to check if one partsDeriv(...) pointer is NULL or not
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_info_set_loop_order(const char *id, size_t loopOrder)
{
    /*

        Set the loop order of a spectrum

    */

    _spec_info_t *info = _spec_info_get_struct(id);

    if (info -> loopOrderMax < loopOrder)
      {
        printf("Cannot set the loop order of a spectrum to %ld if only %ld loops are expected.\n", loopOrder, info -> loopOrderMax);
        exit(1);

        return 1;
      }

    info -> loopOrder = loopOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


size_t spec_info_get_order(const char *id)
{
    /*

        Get a spectrum's order

    */

    _spec_info_t *info = _spec_info_get_struct(id);

    return info -> specOrder;
}


size_t spec_info_get_loop_order(const char *id)
{
    /*

        Get the loop order of a spectrum

    */

    _spec_info_t *info = _spec_info_get_struct(id);

    return info -> loopOrder;
}


intgrt_t *spec_info_get_loop_integrate(const char *id, size_t loopOrder)
{
    /*

        Get the integrate struct of a given loop order from a spec info struct

    */

    _spec_info_t *info = _spec_info_get_struct(id);

    if (info -> loopOrderMax <= loopOrder)
      {
        printf("Cannot get the %ld'th loop order's 'intgrt_t' struct of a spectrum if only %ld loops are expected.\n", loopOrder, info -> loopOrderMax);
        exit(1);

        return NULL;
      }

    return info -> loopIntgrt[loopOrder];
}


size_t spec_info_get_kern_order(const char *id)
{
    /*

        Get the kern order of a spectrum

    */

    _spec_info_t *info = _spec_info_get_struct(id);

    return info -> specOrder + 2 * info -> loopOrder - 1;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Power Spectrum  **/

/* Loop Integration Variables */
static const size_t _pnlOneLoopIntgrtDim = 3;
static const double _pnlOneLoopIntgrtUpperBounds[3] = {8., 1., 2. * M_PI};
static const double _pnlOneLoopIntgrtLowerBounds[3] = {0.0001, -1., 0.};
static const char* _pnlOneLoopIntgrtRoutine = "vegas";
static const intgrt_vegas_t _pnlOneLoopIntgrtVegas = {10000, 1000, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _pnlOneLoopIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};

/* Part Labels */
static const char *_pnlLabelTree = NULL;
static const char *_pnlLabelP22 = NULL;
static const char *_pnlLabelP13 = NULL;
static const char *_pnlLabelCtr = NULL;
static const char *_pnlLabelSn = NULL;
static const char *_pnlLabelFull = NULL;


/* Power Spectrum Derivatives */

/* Loop Integration Variables */
static const size_t _dpnlOneLoopIntgrtDim = 3;
static const double _dpnlOneLoopIntgrtUpperBounds[3] = {8., 1., 2. * M_PI};
static const double _dpnlOneLoopIntgrtLowerBounds[3] = {0.0001, -1., 0.};
static const char *_dpnlOneLoopIntgrtRoutine = "vegas";
static const intgrt_vegas_t _dpnlOneLoopIntgrtVegas = {10000, 1000, 1., 10, 0.5, 0., 0};
static const intgrt_divonne_t _dpnlOneLoopIntgrtDivonne = {1, 1, 1.e-3, 1.e-12, 0, 0, 50000, 47, 1, 1, 5, 0., 10., 0.25, 0};

/* Part Labels */
static const char *_dpnlLabelA2Ga = NULL;
static const char *_dpnlLabelD2Ga = NULL;
static const char *_dpnlLabelA3GaA = NULL;
static const char *_dpnlLabelA3GaB = NULL;
static const char *_dpnlLabelD3GaA = NULL;
static const char *_dpnlLabelD3GaB = NULL;

static const char *_dpnlLabelB1 = NULL;
static const char *_dpnlLabelB2 = NULL;
static const char *_dpnlLabelC2Ga = NULL;
static const char *_dpnlLabelBGam3 = NULL;

static const char *_dpnlLabelF = NULL;
static const char *_dpnlLabelSigv = NULL;

static const char *_dpnlLabelD = NULL;
static const char *_dpnlLabelH = NULL;

static const char *_dpnlLabelC0 = NULL;
static const char *_dpnlLabelC2 = NULL;
static const char *_dpnlLabelC4 = NULL;

static const char *_dpnlLabelPsn = NULL;
static const char *_dpnlLabelBsn1 = NULL;
static const char *_dpnlLabelBsn2 = NULL;

static const char *_dpnlLabelPk = NULL;


/* Bispectrum */

/* Loop Integration Variables */
 // None

/* Part Labels */
static const char *_btrLabelTree = NULL;
static const char *_btrLabelSn = NULL;
static const char *_btrLabelFull = NULL;


/* Bipectrum Derivatives */

/* Loop Integration Variables */
 // None

/* Part Labels */
static const char *_dbtrLabelA2Ga = NULL;
static const char *_dbtrLabelD2Ga = NULL;
static const char *_dbtrLabelA3GaA = NULL;
static const char *_dbtrLabelA3GaB = NULL;
static const char *_dbtrLabelD3GaA = NULL;
static const char *_dbtrLabelD3GaB = NULL;

static const char *_dbtrLabelB1 = NULL;
static const char *_dbtrLabelB2 = NULL;
static const char *_dbtrLabelC2Ga = NULL;
static const char *_dbtrLabelBGam3 = NULL;

static const char *_dbtrLabelF = NULL;
static const char *_dbtrLabelSigv = NULL;

static const char *_dbtrLabelD = NULL;
static const char *_dbtrLabelH = NULL;

static const char *_dbtrLabelC0 = NULL;
static const char *_dbtrLabelC2 = NULL;
static const char *_dbtrLabelC4 = NULL;

static const char *_dbtrLabelPsn = NULL;
static const char *_dbtrLabelBsn1 = NULL;
static const char *_dbtrLabelBsn2 = NULL;

static const char *_dbtrLabelPk = NULL;


/* Trispectrum */

/* Loop Integration Variables */
 // None

/* Part Labels */
static const char *_ttrLabelTree = NULL;
static const char *_ttrLabelSn = NULL;
static const char *_ttrLabelFull = NULL;


/**  Interpolation of the spectra from file  **/

typedef struct
{
    /*

        Interpolate spectra from file

    */

    /* Interpolation struct */
    interp_t *interp;

    /* Interpolation init / eval / free functions */
    interp_t *(*interpInit)(dat_t*, double (*)(void*, void*, void*, size_t*, size_t));
    double (*interpEval)(void*, interp_t*, void*);
    interp_t *(*interpFree)(interp_t*);

} _spec_interp_t;


/* Power spectrum */

/* File */
static char *_pnlFileIn;

/* Label */
static char *_pnlLabelIn;

/* Power spectrum interpolation struct */
static _spec_interp_t **_pnlInterpIn = NULL;

static _spec_interp_t _pnlInterpInPnl;


/* Power spectrum derivatives */

/* File */
static char *_dpnlFileIn;

/* Power spectrum derivatives interpolation struct */
static _spec_interp_t **_dpnlInterpIn = NULL;

static _spec_interp_t _dpnlInterpInA2Ga;
static _spec_interp_t _dpnlInterpInD2Ga;
static _spec_interp_t _dpnlInterpInA3GaA;
static _spec_interp_t _dpnlInterpInA3GaB;
static _spec_interp_t _dpnlInterpInD3GaA;
static _spec_interp_t _dpnlInterpInD3GaB;

static _spec_interp_t _dpnlInterpInB1;
static _spec_interp_t _dpnlInterpInB2;
static _spec_interp_t _dpnlInterpInC2Ga;
static _spec_interp_t _dpnlInterpInBGam3;

static _spec_interp_t _dpnlInterpInF;
static _spec_interp_t _dpnlInterpInSigv;

static _spec_interp_t _dpnlInterpInD;
static _spec_interp_t _dpnlInterpInH;

static _spec_interp_t _dpnlInterpInC0;
static _spec_interp_t _dpnlInterpInC2;
static _spec_interp_t _dpnlInterpInC4;

static _spec_interp_t _dpnlInterpInPsn;
static _spec_interp_t _dpnlInterpInBsn1;
static _spec_interp_t _dpnlInterpInBsn2;

static _spec_interp_t _dpnlInterpInPk;


/**  Integral bins for derivatives  **/

static double **_pnl1LoopBin;

static double **_dpnlK1LoopBin;
static double **_dpnlMu1LoopBin;



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Initialise Local Variables   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_labels(void);
static int _ini_info(void);
static int _ini_bin(void);
static int _ini_interp(void);
static int _ini_func(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_ini(void)
{
    /*


        Initialise spectra modules

    */

    /* Initialise labels */
    _ini_labels();

    /* Initialise info */
    _ini_info();

    /* Initialise integration bins */
    _ini_bin();

    /* Initialise the interpolation of spectra */
    _ini_interp();

    /* Initialise spectra functions */
    _ini_func();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------    Initialise Labels    ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_labels_pnl(void);
static int _ini_labels_dpnl(void);

static int _ini_labels_btr(void);
static int _ini_labels_dbtr(void);

static int _ini_labels_ttr(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_labels(void)
{
    /*

        Initialise the labels

    */

    /* Non-Linear Power Spectrum */
    _ini_labels_pnl();

    /* Non-Linear Power Spectrum Derivatives */
    _ini_labels_dpnl();


    /* Tree-Level Bispectrum */
    _ini_labels_btr();

    /* Tree-Level Bispectrum Derivatives */
    _ini_labels_dbtr();


    /* Tree-Level Trispectrum */
    _ini_labels_ttr();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_labels_pnl(void)
{
    /*

        Initialise the non-linear power spectrum labels

    */

    _pnlLabelTree = misc_scat(5, "'", _idSpecPnl_, "-", _idTree_, "'");
    _pnlLabelP22 = misc_scat(5, "'", _idSpecPnl_, "-", "p22", "'");
    _pnlLabelP13 = misc_scat(5, "'", _idSpecPnl_, "-", "p13", "'");
    _pnlLabelSn = misc_scat(5, "'", _idSpecPnl_, "-", _idSn_, "'");
    _pnlLabelCtr = misc_scat(5, "'", _idSpecPnl_, "-", _idCtr_, "'");
    _pnlLabelFull = misc_scat(5, "'", _idSpecPnl_, "-", _idFull_, "'");

    return 0;
}


static int _ini_labels_dpnl(void)
{
    /*

        Initialise the non-linear power spectrum derivatives labels

    */

    _dpnlLabelA2Ga = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsA2Ga_, "'");
    _dpnlLabelD2Ga = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsD2Ga_, "'");
    _dpnlLabelA3GaA = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsA3GaA_, "'");
    _dpnlLabelA3GaB = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsA3GaB_, "'");
    _dpnlLabelD3GaA = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsD3GaA_, "'");
    _dpnlLabelD3GaB = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsD3GaB_, "'");

    _dpnlLabelB1 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsB1_, "'");
    _dpnlLabelB2 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsB2_, "'");
    _dpnlLabelC2Ga = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsC2Ga_, "'");
    _dpnlLabelBGam3 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsBGam3_, "'");

    _dpnlLabelF = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsF_, "'");
    _dpnlLabelSigv = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsSigv_, "'");

    _dpnlLabelD = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsD_, "'");
    _dpnlLabelH = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsH_, "'");

    _dpnlLabelC0 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsC0_, "'");
    _dpnlLabelC2 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsC2_, "'");
    _dpnlLabelC4 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsC4_, "'");

    _dpnlLabelPsn = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsPsn_, "'");
    _dpnlLabelBsn1 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsBsn1_, "'");
    _dpnlLabelBsn2 = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsBsn2_, "'");

    _dpnlLabelPk = misc_scat(5, "'", _idSpecDPnl_, "-", _idParamsPk_, "'");

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_labels_btr(void)
{
    /*

        Initialise the tree-level bispectrum

    */

    _btrLabelTree = misc_scat(5, "'", _idSpecBtr_, "-", _idTree_, "'");
    _btrLabelSn = misc_scat(5, "'", _idSpecBtr_, "-", _idSn_, "'");
    _btrLabelFull = misc_scat(5, "'", _idSpecBtr_, "-", _idFull_, "'");

    return 0;
}


static int _ini_labels_dbtr(void)
{
    /*

        Initialise the tree-level bispectrum derivatives labels

    */

    _dbtrLabelA2Ga = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsA2Ga_, "'");
    _dbtrLabelD2Ga = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsD2Ga_, "'");
    _dbtrLabelA3GaA = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsA3GaA_, "'");
    _dbtrLabelA3GaB = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsA3GaB_, "'");
    _dbtrLabelD3GaA = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsD3GaA_, "'");
    _dbtrLabelD3GaB = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsD3GaB_, "'");

    _dbtrLabelB1 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsB1_, "'");
    _dbtrLabelB2 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsB2_, "'");
    _dbtrLabelC2Ga = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsC2Ga_, "'");
    _dbtrLabelBGam3 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsBGam3_, "'");

    _dbtrLabelF = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsF_, "'");
    _dbtrLabelSigv = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsSigv_, "'");

    _dbtrLabelD = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsD_, "'");
    _dbtrLabelH = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsH_, "'");

    _dbtrLabelC0 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsC0_, "'");
    _dbtrLabelC2 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsC2_, "'");
    _dbtrLabelC4 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsC4_, "'");

    _dbtrLabelPsn = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsPsn_, "'");
    _dbtrLabelBsn1 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsBsn1_, "'");
    _dbtrLabelBsn2 = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsBsn2_, "'");

    _dbtrLabelPk = misc_scat(5, "'", _idSpecDBtr_, "-", _idParamsPk_, "'");

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_labels_ttr(void)
{
    /*

        Initialise the tree-level trispectrum

    */

    _ttrLabelTree = misc_scat(5, "'", _idSpecTtr_, "-", _idTree_, "'");
    _ttrLabelSn = misc_scat(5, "'", _idSpecTtr_, "-", _idSn_, "'");
    _ttrLabelFull = misc_scat(5, "'", _idSpecTtr_, "-", _idFull_, "'");

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_info_pnl(void);
static int _ini_info_dpnl(void);

static int _ini_info_btr(void);
static int _ini_info_dbtr(void);

static int _ini_info_ttr(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info(void)
{
    /*

        Initialise the info structs

    */

    /* Power spectrum */
    _ini_info_pnl();

    /* Power spectrum derivatives */
    _ini_info_dpnl();


    /* Bispectrum */
    _ini_info_btr();

    /* Bispectrum derivatives */
    _ini_info_dbtr();


    /* Trispectrum */
    _ini_info_ttr();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info_pnl(void)
{
    /*

        Initialise the power spectrum info struct

    */


    _pnlInfo = malloc(sizeof(_spec_info_t));

    _pnlInfo -> id = _idSpecPnl_;

    _pnlInfo -> specOrder = 2;


    /**  Loop Order of the Power Spectrum  **/

    _pnlInfo -> loopOrder = 1;
    _pnlInfo -> loopOrderMax = 1;

    _pnlInfo -> loopIntgrt = malloc(sizeof(intgrt_t*) * _pnlInfo -> loopOrderMax);

    /* 1 - loop */
    _pnlInfo -> loopIntgrt[0] = integrate_new();

    integrate_set_dim(_pnlInfo -> loopIntgrt[0], (size_t) _pnlOneLoopIntgrtDim);
    integrate_set_bounds_upper(_pnlInfo -> loopIntgrt[0], (double*) _pnlOneLoopIntgrtUpperBounds);
    integrate_set_bounds_lower(_pnlInfo -> loopIntgrt[0], (double*) _pnlOneLoopIntgrtLowerBounds);

    integrate_set_routine(_pnlInfo -> loopIntgrt[0], (const char*) _pnlOneLoopIntgrtRoutine);

    integrate_set_vegas(_pnlInfo -> loopIntgrt[0], (intgrt_vegas_t*) &_pnlOneLoopIntgrtVegas);
    integrate_set_divonne(_pnlInfo -> loopIntgrt[0], (intgrt_divonne_t*) &_pnlOneLoopIntgrtDivonne);


    /**  Parts of the Power Spectrum  **/

    _pnlInfo -> partsSize = 6;
    _pnlInfo -> partsLabels = malloc(sizeof(char*) * _pnlInfo -> partsSize);
    _pnlInfo -> partsHasErr = malloc(sizeof(bool*) * _pnlInfo -> partsSize);
    _pnlInfo -> partsExist = malloc(sizeof(bool*) * _pnlInfo -> partsSize);

    _pnlInfo -> partsDerivMult = NULL;
    _pnlInfo -> partsDerivVals = NULL;
    _pnlInfo -> partsDerivLog = NULL;

    /* Tree level */
    _pnlInfo -> partsLabels[0] = (char*) _pnlLabelTree;
    _pnlInfo -> partsHasErr[0] = (bool*) &_false_;
    _pnlInfo -> partsExist[0] = (bool*) &_true_;

    /* One-loop 22 contribution */
    _pnlInfo -> partsLabels[1] = (char*) _pnlLabelP22;
    _pnlInfo -> partsHasErr[1] = (bool*) &_true_;
    _pnlInfo -> partsExist[1] = (bool*) &_true_;

    /* One-loop 13 + 31 contributions */
    _pnlInfo -> partsLabels[2] = (char*) _pnlLabelP13;
    _pnlInfo -> partsHasErr[2] = (bool*) &_true_;
    _pnlInfo -> partsExist[2] = (bool*) &_true_;

    /* Counter terms */
    _pnlInfo -> partsLabels[3] = (char*) _pnlLabelCtr;
    _pnlInfo -> partsHasErr[3] = (bool*) &_false_;
    _pnlInfo -> partsExist[3] = (bool*) &_fidInclCtr_;

    /* Shot noise */
    _pnlInfo -> partsLabels[4] = (char*) _pnlLabelSn;
    _pnlInfo -> partsHasErr[4] = (bool*) &_false_;
    _pnlInfo -> partsExist[4] = (bool*) &_fidInclSn_;

    /* Non-linear power spectrum */
    _pnlInfo -> partsLabels[5] = (char*) _pnlLabelFull;
    _pnlInfo -> partsHasErr[5] = (bool*) &_true_;
    _pnlInfo -> partsExist[5] = (bool*) &_true_;


    return 0;
}


static int _ini_info_dpnl(void)
{
    /*

        Initialise the power spectrum derivatives info struct

    */


    _dpnlInfo = malloc(sizeof(_spec_info_t));

    _dpnlInfo -> id = _idSpecDPnl_;

    _dpnlInfo -> specOrder = 2;


    /**  Loop Order of the Power Spectrum  **/

    _dpnlInfo -> loopOrder = 1;
    _dpnlInfo -> loopOrderMax = 1;

    _dpnlInfo -> loopIntgrt = malloc(sizeof(intgrt_t*) * _dpnlInfo -> loopOrderMax);

    /* 1 - loop */
    _dpnlInfo -> loopIntgrt[0] = integrate_new();

    integrate_set_dim(_dpnlInfo -> loopIntgrt[0], (size_t) _dpnlOneLoopIntgrtDim);
    integrate_set_bounds_upper(_dpnlInfo -> loopIntgrt[0], (double*) _dpnlOneLoopIntgrtUpperBounds);
    integrate_set_bounds_lower(_dpnlInfo -> loopIntgrt[0], (double*) _dpnlOneLoopIntgrtLowerBounds);

    integrate_set_routine(_dpnlInfo -> loopIntgrt[0], (const char*) _dpnlOneLoopIntgrtRoutine);

    integrate_set_vegas(_dpnlInfo -> loopIntgrt[0], (intgrt_vegas_t*) &_dpnlOneLoopIntgrtVegas);
    integrate_set_divonne(_dpnlInfo -> loopIntgrt[0], (intgrt_divonne_t*) &_dpnlOneLoopIntgrtDivonne);


    /**  Parts of the Power Spectrum  **/

    // TODO: Maybe have partsExist variables for each loop level...

    _dpnlInfo -> partsSize = __ID_PARAMS_SIZE__;
    _dpnlInfo -> partsLabels = malloc(sizeof(char*) * _dpnlInfo -> partsSize);
    _dpnlInfo -> partsHasErr = malloc(sizeof(bool*) * _dpnlInfo -> partsSize);
    _dpnlInfo -> partsExist = malloc(sizeof(bool*) * _dpnlInfo -> partsSize);

    _dpnlInfo -> partsDerivMult = malloc(sizeof(int*) * _dpnlInfo -> partsSize);
    _dpnlInfo -> partsDerivVals = malloc(sizeof(dat_t*) * _dpnlInfo -> partsSize);
    _dpnlInfo -> partsDerivLog = malloc(sizeof(bool) * _dpnlInfo -> partsSize);

    /* a2Ga */
    _dpnlInfo -> partsLabels[0] = (char*) _dpnlLabelA2Ga;
    _dpnlInfo -> partsHasErr[0] = (bool*) &_true_;
    _dpnlInfo -> partsExist[0] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[0] = (int*) _multParamsA2Ga_;
    _dpnlInfo -> partsDerivVals[0] = NULL;
    _dpnlInfo -> partsDerivLog[0] = _false_;

    /* da2Ga */
    _dpnlInfo -> partsLabels[1] = (char*) _dpnlLabelD2Ga;
    _dpnlInfo -> partsHasErr[1] = (bool*) &_true_;
    _dpnlInfo -> partsExist[1] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[1] = (int*) _multParamsD2Ga_;
    _dpnlInfo -> partsDerivVals[1] = NULL;
    _dpnlInfo -> partsDerivLog[1] = _false_;

    /* a3GaA */
    _dpnlInfo -> partsLabels[2] = (char*) _dpnlLabelA3GaA;
    _dpnlInfo -> partsHasErr[2] = (bool*) &_true_;
    _dpnlInfo -> partsExist[2] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[2] = (int*) _multParamsA3GaA_;
    _dpnlInfo -> partsDerivVals[2] = NULL;
    _dpnlInfo -> partsDerivLog[2] = _false_;

    /* a3GaB */
    _dpnlInfo -> partsLabels[3] = (char*) _dpnlLabelA3GaB;
    _dpnlInfo -> partsHasErr[3] = (bool*) &_true_;
    _dpnlInfo -> partsExist[3] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[3] = (int*) _multParamsA3GaB_;
    _dpnlInfo -> partsDerivVals[3] = NULL;
    _dpnlInfo -> partsDerivLog[3] = _false_;

    /* d3GaA */
    _dpnlInfo -> partsLabels[4] = (char*) _dpnlLabelD3GaA;
    _dpnlInfo -> partsHasErr[4] = (bool*) &_true_;
    _dpnlInfo -> partsExist[4] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[4] = (int*) _multParamsD3GaA_;
    _dpnlInfo -> partsDerivVals[4] = NULL;
    _dpnlInfo -> partsDerivLog[4] = _false_;

    /* d3GaB */
    _dpnlInfo -> partsLabels[5] = (char*) _dpnlLabelD3GaB;
    _dpnlInfo -> partsHasErr[5] = (bool*) &_true_;
    _dpnlInfo -> partsExist[5] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[5] = (int*) _multParamsD3GaB_;
    _dpnlInfo -> partsDerivVals[5] = NULL;
    _dpnlInfo -> partsDerivLog[5] = _false_;


    /* b1 */
    _dpnlInfo -> partsLabels[6] = (char*) _dpnlLabelB1;
    _dpnlInfo -> partsHasErr[6] = (bool*) &_true_;
    _dpnlInfo -> partsExist[6] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[6] = (int*) _multParamsB1_;
    _dpnlInfo -> partsDerivVals[6] = NULL;
    _dpnlInfo -> partsDerivLog[6] = _false_;

    /* b2 */
    _dpnlInfo -> partsLabels[7] = (char*) _dpnlLabelB2;
    _dpnlInfo -> partsHasErr[7] = (bool*) &_true_;
    _dpnlInfo -> partsExist[7] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[7] = (int*) _multParamsB2_;
    _dpnlInfo -> partsDerivVals[7] = NULL;
    _dpnlInfo -> partsDerivLog[7] = _false_;

    /* c2Ga */
    _dpnlInfo -> partsLabels[8] = (char*) _dpnlLabelC2Ga;
    _dpnlInfo -> partsHasErr[8] = (bool*) &_true_;
    _dpnlInfo -> partsExist[8] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[8] = (int*) _multParamsC2Ga_;
    _dpnlInfo -> partsDerivVals[8] = NULL;
    _dpnlInfo -> partsDerivLog[8] = _false_;

    /* bGam3 */
    _dpnlInfo -> partsLabels[9] = (char*) _dpnlLabelBGam3;
    _dpnlInfo -> partsHasErr[9] = (bool*) &_true_;
    _dpnlInfo -> partsExist[9] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[9] = (int*) _multParamsBGam3_;
    _dpnlInfo -> partsDerivVals[9] = NULL;
    _dpnlInfo -> partsDerivLog[9] = _false_;


    /* f */
    _dpnlInfo -> partsLabels[10] = (char*) _dpnlLabelF;
    _dpnlInfo -> partsHasErr[10] = (bool*) &_true_;
    _dpnlInfo -> partsExist[10] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[10] = (int*) _multParamsF_;
    _dpnlInfo -> partsDerivVals[10] = NULL;
    _dpnlInfo -> partsDerivLog[10] = _false_;

    /* sigv */
    _dpnlInfo -> partsLabels[11] = (char*) _dpnlLabelSigv;
    _dpnlInfo -> partsHasErr[11] = (bool*) &_true_;
    _dpnlInfo -> partsExist[11] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[11] = (int*) _multParamsSigv_;
    _dpnlInfo -> partsDerivVals[11] = NULL;
    _dpnlInfo -> partsDerivLog[11] = _false_;


    /* d */
    _dpnlInfo -> partsLabels[12] = (char*) _dpnlLabelD;
    _dpnlInfo -> partsHasErr[12] = (bool*) &_true_;
    _dpnlInfo -> partsExist[12] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[12] = (int*) _multParamsD_;
    _dpnlInfo -> partsDerivVals[12] = NULL;
    _dpnlInfo -> partsDerivLog[12] = _false_;

    /* h */
    _dpnlInfo -> partsLabels[13] = (char*) _dpnlLabelH;
    _dpnlInfo -> partsHasErr[13] = (bool*) &_true_;
    _dpnlInfo -> partsExist[13] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[13] = (int*) _multParamsH_;
    _dpnlInfo -> partsDerivVals[13] = NULL;
    _dpnlInfo -> partsDerivLog[13] = _false_;


    /* c0 */
    _dpnlInfo -> partsLabels[14] = (char*) _dpnlLabelC0;
    _dpnlInfo -> partsHasErr[14] = (bool*) &_false_;
    _dpnlInfo -> partsExist[14] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[14] = (int*) _multParamsC0_;
    _dpnlInfo -> partsDerivVals[14] = NULL;
    _dpnlInfo -> partsDerivLog[14] = _false_;

    /* c2 */
    _dpnlInfo -> partsLabels[15] = (char*) _dpnlLabelC2;
    _dpnlInfo -> partsHasErr[15] = (bool*) &_false_;
    _dpnlInfo -> partsExist[15] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[15] = (int*) _multParamsC2_;
    _dpnlInfo -> partsDerivVals[15] = NULL;
    _dpnlInfo -> partsDerivLog[15] = _false_;

    /* c4 */
    _dpnlInfo -> partsLabels[16] = (char*) _dpnlLabelC4;
    _dpnlInfo -> partsHasErr[16] = (bool*) &_false_;
    _dpnlInfo -> partsExist[16] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[16] = (int*) _multParamsC4_;
    _dpnlInfo -> partsDerivVals[16] = NULL;
    _dpnlInfo -> partsDerivLog[16] = _false_;


    /* Psn */
    _dpnlInfo -> partsLabels[17] = (char*) _dpnlLabelPsn;
    _dpnlInfo -> partsHasErr[17] = (bool*) &_false_;
    _dpnlInfo -> partsExist[17] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[17] = (int*) _multParamsPsn_;
    _dpnlInfo -> partsDerivVals[17] = NULL;
    _dpnlInfo -> partsDerivLog[17] = _false_;

    /* Bsn1 */
    _dpnlInfo -> partsLabels[18] = (char*) _dpnlLabelBsn1;
    _dpnlInfo -> partsHasErr[18] = (bool*) &_false_;
    _dpnlInfo -> partsExist[18] = (bool*) &_false_;

    _dpnlInfo -> partsDerivMult[18] = (int*) _multParamsBsn1_;
    _dpnlInfo -> partsDerivVals[18] = NULL;
    _dpnlInfo -> partsDerivLog[18] = _false_;

    /* Bsn2 */
    _dpnlInfo -> partsLabels[19] = (char*) _dpnlLabelBsn2;
    _dpnlInfo -> partsHasErr[19] = (bool*) &_false_;
    _dpnlInfo -> partsExist[19] = (bool*) &_false_;

    _dpnlInfo -> partsDerivMult[19] = (int*) _multParamsBsn2_;
    _dpnlInfo -> partsDerivVals[19] = NULL;
    _dpnlInfo -> partsDerivLog[19] = _false_;


    /* Pk */
    _dpnlInfo -> partsLabels[20] = (char*) _dpnlLabelPk;
    _dpnlInfo -> partsHasErr[20] = (bool*) &_true_;
    _dpnlInfo -> partsExist[20] = (bool*) &_true_;

    _dpnlInfo -> partsDerivMult[20] = (int*) _multParamsPk_;
    _dpnlInfo -> partsDerivVals[20] = NULL;
    _dpnlInfo -> partsDerivLog[20] = _false_;


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info_btr(void)
{
    /*

        Initialise the bispectrum info struct

    */


    _btrInfo = malloc(sizeof(_spec_info_t));

    _btrInfo -> id = _idSpecBtr_;

    _btrInfo -> specOrder = 3;


    /**  Loop Order of the Bispectrum  **/

    _btrInfo -> loopOrder = 0;
    _btrInfo -> loopOrderMax = 0;

    _btrInfo -> loopIntgrt = malloc(sizeof(intgrt_t*) * _btrInfo -> loopOrderMax);


    /**  Parts of the Bispectrum  **/

    _btrInfo -> partsSize = 3;
    _btrInfo -> partsLabels = malloc(sizeof(char*) * _btrInfo -> partsSize);
    _btrInfo -> partsHasErr = malloc(sizeof(bool*) * _btrInfo -> partsSize);
    _btrInfo -> partsExist = malloc(sizeof(bool*) * _btrInfo -> partsSize);

    _btrInfo -> partsDerivMult = NULL;
    _btrInfo -> partsDerivVals = NULL;
    _btrInfo -> partsDerivLog = NULL;

    /* Tree level */
    _btrInfo -> partsLabels[0] = (char*) _btrLabelTree;
    _btrInfo -> partsHasErr[0] = (bool*) &_false_;
    _btrInfo -> partsExist[0] = (bool*) &_true_;

    /* Shot noise */
    _btrInfo -> partsLabels[1] = (char*) _btrLabelSn;
    _btrInfo -> partsHasErr[1] = (bool*) &_false_;
    _btrInfo -> partsExist[1] = (bool*) &_fidInclSn_;

    /* Non-linear power spectrum */
    _btrInfo -> partsLabels[2] = (char*) _btrLabelFull;
    _btrInfo -> partsHasErr[2] = (bool*) &_false_;
    _btrInfo -> partsExist[2] = (bool*) &_true_;


    return 0;
}


static int _ini_info_dbtr(void)
{
    /*

        Initialise the bispectrum derivatives info struct

    */


    _dbtrInfo = malloc(sizeof(_spec_info_t));

    _dbtrInfo -> id = _idSpecDBtr_;

    _dbtrInfo -> specOrder = 3;


    /**  Loop Order of the Bispectrum  **/

    _dbtrInfo -> loopOrder = 0;
    _dbtrInfo -> loopOrderMax = 0;

    _dbtrInfo -> loopIntgrt = malloc(sizeof(intgrt_t*) * _dbtrInfo -> loopOrderMax);


    /**  Parts of the Bispectrum  **/

    _dbtrInfo -> partsSize = __ID_PARAMS_SIZE__;
    _dbtrInfo -> partsLabels = malloc(sizeof(char*) * _dbtrInfo -> partsSize);
    _dbtrInfo -> partsHasErr = malloc(sizeof(bool*) * _dbtrInfo -> partsSize);
    _dbtrInfo -> partsExist = malloc(sizeof(bool*) * _dbtrInfo -> partsSize);

    _dbtrInfo -> partsDerivMult = malloc(sizeof(int*) * _dbtrInfo -> partsSize);
    _dbtrInfo -> partsDerivVals = malloc(sizeof(dat_t*) * _dbtrInfo -> partsSize);
    _dbtrInfo -> partsDerivLog = malloc(sizeof(bool) * _dbtrInfo -> partsSize);

    /* a2Ga */
    _dbtrInfo -> partsLabels[0] = (char*) _dbtrLabelA2Ga;
    _dbtrInfo -> partsHasErr[0] = (bool*) &_false_;
    _dbtrInfo -> partsExist[0] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[0] = (int*) _multParamsA2Ga_;
    _dbtrInfo -> partsDerivVals[0] = NULL;
    _dbtrInfo -> partsDerivLog[0] = _false_;

    /* da2Ga */
    _dbtrInfo -> partsLabels[1] = (char*) _dbtrLabelD2Ga;
    _dbtrInfo -> partsHasErr[1] = (bool*) &_false_;
    _dbtrInfo -> partsExist[1] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[1] = (int*) _multParamsD2Ga_;
    _dbtrInfo -> partsDerivVals[1] = NULL;
    _dbtrInfo -> partsDerivLog[1] = _false_;

    /* a3GaA */
    _dbtrInfo -> partsLabels[2] = (char*) _dbtrLabelA3GaA;
    _dbtrInfo -> partsHasErr[2] = (bool*) &_false_;
    _dbtrInfo -> partsExist[2] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[2] = (int*) _multParamsA3GaA_;
    _dbtrInfo -> partsDerivVals[2] = NULL;
    _dbtrInfo -> partsDerivLog[2] = _false_;

    /* a3GaB */
    _dbtrInfo -> partsLabels[3] = (char*) _dbtrLabelA3GaB;
    _dbtrInfo -> partsHasErr[3] = (bool*) &_false_;
    _dbtrInfo -> partsExist[3] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[3] = (int*) _multParamsA3GaB_;
    _dbtrInfo -> partsDerivVals[3] = NULL;
    _dbtrInfo -> partsDerivLog[3] = _false_;

    /* d3GaA */
    _dbtrInfo -> partsLabels[4] = (char*) _dbtrLabelD3GaA;
    _dbtrInfo -> partsHasErr[4] = (bool*) &_false_;
    _dbtrInfo -> partsExist[4] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[4] = (int*) _multParamsD3GaA_;
    _dbtrInfo -> partsDerivVals[4] = NULL;
    _dbtrInfo -> partsDerivLog[4] = _false_;

    /* d3GaB */
    _dbtrInfo -> partsLabels[5] = (char*) _dbtrLabelD3GaB;
    _dbtrInfo -> partsHasErr[5] = (bool*) &_false_;
    _dbtrInfo -> partsExist[5] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[5] = (int*) _multParamsD3GaB_;
    _dbtrInfo -> partsDerivVals[5] = NULL;
    _dbtrInfo -> partsDerivLog[5] = _false_;


    /* b1 */
    _dbtrInfo -> partsLabels[6] = (char*) _dbtrLabelB1;
    _dbtrInfo -> partsHasErr[6] = (bool*) &_false_;
    _dbtrInfo -> partsExist[6] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[6] = (int*) _multParamsB1_;
    _dbtrInfo -> partsDerivVals[6] = NULL;
    _dbtrInfo -> partsDerivLog[6] = _false_;

    /* b2 */
    _dbtrInfo -> partsLabels[7] = (char*) _dbtrLabelB2;
    _dbtrInfo -> partsHasErr[7] = (bool*) &_false_;
    _dbtrInfo -> partsExist[7] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[7] = (int*) _multParamsB2_;
    _dbtrInfo -> partsDerivVals[7] = NULL;
    _dbtrInfo -> partsDerivLog[7] = _false_;

    /* c2Ga */
    _dbtrInfo -> partsLabels[8] = (char*) _dbtrLabelC2Ga;
    _dbtrInfo -> partsHasErr[8] = (bool*) &_false_;
    _dbtrInfo -> partsExist[8] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[8] = (int*) _multParamsC2Ga_;
    _dbtrInfo -> partsDerivVals[8] = NULL;
    _dbtrInfo -> partsDerivLog[8] = _false_;

    /* bGam3 */
    _dbtrInfo -> partsLabels[9] = (char*) _dbtrLabelBGam3;
    _dbtrInfo -> partsHasErr[9] = (bool*) &_false_;
    _dbtrInfo -> partsExist[9] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[9] = (int*) _multParamsBGam3_;
    _dbtrInfo -> partsDerivVals[9] = NULL;
    _dbtrInfo -> partsDerivLog[9] = _false_;


    /* f */
    _dbtrInfo -> partsLabels[10] = (char*) _dbtrLabelF;
    _dbtrInfo -> partsHasErr[10] = (bool*) &_false_;
    _dbtrInfo -> partsExist[10] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[10] = (int*) _multParamsF_;
    _dbtrInfo -> partsDerivVals[10] = NULL;
    _dbtrInfo -> partsDerivLog[10] = _false_;

    /* sigv */
    _dbtrInfo -> partsLabels[11] = (char*) _dbtrLabelSigv;
    _dbtrInfo -> partsHasErr[11] = (bool*) &_false_;
    _dbtrInfo -> partsExist[11] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[11] = (int*) _multParamsSigv_;
    _dbtrInfo -> partsDerivVals[11] = NULL;
    _dbtrInfo -> partsDerivLog[11] = _false_;


    /* d */
    _dbtrInfo -> partsLabels[12] = (char*) _dbtrLabelD;
    _dbtrInfo -> partsHasErr[12] = (bool*) &_false_;
    _dbtrInfo -> partsExist[12] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[12] = (int*) _multParamsD_;
    _dbtrInfo -> partsDerivVals[12] = NULL;
    _dbtrInfo -> partsDerivLog[12] = _false_;

    /* h */
    _dbtrInfo -> partsLabels[13] = (char*) _dbtrLabelH;
    _dbtrInfo -> partsHasErr[13] = (bool*) &_false_;
    _dbtrInfo -> partsExist[13] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[13] = (int*) _multParamsH_;
    _dbtrInfo -> partsDerivVals[13] = NULL;
    _dbtrInfo -> partsDerivLog[13] = _false_;


    /* c0 */
    _dbtrInfo -> partsLabels[14] = (char*) _dbtrLabelC0;
    _dbtrInfo -> partsHasErr[14] = (bool*) &_false_;
    _dbtrInfo -> partsExist[14] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[14] = (int*) _multParamsC0_;
    _dbtrInfo -> partsDerivVals[14] = NULL;
    _dbtrInfo -> partsDerivLog[14] = _false_;

    /* c2 */
    _dbtrInfo -> partsLabels[15] = (char*) _dbtrLabelC2;
    _dbtrInfo -> partsHasErr[15] = (bool*) &_false_;
    _dbtrInfo -> partsExist[15] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[15] = (int*) _multParamsC2_;
    _dbtrInfo -> partsDerivVals[15] = NULL;
    _dbtrInfo -> partsDerivLog[15] = _false_;

    /* c4 */
    _dbtrInfo -> partsLabels[16] = (char*) _dbtrLabelC4;
    _dbtrInfo -> partsHasErr[16] = (bool*) &_false_;
    _dbtrInfo -> partsExist[16] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[16] = (int*) _multParamsC4_;
    _dbtrInfo -> partsDerivVals[16] = NULL;
    _dbtrInfo -> partsDerivLog[16] = _false_;


    /* Psn */
    _dbtrInfo -> partsLabels[17] = (char*) _dbtrLabelPsn;
    _dbtrInfo -> partsHasErr[17] = (bool*) &_false_;
    _dbtrInfo -> partsExist[17] = (bool*) &_false_;

    _dbtrInfo -> partsDerivMult[17] = (int*) _multParamsPsn_;
    _dbtrInfo -> partsDerivVals[17] = NULL;
    _dbtrInfo -> partsDerivLog[17] = _false_;

    /* Bsn1 */
    _dbtrInfo -> partsLabels[18] = (char*) _dbtrLabelBsn1;
    _dbtrInfo -> partsHasErr[18] = (bool*) &_false_;
    _dbtrInfo -> partsExist[18] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[18] = (int*) _multParamsBsn1_;
    _dbtrInfo -> partsDerivVals[18] = NULL;
    _dbtrInfo -> partsDerivLog[18] = _false_;

    /* Bsn2 */
    _dbtrInfo -> partsLabels[19] = (char*) _dbtrLabelBsn2;
    _dbtrInfo -> partsHasErr[19] = (bool*) &_false_;
    _dbtrInfo -> partsExist[19] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[19] = (int*) _multParamsBsn2_;
    _dbtrInfo -> partsDerivVals[19] = NULL;
    _dbtrInfo -> partsDerivLog[19] = _false_;


    /* Pk */
    _dbtrInfo -> partsLabels[20] = (char*) _dbtrLabelPk;
    _dbtrInfo -> partsHasErr[20] = (bool*) &_false_;
    _dbtrInfo -> partsExist[20] = (bool*) &_true_;

    _dbtrInfo -> partsDerivMult[20] = (int*) _multParamsPk_;
    _dbtrInfo -> partsDerivVals[20] = NULL;
    _dbtrInfo -> partsDerivLog[20] = _false_;


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info_ttr(void)
{
    /*

        Initialise the trispectrum info struct

    */


    _ttrInfo = malloc(sizeof(_spec_info_t));

    _ttrInfo -> id = _idSpecBtr_;

    _ttrInfo -> specOrder = 4;


    /**  Loop Order of the Trispectrum  **/

    _ttrInfo -> loopOrder = 0;
    _ttrInfo -> loopOrderMax = 0;

    _ttrInfo -> loopIntgrt = malloc(sizeof(intgrt_t*) * _ttrInfo -> loopOrderMax);


    /**  Parts of the Trispectrum  **/

    _ttrInfo -> partsSize = 3;
    _ttrInfo -> partsLabels = malloc(sizeof(char*) * _ttrInfo -> partsSize);
    _ttrInfo -> partsHasErr = malloc(sizeof(bool*) * _ttrInfo -> partsSize);
    _ttrInfo -> partsExist = malloc(sizeof(bool*) * _ttrInfo -> partsSize);

    _ttrInfo -> partsDerivMult = NULL;
    _ttrInfo -> partsDerivVals = NULL;
    _ttrInfo -> partsDerivLog = NULL;

    /* Tree level */
    _ttrInfo -> partsLabels[0] = (char*) _ttrLabelTree;
    _ttrInfo -> partsHasErr[0] = (bool*) &_false_;
    _ttrInfo -> partsExist[0] = (bool*) &_true_;

    /* Shot noise */
    _ttrInfo -> partsLabels[1] = (char*) _ttrLabelSn;
    _ttrInfo -> partsHasErr[1] = (bool*) &_false_;
    _ttrInfo -> partsExist[1] = (bool*) &_fidInclSn_;

    /* Non-linear power spectrum */
    _ttrInfo -> partsLabels[2] = (char*) _ttrLabelFull;
    _ttrInfo -> partsHasErr[2] = (bool*) &_false_;
    _ttrInfo -> partsExist[2] = (bool*) &_true_;


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_bin(void)
{
    /*

        Initialise integration bins for the derivatives

    */

    _pnl1LoopBin = NULL;
    _dpnlK1LoopBin = NULL;
    _dpnlMu1LoopBin = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_interp_pnl(void);
static int _ini_interp_dpnl(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_interp(void)
{
    /*

        Initialise the interpolation parameters for the spectra

    */

    /* Power spectrum */
    _ini_interp_pnl();

    /* Power spectrum derivatives */
    _ini_interp_dpnl();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_interp_pnl(void)
{
    /*

        Initialise the interpolation of the power spectrum

    */

    /* File */
    _pnlFileIn = misc_scat(2, _pnlInfo -> id, _outExt);

    /* Label */
    _pnlLabelIn = misc_scat(1, _pnlLabelFull);

    /* Interpolation struct */
    _pnlInterpIn = malloc(sizeof(_spec_interp_t*));

    _pnlInterpIn[0] = &_pnlInterpInPnl;

    /* Interpolation functions */
    _pnlInterpIn[0] -> interp = NULL;

    _pnlInterpIn[0] -> interpInit = interpolate_interp_init_splinter;
    _pnlInterpIn[0] -> interpEval = interpolate_interp_eval_splinter;
    _pnlInterpIn[0] -> interpFree = interpolate_interp_free_splinter;

    return 0;
}


static int _ini_interp_dpnl(void)
{
    /*

        Initialise the interpolation of the power spectrum derivatives

    */

    /* File */
    _dpnlFileIn = misc_scat(2, _dpnlInfo -> id, _outExt);

    /* Interpolation struct */
    _dpnlInterpIn = malloc(sizeof(_spec_interp_t*) * __ID_PARAMS_SIZE__);

    _dpnlInterpIn[0] = &_dpnlInterpInA2Ga;
    _dpnlInterpIn[1] = &_dpnlInterpInD2Ga;
    _dpnlInterpIn[2] = &_dpnlInterpInA3GaA;
    _dpnlInterpIn[3] = &_dpnlInterpInA3GaB;
    _dpnlInterpIn[4] = &_dpnlInterpInD3GaA;
    _dpnlInterpIn[5] = &_dpnlInterpInD3GaB;

    _dpnlInterpIn[6] = &_dpnlInterpInB1;
    _dpnlInterpIn[7] = &_dpnlInterpInB2;
    _dpnlInterpIn[8] = &_dpnlInterpInC2Ga;
    _dpnlInterpIn[9] = &_dpnlInterpInBGam3;

    _dpnlInterpIn[10] = &_dpnlInterpInF;
    _dpnlInterpIn[11] = &_dpnlInterpInSigv;

    _dpnlInterpIn[12] = &_dpnlInterpInD;
    _dpnlInterpIn[13] = &_dpnlInterpInH;

    _dpnlInterpIn[14] = &_dpnlInterpInC0;
    _dpnlInterpIn[15] = &_dpnlInterpInC2;
    _dpnlInterpIn[16] = &_dpnlInterpInC4;

    _dpnlInterpIn[17] = &_dpnlInterpInPsn;
    _dpnlInterpIn[18] = &_dpnlInterpInBsn1;
    _dpnlInterpIn[19] = &_dpnlInterpInBsn2;

    _dpnlInterpIn[20] = &_dpnlInterpInPk;

    /* Interpolation functions */
    for (size_t i = 0; i < __ID_PARAMS_SIZE__; i++)
      {
        _dpnlInterpIn[i] -> interp = NULL;

        _dpnlInterpIn[i] -> interpInit = interpolate_interp_init_splinter;
        _dpnlInterpIn[i] -> interpEval = interpolate_interp_eval_splinter;
        _dpnlInterpIn[i] -> interpFree = interpolate_interp_free_splinter;
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_func_pnl(void);
static int _ini_func_dpnl(void);

static int _ini_func_btr(void);
static int _ini_func_dbtr(void);

static int _ini_func_ttr(void);

static int _ini_func_poly(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func(void)
{
    /*

        Initialise spectra functions

    */

    /* Power Spectrum */
    _ini_func_pnl();

    /* Power Spectrum Derivatives */
    _ini_func_dpnl();


    /* Bispectrum */
    _ini_func_btr();

    /* Bipectrum Derivatives */
    _ini_func_dbtr();


    /* Trispectrum */
    _ini_func_ttr();


    /* Poly Spectrum */
    _ini_func_poly();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _spec_ptr(void *var, void *params);
static double _spec_pnl(void *var, void *params);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func_pnl(void)
{
    /*

        Initialise the power spectrum functions

    */

    /* Tree-level power spectrum (not obtained via interpolation) */
    _specPtr_ = _spec_ptr;

    /* Non-linear power spectrum (obtained via interpolation) */
    _specPnl_ = _spec_pnl;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _spec_dpnl_a2ga(void *var, void *params);
static double _spec_dpnl_d2ga(void *var, void *params);
static double _spec_dpnl_a3gaa(void *var, void *params);
static double _spec_dpnl_a3gab(void *var, void *params);
static double _spec_dpnl_d3gaa(void *var, void *params);
static double _spec_dpnl_d3gab(void *var, void *params);

static double _spec_dpnl_b1(void *var, void *params);
static double _spec_dpnl_b2(void *var, void *params);
static double _spec_dpnl_c2ga(void *var, void *params);
static double _spec_dpnl_bgam3(void *var, void *params);

static double _spec_dpnl_f(void *var, void *params);
static double _spec_dpnl_sigv(void *var, void *params);

static double _spec_dpnl_d(void *var, void *params);
static double _spec_dpnl_h(void *var, void *params);

static double _spec_dpnl_c0(void *var, void *params);
static double _spec_dpnl_c2(void *var, void *params);
static double _spec_dpnl_c4(void *var, void *params);

static double _spec_dpnl_psn(void *var, void *params);
static double _spec_dpnl_bsn1(void *var, void *params);
static double _spec_dpnl_bsn2(void *var, void *params);

static double _spec_dpnl_pk(void *var, void *params);

static double (*_spec_dpnl(const char *label))(void*, void*);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func_dpnl(void)
{
    /*

        Initialise the power specturm derivatives functions

    */

    /* Bootstrap parameters */
    _specDPnlA2Ga_ = _spec_dpnl_a2ga;
    _specDPnlD2Ga_ = _spec_dpnl_d2ga;
    _specDPnlA3GaA_ = _spec_dpnl_a3gaa;
    _specDPnlA3GaB_ = _spec_dpnl_a3gab;
    _specDPnlD3GaA_ = _spec_dpnl_d3gaa;
    _specDPnlD3GaB_ = _spec_dpnl_d3gab;

    /* Bias parameters */
    _specDPnlB1_ = _spec_dpnl_b1;
    _specDPnlB2_ = _spec_dpnl_b2;
    _specDPnlC2Ga_ = _spec_dpnl_c2ga;
    _specDPnlBGam3_ = _spec_dpnl_bgam3;

    /* Redshift space distortion parameters */
    _specDPnlF_ = _spec_dpnl_f;
    _specDPnlSigv_ = _spec_dpnl_sigv;

    /* AP parameters */
    _specDPnlD_ = _spec_dpnl_d;
    _specDPnlH_ = _spec_dpnl_h;

    /* Counter term parameters */
    _specDPnlC0_ = _spec_dpnl_c0;
    _specDPnlC2_ = _spec_dpnl_c2;
    _specDPnlC4_ = _spec_dpnl_c4;

    /* Shot noise */
    _specDPnlPsn_ = _spec_dpnl_psn;
    _specDPnlBsn1_ = _spec_dpnl_bsn1;
    _specDPnlBsn2_ = _spec_dpnl_bsn2;

    /* Linear power spectrum */
    _specDPnlPk_ = _spec_dpnl_pk;

    /* All derivatives */
    _specDPnl_ = _spec_dpnl;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _spec_btr(void *var, void *params);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func_btr(void)
{
    /*

        Initialise the bispectrum functions

    */

    /* Tree-level bispectrum */
    _specBtr_ = _spec_btr;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _spec_dbtr_a2ga(void *var, void *params);
static double _spec_dbtr_d2ga(void *var, void *params);
static double _spec_dbtr_a3gaa(void *var, void *params);
static double _spec_dbtr_a3gab(void *var, void *params);
static double _spec_dbtr_d3gaa(void *var, void *params);
static double _spec_dbtr_d3gab(void *var, void *params);

static double _spec_dbtr_b1(void *var, void *params);
static double _spec_dbtr_b2(void *var, void *params);
static double _spec_dbtr_c2ga(void *var, void *params);
static double _spec_dbtr_bgam3(void *var, void *params);

static double _spec_dbtr_f(void *var, void *params);
static double _spec_dbtr_sigv(void *var, void *params);

static double _spec_dbtr_d(void *var, void *params);
static double _spec_dbtr_h(void *var, void *params);

static double _spec_dbtr_c0(void *var, void *params);
static double _spec_dbtr_c2(void *var, void *params);
static double _spec_dbtr_c4(void *var, void *params);

static double _spec_dbtr_psn(void *var, void *params);
static double _spec_dbtr_bsn1(void *var, void *params);
static double _spec_dbtr_bsn2(void *var, void *params);

static double _spec_dbtr_pk(void *var, void *params);

static double (*_spec_dbtr(const char *label))(void*, void*);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func_dbtr(void)
{
    /*

        Initialise the bispectrum derivatives functions

    */

    /* Bootstrap parameters */
    _specDBtrA2Ga_ = _spec_dbtr_a2ga;
    _specDBtrD2Ga_ = _spec_dbtr_d2ga;
    _specDBtrA3GaA_ = _spec_dbtr_a3gaa;
    _specDBtrA3GaB_ = _spec_dbtr_a3gab;
    _specDBtrD3GaA_ = _spec_dbtr_d3gaa;
    _specDBtrD3GaB_ = _spec_dbtr_d3gab;

    /* Bias parameters */
    _specDBtrB1_ = _spec_dbtr_b1;
    _specDBtrB2_ = _spec_dbtr_b2;
    _specDBtrC2Ga_ = _spec_dbtr_c2ga;
    _specDBtrBGam3_ = _spec_dbtr_bgam3;

    /* Redshift space distortion parameters */
    _specDBtrF_ = _spec_dbtr_f;
    _specDBtrSigv_ = _spec_dbtr_sigv;

    /* AP parameters */
    _specDBtrD_ = _spec_dbtr_d;
    _specDBtrH_ = _spec_dbtr_h;

    /* Counter term parameters */
    _specDBtrC0_ = _spec_dbtr_c0;
    _specDBtrC2_ = _spec_dbtr_c2;
    _specDBtrC4_ = _spec_dbtr_c4;

    /* Shot noise */
    _specDBtrPsn_ = _spec_dbtr_psn;
    _specDBtrBsn1_ = _spec_dbtr_bsn1;
    _specDBtrBsn2_ = _spec_dbtr_bsn2;

    /* Linear power spectrum */
    _specDBtrPk_ = _spec_dbtr_pk;

    /* All derivatives */
    _specDBtr_ = _spec_dbtr;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _spec_ttr(void *var, void *params);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func_ttr(void)
{
    /*

        Initialise the trispectrum functions

    */

    /* Tree-level trispectrum */
    _specTtr_ = _spec_ttr;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double (*_spec_poly(const char *specLabel, const char *label))(void*, void*);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func_poly(void)
{
    /*

        Initialise the poly spectrum function

    */

    /* Poly Spectrum */
    _specPoly_ = _spec_poly;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Free Local Variables   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_labels(void);
static int _free_info(void);
static int _free_interp(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_free(void)
{
    /*

        Free local variables

    */

    /* Free labels */
    _free_labels();

    /* Free info */
    _free_info();

    /* Free interp */
    _free_interp();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------------    Free Labels    ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_labels_pnl(void);
static int _free_labels_dpnl(void);

static int _free_labels_btr(void);
static int _free_labels_dbtr(void);

static int _free_labels_ttr(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_labels(void)
{
    /*

        Free the labels

    */

    /* Non-Linear Power Spectrum */
    _free_labels_pnl();

    /* Non-Linear Power Spectrum Derivatives */
    _free_labels_dpnl();


    /* Tree-Level Bispectrum */
    _free_labels_btr();

    /* Tree-Level Bispectrum Derivatives */
    _free_labels_dbtr();


    /* Tree-Level Trispectrum */
    _free_labels_ttr();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_labels_pnl(void)
{
    /*

        Initialise the non-linear power spectrum labels

    */

    /* Free labels */
    free((char*) _pnlLabelTree);
    free((char*) _pnlLabelP22);
    free((char*) _pnlLabelP13);
    free((char*) _pnlLabelSn);
    free((char*) _pnlLabelCtr);
    free((char*) _pnlLabelFull);

    /* Reset labels */
    _pnlLabelTree = NULL;
    _pnlLabelP22 = NULL;
    _pnlLabelP13 = NULL;
    _pnlLabelSn = NULL;
    _pnlLabelCtr = NULL;
    _pnlLabelFull = NULL;

    return 0;
}


static int _free_labels_dpnl(void)
{
    /*

        Initialise the non-linear power spectrum derivatives labels

    */

    /* Free labels */
    free((char*) _dpnlLabelA2Ga);
    free((char*) _dpnlLabelD2Ga);
    free((char*) _dpnlLabelA3GaA);
    free((char*) _dpnlLabelA3GaB);
    free((char*) _dpnlLabelD3GaA);
    free((char*) _dpnlLabelD3GaB);

    free((char*) _dpnlLabelB1);
    free((char*) _dpnlLabelB2);
    free((char*) _dpnlLabelC2Ga);
    free((char*) _dpnlLabelBGam3);

    free((char*) _dpnlLabelF);
    free((char*) _dpnlLabelSigv);

    free((char*) _dpnlLabelD);
    free((char*) _dpnlLabelH);

    free((char*) _dpnlLabelC0);
    free((char*) _dpnlLabelC2);
    free((char*) _dpnlLabelC4);

    free((char*) _dpnlLabelPsn);
    free((char*) _dpnlLabelBsn1);
    free((char*) _dpnlLabelBsn2);

    free((char*) _dpnlLabelPk);

    /* Reset labels */
    _dpnlLabelA2Ga = NULL;
    _dpnlLabelD2Ga = NULL;
    _dpnlLabelA3GaA = NULL;
    _dpnlLabelA3GaB = NULL;
    _dpnlLabelD3GaA = NULL;
    _dpnlLabelD3GaB = NULL;

    _dpnlLabelB1 = NULL;
    _dpnlLabelB2 = NULL;
    _dpnlLabelC2Ga = NULL;
    _dpnlLabelBGam3 = NULL;

    _dpnlLabelF = NULL;
    _dpnlLabelSigv = NULL;

    _dpnlLabelD = NULL;
    _dpnlLabelH = NULL;

    _dpnlLabelC0 = NULL;
    _dpnlLabelC2 = NULL;
    _dpnlLabelC4 = NULL;

    _dpnlLabelPsn = NULL;
    _dpnlLabelBsn1 = NULL;
    _dpnlLabelBsn2 = NULL;

    _dpnlLabelPk = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _free_labels_btr(void)
{
    /*

        Initialise the tree-level bispectrum

    */

    /* Free labels */
    free((char*) _btrLabelTree);
    free((char*) _btrLabelSn);
    free((char*) _btrLabelFull);

    /* Reset labels */
    _btrLabelTree = NULL;
    _btrLabelSn = NULL;
    _btrLabelFull = NULL;

    return 0;
}


static int _free_labels_dbtr(void)
{
    /*

        Initialise the tree-level bispectrum derivatives labels

    */

    /* Free labels */
    free((char*) _dbtrLabelA2Ga);
    free((char*) _dbtrLabelD2Ga);
    free((char*) _dbtrLabelA3GaA);
    free((char*) _dbtrLabelA3GaB);
    free((char*) _dbtrLabelD3GaA);
    free((char*) _dbtrLabelD3GaB);

    free((char*) _dbtrLabelB1);
    free((char*) _dbtrLabelB2);
    free((char*) _dbtrLabelC2Ga);
    free((char*) _dbtrLabelBGam3);

    free((char*) _dbtrLabelF);
    free((char*) _dbtrLabelSigv);

    free((char*) _dbtrLabelD);
    free((char*) _dbtrLabelH);

    free((char*) _dbtrLabelC0);
    free((char*) _dbtrLabelC2);
    free((char*) _dbtrLabelC4);

    free((char*) _dbtrLabelPsn);
    free((char*) _dbtrLabelBsn1);
    free((char*) _dbtrLabelBsn2);

    free((char*) _dbtrLabelPk);

    /* Reset labels */
    _dbtrLabelA2Ga = NULL;
    _dbtrLabelD2Ga = NULL;
    _dbtrLabelA3GaA = NULL;
    _dbtrLabelA3GaB = NULL;
    _dbtrLabelD3GaA = NULL;
    _dbtrLabelD3GaB = NULL;

    _dbtrLabelB1 = NULL;
    _dbtrLabelB2 = NULL;
    _dbtrLabelC2Ga = NULL;
    _dbtrLabelBGam3 = NULL;

    _dbtrLabelF = NULL;
    _dbtrLabelSigv = NULL;

    _dbtrLabelD = NULL;
    _dbtrLabelH = NULL;

    _dbtrLabelC0 = NULL;
    _dbtrLabelC2 = NULL;
    _dbtrLabelC4 = NULL;

    _dbtrLabelPsn = NULL;
    _dbtrLabelBsn1 = NULL;
    _dbtrLabelBsn2 = NULL;

    _dbtrLabelPk = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _free_labels_ttr(void)
{
    /*

        Initialise the tree-level tr

    */

    /* Free labels */
    free((char*) _ttrLabelTree);
    free((char*) _ttrLabelSn);
    free((char*) _ttrLabelFull);

    /* Reset labels */
    _ttrLabelTree = NULL;
    _ttrLabelSn = NULL;
    _ttrLabelFull = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_info(void)
{
    /*

        Free the info structs

    */

    /* Power spectrum */
    _pnlInfo = _spec_info_free(_pnlInfo);

    /* Power spectrum derivatives */
    _dpnlInfo = _spec_info_free(_dpnlInfo);


    /* Bispectrum */
    _btrInfo = _spec_info_free(_btrInfo);

    /* Bispectrum derivatives */
    _dbtrInfo = _spec_info_free(_dbtrInfo);


    /* Bispectrum */
    _ttrInfo = _spec_info_free(_ttrInfo);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_interp_pnl(void);
static int _free_interp_dpnl(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_interp(void)
{
    /*

        Free the interpolation parameters for the spectra

    */

    /* Power spectrum */
    _free_interp_pnl();

    /* Power spectrum derivatives */
    _free_interp_dpnl();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _free_interp_pnl(void)
{
    /*

        Free the interpolation structs for the non-linear power spectrum

    */

    free(_pnlFileIn);
    free(_pnlLabelIn);

    _pnlInterpIn[0] -> interp = interpolate_interp_free(_pnlInterpIn[0] -> interpFree, _pnlInterpIn[0] -> interp);

    free(_pnlInterpIn);

    return 0;
}


static int _free_interp_dpnl(void)
{
    /*

        Free the interpolation structs for the non-linear power spectrum derivatives

    */

    free(_dpnlFileIn);

    for (size_t i = 0; i < __ID_PARAMS_SIZE__; i++)
        _dpnlInterpIn[i] -> interp = interpolate_interp_free(_dpnlInterpIn[i] -> interpFree, _dpnlInterpIn[i] -> interp);

    free(_dpnlInterpIn);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SPECTRA FROM FILE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Setters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_in_pnl_set_file(const char *file)
{
    /*

        Set the input file of the non-linear power spectrum

    */

    _pnlFileIn = realloc(_pnlFileIn, strlen(file) + 1);
    snprintf(_pnlFileIn, strlen(file) + 1, "%s", file);

    return 0;
}


int spec_in_pnl_set_label(const char *label)
{
    /*

        Set the label of the non-linear power spectrum in the input file

    */

    _pnlLabelIn = realloc(_pnlLabelIn, strlen(label) + 1);
    snprintf(_pnlLabelIn, strlen(label) + 1, "%s", label);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_in_dpnl_set_file(const char *file)
{
    /*

        Set the input file of the non-linear power spectrum derivatives

    */

    _dpnlFileIn = realloc(_dpnlFileIn, strlen(file) + 1);
    snprintf(_dpnlFileIn, strlen(file) + 1, "%s", file);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Extrapolation functions   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_pnl_extrapolate(void *var, void *interp, void *params, size_t *indices, size_t size)
{
    /*

        Extrapolate the non-linear power spectrum in k < kmin:

                    Pnl ~ Ptree for k < kmin

    */

    (void) var;

    interp_t *interpPnl = (interp_t*) interp;
    spec_arg_t *specArg = (spec_arg_t*) params;
    kern_t *kern = specArg -> kern;

    double k = kernels_qget_k(kern, 0);

    /* Can only extrapolate in k < kmin */
    if (size > 1 || *indices != 1)
      {
        printf("Can only extrapolate the power spectrum in 'k', but not in 'z' or 'mu'.\n");
        exit(1);

        return NAN;
      }

    if (k > interpPnl -> xBounds[1][1])
      {
        printf("Can only extrapolate the power spectrum in 'k < kmin', but not in 'k > kmax'.\n");
        exit(1);

        return NAN;
      }

    kernels_qset_k(kern, 0, interpPnl -> xBounds[1][0]);
    double pnlBound = _specPnl_(specArg, NULL);
    double ptrBound = _specPtr_(specArg, NULL);

    kernels_qset_k(kern, 0, k);
    double pnl = pnlBound / ptrBound * _specPtr_(specArg, NULL);

    return pnl;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Setup the Interpolation   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_in_pnl_setup(void)
{
    /*

        Interpolate the non-linear power spectrum

    */

    /* Path to the file */
    char *path = misc_scat(3, __FISHERLSSDIR__, _outDir, _pnlFileIn);

    /* Read the data */
    dat_t *datPnl = dat_input(path, _pnlLabelIn, _true_);

    /* Interpolate the data */
    _pnlInterpInPnl.interp = interpolate_interp_free(_pnlInterpInPnl.interpFree, _pnlInterpInPnl.interp);
    _pnlInterpInPnl.interp = interpolate_interp_init(_pnlInterpInPnl.interpInit, datPnl, _spec_pnl_extrapolate);

    /* Free memory */
    free(path);
    datPnl = dat_free(datPnl);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static dat_t *_in_dpnl_setup_dat(dat_t *datDPnl, char *label);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_in_dpnl_setup(void)
{
    /*

        Interpolate the non-linear power spectrum derivatives

    */

    /* Path to the file */
    char *path = misc_scat(3, __FISHERLSSDIR__, _outDir, _dpnlFileIn);

    /* Read the data (no specific label -> read all of the data) */
    dat_t *datDPnl = dat_input(path, NULL, _true_);

    /* Interpolate the data */
    for (size_t i = 0; i < _dpnlInfo -> partsSize; i++)
      {
        /* Get the data with current label and additional x dimension depending on multiplicity */
        dat_t *datDPnlSub = _in_dpnl_setup_dat(datDPnl, _dpnlInfo -> partsLabels[i]);

        /* Free any previous interpolation structs */
        _dpnlInterpIn[i] -> interp = interpolate_interp_free(_dpnlInterpIn[i] -> interpFree, _dpnlInterpIn[i] -> interp);

        /* Interpolate the data if it was found in the file */
        if (datDPnlSub != NULL)
          {
            _dpnlInterpIn[i] -> interp = interpolate_interp_init(_dpnlInterpIn[i] -> interpInit, datDPnlSub, NULL);
          }

        datDPnlSub = dat_free(datDPnlSub);
      }

    /* Free memory */
    free(path);
    datDPnl = dat_free(datDPnl);

    return 0;
}


static dat_t *_in_dpnl_setup_dat(dat_t *datDPnl, char *label)
{
    /*

        Setup the data struct for the power spectrum derivative with given label and multiplicities

    */

    dat_t *dat = NULL;

    /* Dimensions */
    size_t xDim = datDPnl -> xDim;
    size_t yDim = 0;
    size_t size = datDPnl -> size;


    /**  Check if more x dimensions must be added  **/

    bool addX = _true_;
    size_t addXDim = 0;
    size_t *addXSizes = NULL;

    for (size_t i = 0; i < datDPnl -> yDim; i++)
      {
        /* Skip errors */
        if (misc_sin(datDPnl -> yLabels[i], _idIntError_))
            continue;

        /* Correct label */
        if (misc_sin(datDPnl -> yLabels[i], label))
          {
            /* yDim can only be 1 or 0 */
            if (yDim == 0)
                yDim = 1;

            /* Copy the correct y label */
            char *yLabel = misc_scat(1, datDPnl -> yLabels[i]);
            misc_ssrm(yLabel, label);

            /* Flag whether a double was found in the label or not */
            bool dInStr;

            /* Index of the xSizes array */
            size_t ind = 0;

            /* Check if there are any doubles in yLabel */
            do
              {
                /* Get a potential double and update yLabel */
                misc_stod(yLabel, &yLabel, &dInStr);

                /* Double found */
                if (dInStr)
                  {
                    /* Add another x dimension */
                    if (addX)
                      {
                        addXDim++;
                        addXSizes = realloc(addXSizes, sizeof(size_t) * addXDim);
                        addXSizes[addXDim - 1] = 0;
                      }

                    /* Update the xSize */
                    addXSizes[ind] += 1;
                    ind += 1;
                  }
              }

            while (dInStr);

            /* Every correct label has the same amount of doubles in it (e.g. Pk has k or f has z...) */
            addX = _false_;

            /* Free memory */
            free(yLabel);
          }
      }

    /* Add the new x dimensions */
    if (addXDim != 0)
      {
        xDim += addXDim;

        for (size_t i = 1; i < addXDim; i++)
          {
            addXSizes[0] *= addXSizes[i];
          }

        size *= addXSizes[0];
      }

    /* Free memory */
    free(addXSizes);

    /* Label not found */
    if (yDim == 0)
      {
        return dat;
      }


    /* Allocate memory to dat */
    dat = dat_new(xDim, yDim, size);


    /**  Insert the data into dat struct  **/

    size_t col = 0;
    size_t row = 0;

    for (size_t i = 0; i < datDPnl -> yDim; i++)
      {
        /* Skip errors */
        if (misc_sin(datDPnl -> yLabels[i], _idIntError_))
            continue;

        /* Correct label */
        if (misc_sin(datDPnl -> yLabels[i], label))
          {
            /* Copy the correct y label */
            char *yLabel = misc_scat(1, datDPnl -> yLabels[i]);
            misc_ssrm(yLabel, label);

            /* Get the extra x data */
            double *nums = malloc(sizeof(double) * (dat -> xDim - datDPnl -> xDim));

            for (size_t j = 0; j < dat -> xDim - datDPnl -> xDim; j++)
                nums[j] = misc_stod(yLabel, &yLabel, NULL);

            /* Copy the data from datDPnl */
            for (size_t n = 0; n < datDPnl -> size; n++)
              {
                col = 0;

                /* x data */
                for (size_t j = 0; j < datDPnl -> xDim; j++)
                  {
                    dat_set_value(dat, col, row, 'x', dat_get_value(datDPnl, j, n, 'x'));

                    col++;
                  }

                /* Extra x data */
                for (size_t j = 0; j < dat -> xDim - datDPnl -> xDim; j++)
                  {
                    dat_set_value(dat, col, row, 'x', nums[j]);

                    col++;
                  }

                /* y data */
                dat_set_value(dat, 0, row, 'y', dat_get_value(datDPnl, i, n, 'y'));
                row++;
              }

            /* Free memory */
            free(yLabel);
            free(nums);
          }
      }

    return dat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SPECTRA STRUCTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Derivative Struct   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_deriv_t *spec_deriv_new(void)
{
    /*

        Create a new derivative struct.

    */

    spec_deriv_t *deriv = malloc(sizeof(spec_deriv_t));

    /* Some initial values */
    deriv -> z = 0.;
    deriv -> k = 0.;
    deriv -> mu = 0.;

    deriv -> log = _false_;

    return deriv;
}


spec_deriv_t *spec_deriv_free(spec_deriv_t *deriv)
{
    /*

        Free deriv

    */

    free(deriv);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_deriv_set_var(spec_deriv_t *deriv, const char *label, double value)
{
    /*

        Set the variable of given label for the deriv

        TODO: Maybe check for validity of the value... (ie k in (0, inf), mu and nu in (-1, 1))

    */

    /* Redshift z */
    if (!strcmp(label, _idVarZ_))
      {
        deriv -> z = value;

        return 0;
      }

    /* Mode k */
    if (!strcmp(label, _idVarK_))
      {
        deriv -> k = value;

        return 0;
      }

    /* Angle mu */
    if (!strcmp(label, _idVarMu_))
      {
        deriv -> mu = value;

        return 0;
      }

    /* Label not found */
    if (deriv -> id != NULL)
        printf("Can only set the variable of the derivative parameter '%s' if its label is '%s', '%s' or '%s'!\n", deriv -> id, _idVarZ_, _idVarK_, _idVarMu_);

    else
        printf("Can only set the variable of a derivative parameter if its label is '%s', '%s' or '%s'!\n", _idVarZ_, _idVarK_, _idVarMu_);

    exit(1);

    return 1;
}


int spec_deriv_set_z(spec_deriv_t *deriv, double z)
{
    /*

        Set the redshift z for the deriv struct

    */

    deriv -> z = z;

    return 0;
}


int spec_deriv_set_k(spec_deriv_t *deriv, double k)
{
    /*

        Set the mode k for the deriv struct

    */

    deriv -> k = k;

    return 0;
}


int spec_deriv_set_mu(spec_deriv_t *deriv, double mu)
{
    /*

        Set the angle mu for the deriv struct

    */

    deriv -> mu = mu;

    return 0;
}


int spec_deriv_set_log(spec_deriv_t *deriv, bool log)
{
    /*

        Set the flag for calculating the logarithmic derivative for the deriv struct

    */

    deriv -> log = log;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Spectra Arguments Struct   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_arg_t *spec_arg_new(const char *id)
{
    /*

        Create a new spec_arg_t struct

    */

    spec_arg_t *specArg = malloc(sizeof(spec_arg_t));

    /* Do not allocate any memory if id is NULL (eg use this if kern should be created by hand) */
    if (id == NULL)
      {
        specArg -> id = NULL;
        specArg -> kern = NULL;
        specArg -> deriv = NULL;

        specArg -> dz = 0.;
        specArg -> dk = 0.;
        specArg -> dmu = 0.;

        specArg -> binVolume = 1.;

        return specArg;
      }

    _spec_info_t *info = _spec_info_get_struct(id);

    specArg -> id = info -> id;

    specArg -> kern = kernels_new(info -> specOrder, info -> loopOrder);
    specArg -> deriv = (spec_info_isderiv(info -> id)) ? spec_deriv_new() : NULL;

    specArg -> dz = 0.;
    specArg -> dk = 0.;
    specArg -> dmu = 0.;

    specArg -> binVolume = 1.;

    return specArg;
}


spec_arg_t *spec_arg_free(spec_arg_t *specArg)
{
    /*

        Free specArg

    */

    /* Check for NULL */
    if (specArg == NULL)
        return NULL;

    /* Free contents */
    specArg -> kern = kernels_free(specArg -> kern);
    specArg -> deriv = spec_deriv_free(specArg -> deriv);

    /* Free specArg itself */
    free(specArg);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_arg_set_dz(spec_arg_t *specArg, double dz)
{
    /*

        Set the redshift bin width dz in specArg

    */

    specArg -> dz = dz;

    return 0;
}


int spec_arg_set_dk(spec_arg_t *specArg, double dk)
{
    /*

        Set the mode bin width dk in specArg

    */

    specArg -> dk = dk;

    return 0;
}


int spec_arg_set_dmu(spec_arg_t *specArg, double dmu)
{
    /*

        Set the angle bin width dmu in specArg

    */

    specArg -> dmu = dmu;

    return 0;
}


int spec_arg_set_volume(spec_arg_t *specArg, double binVolume)
{
    /*

        Set the bin volume in specArg

    */

    specArg -> binVolume = binVolume;

    return 0;
}


int spec_arg_set_kern(spec_arg_t *specArg, kern_t *kern)
{
    /*

        Set the kern struct for specArg iff specArg -> id is NULL

    */

    if (specArg -> id != NULL)
      {
        printf("Cannot set the 'kern_t' struct of a 'spec_arg_t' struct if the ID is not NULL.\n");
        exit(1);

        return 1;
      }

    specArg -> kern = kern;

    return 0;
}


int spec_arg_set_deriv(spec_arg_t *specArg, spec_deriv_t *deriv)
{
    /*

        Set the deriv struct for specArg iff specArg -> id is NULL

    */

    if (specArg -> id != NULL)
      {
        printf("Cannot set the 'spec_deriv_t' struct of a 'spec_arg_t' struct if the ID is not NULL.\n");
        exit(1);

        return 1;
      }

    specArg -> deriv = deriv;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


spec_deriv_t *spec_arg_get_deriv(spec_arg_t *specArg)
{
    /*

        Get deriv struct from specArg

    */

    return specArg -> deriv;
}


kern_t *spec_arg_get_kern(spec_arg_t *specArg)
{
    /*

        Get the kernel struct from specArg

    */

    return specArg -> kern;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Spectra Output Struct   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_out_t *spec_out_new(void)
{
    /*

        Create a new spec_out_t struct

    */

    spec_out_t *out = malloc(sizeof(spec_out_t));

    out -> size = 0;
    out -> labels = NULL;

    out -> err = NULL;
    out -> log = NULL;

    out -> file = NULL;

    out -> binary = _true_;
    out -> precision = 0;

    return out;
}


spec_out_t *spec_out_free(spec_out_t *out)
{
    /*

        Free out

    */

    /* Check for NULL */
    if (out == NULL)
        return NULL;

    /* Free contents */
    if (out -> labels != NULL)
      {
        for (size_t i = 0; i < out ->size; i++)
            free(out -> labels[i]);
      }

    free(out -> labels);

    free(out -> err);
    free(out -> log);

    free(out -> file);

    /* Free out itself */
    free(out);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


spec_out_t *spec_out_cp(spec_out_t *out)
{
    /*

        Copy out

    */

    spec_out_t *outCp = malloc(sizeof(spec_out_t));

    outCp -> size = out -> size;

    outCp -> labels = malloc(sizeof(char*) * out -> size);

    outCp -> err = malloc(sizeof(bool) * out -> size);
    outCp -> log = malloc(sizeof(bool) * out -> size);

    for (size_t i = 0; i < out -> size; i++)
      {
        outCp -> labels[i] = misc_scat(1, out -> labels[i]);

        outCp -> err[i] = out -> err[i];
        outCp -> log[i] = out -> log[i];
      }

    outCp -> file = misc_scat(1, out -> file);

    outCp -> binary = out -> binary;
    outCp -> precision = out -> precision;

    return outCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_out_add_label(spec_out_t *out, const char *label)
{
    /*

        Add an output label to out

    */

    /* Label already contained in out -> labels */
    if (misc_ssearch(out -> labels, out -> size, label, NULL))
        return 0;

    /* Increment the output size and reallocate memory */
    out -> size += 1;
    out -> labels = realloc(out -> labels, sizeof(char*) * out -> size);
    out -> err = realloc(out -> err, sizeof(bool) * out -> size);
    out -> log = realloc(out -> log, sizeof(bool) * out -> size);

    /* Add the label */
    out -> labels[out -> size - 1] = misc_scat(1, label);

    /* Set the flags to _false_ */
    out -> err[out -> size - 1] = _false_;
    out -> log[out -> size - 1] = _false_;

    return 0;
}


int spec_out_rm_label(spec_out_t *out, const char *label)
{
    /*

        Remove an output label from out

    */

    bool success;

    /* Get the index */
    size_t index = misc_ssearch_index(out -> labels, out -> size, label, &success);

    /* Not in the array */
    if (!success)
        return 0;

    /* Decrease the size of the array */
    out -> size -= 1;

    /* Push the string to the end of the array */
    for (size_t i = index; i < out -> size; i++)
      {
        misc_swap(&(out -> labels)[i], &(out -> labels)[i+1], "s");
        misc_swap(&(out -> err)[i], &(out -> err)[i+1], "b");
        misc_swap(&(out -> log)[i], &(out -> log)[i+1], "b");
      }

    /* Free the last label */
    free(out -> labels[out -> size]);

    /* Reallocate memory */
    out -> labels = realloc(out -> labels, sizeof(char*) * out -> size);
    out -> err = realloc(out -> err, sizeof(bool*) * out -> size);
    out -> log = realloc(out -> log, sizeof(bool*) * out -> size);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_out_set_err(spec_out_t *out, const char *label, bool err)
{
    /*

        Set the error flag for the output label in out

    */

    bool success;

    /* Get the index */
    size_t index = misc_ssearch_index(out -> labels, out -> size, label, &success);

    /* Not in the array */
    if (!success)
      {
        printf("Cannot set the integration error flag for non-existent label '%s' in a spec_out_t struct!\n", label);

        return 0;
      }

    out -> err[index] = err;

    return 0;
}


int spec_out_set_log(spec_out_t *out, const char *label, bool log)
{
    /*

        Set the logarithimic derivative flag for the output label in out

    */

    bool success;

    /* Get the index */
    size_t index = misc_ssearch_index(out -> labels, out -> size, label, &success);

    /* Not in the array */
    if (!success)
      {
        printf("Cannot set the logarithmic derivative flag for non-existent label '%s' in a spec_out_t struct!\n", label);

        return 0;
      }

    out -> log[index] = log;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_out_set_file(spec_out_t *out, const char *file)
{
    /*

        Set the file name in out (if file is NULL nothing will be output to file)

    */

    misc_scp(&(out -> file), file);

    return 0;
}


int spec_out_set_binary(spec_out_t *out, bool binary)
{
    /*

        Set the binary flag for out

    */

    out -> binary = binary;

    return 0;
}


int spec_out_set_precision(spec_out_t *out, int precision)
{
    /*

        Set the precision in the output file for out -> also set the binary flag to _false_

    */

    out -> binary = _false_;
    out -> precision = precision;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Spectra Struct   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_dat_t *spec_dat_new(const char *id)
{
    /*

        New spec_dat_t struct

    */

    spec_dat_t *specDat = malloc(sizeof(spec_dat_t));

    /* Info struct */
    _spec_info_t *info = _spec_info_get_struct(id);

    /* Id */
    specDat -> id = info -> id;

    /* Output */
    specDat -> out = NULL;

    /* Print */
    specDat -> flags = NULL;

    return specDat;
}


/*  ------------------------------------------------------------------------------------------------------  */


spec_dat_t *spec_dat_free(spec_dat_t *specDat)
{
    /*

        Free a spec_dat_t struct

    */

    /* Check for NULL */
    if (specDat == NULL)
        return NULL;

    /* Free contents */

    /* Out */
    specDat -> out = spec_out_free(specDat -> out);

    /* Flags */
    specDat -> flags = print_flags_free(specDat -> flags);

    /* Free spec itself */
    free(specDat);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_dat_set_out(spec_dat_t *specDat, spec_out_t *out)
{
    /*

        Set the output parameters for spec

    */

    specDat -> out = spec_out_free(specDat -> out);
    specDat -> out = spec_out_cp(out);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dat_set_flags(spec_dat_t *specDat, print_flags_t *flags)
{
    /*

        Set the print parameters for spec

    */

    specDat -> flags = print_flags_free(specDat -> flags);
    specDat -> flags = print_flags_cp(flags);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_out_t *spec_dat_get_out(spec_dat_t *specDat)
{
    /*

        Get the out struct from specDat

    */

    return specDat -> out;
}


print_flags_t *spec_dat_get_flags(spec_dat_t *specDat)
{
    /*

        Get the print struct from specDat

    */

    return specDat -> flags;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Default Settings   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_dat_t *spec_dat_pnl_default(void)
{
    /*

        Default parameters for the non-linear power spectrum function

    */

    spec_dat_t *specDat = spec_dat_new(_idSpecPnl_);


    /* Output */

    spec_out_t *out = spec_out_new();

    /* Labels and errors */
    for (size_t i = 0; i < _pnlInfo -> partsSize; i++)
      {
        spec_out_add_label(out, _pnlInfo -> partsLabels[i]);
        spec_out_set_err(out, _pnlInfo -> partsLabels[i], *_pnlInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _pnlInfo -> id, _outExt);
    spec_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    spec_out_set_precision(out, 16);
    spec_out_set_binary(out, _true_);

    spec_dat_set_out(specDat, out);

    out = spec_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    spec_dat_set_flags(specDat, flags);

    flags = print_flags_free(flags);


    return specDat;
}


/*  ------------------------------------------------------------------------------------------------------  */


spec_dat_t *spec_dat_dpnl_default(void)
{
    /*

        Default parameters for the non-linear power spectrum derivatives function

    */

    spec_dat_t *specDat = spec_dat_new(_idSpecDPnl_);


    /* Output */

    spec_out_t *out = spec_out_new();

    /* Labels and errors */
    for (size_t i = 0; i < _dpnlInfo -> partsSize; i++)
      {
        spec_out_add_label(out, _dpnlInfo -> partsLabels[i]);
        spec_out_set_err(out, _dpnlInfo -> partsLabels[i], *_dpnlInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _dpnlInfo -> id, _outExt);
    spec_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    spec_out_set_precision(out, 16);
    spec_out_set_binary(out, _true_);

    spec_dat_set_out(specDat, out);

    out = spec_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    spec_dat_set_flags(specDat, flags);

    flags = print_flags_free(flags);


    return specDat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_dat_t *spec_dat_btr_default(void)
{
    /*

        Default parameters for the tree-level bispectrum function

    */

    spec_dat_t *specDat = spec_dat_new(_idSpecBtr_);


    /* Output */

    spec_out_t *out = spec_out_new();

    /* Labels and errors */
    for (size_t i = 0; i < _btrInfo -> partsSize; i++)
      {
        spec_out_add_label(out, _btrInfo -> partsLabels[i]);
        spec_out_set_err(out, _btrInfo -> partsLabels[i], *_btrInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _btrInfo -> id, _outExt);
    spec_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    spec_out_set_precision(out, 16);
    spec_out_set_binary(out, _true_);

    spec_dat_set_out(specDat, out);

    out = spec_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    spec_dat_set_flags(specDat, flags);

    flags = print_flags_free(flags);


    return specDat;
}


/*  ------------------------------------------------------------------------------------------------------  */


spec_dat_t *spec_dat_dbtr_default(void)
{
    /*

        Default parameters for the tree-level bispectrum function

    */

    spec_dat_t *specDat = spec_dat_new(_idSpecDBtr_);


    /* Output */

    spec_out_t *out = spec_out_new();

    /* Labels and errors */
    for (size_t i = 0; i < _dbtrInfo -> partsSize; i++)
      {
        spec_out_add_label(out, _dbtrInfo -> partsLabels[i]);
        spec_out_set_err(out, _dbtrInfo -> partsLabels[i], *_dbtrInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _dbtrInfo -> id, _outExt);
    spec_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    spec_out_set_precision(out, 16);
    spec_out_set_binary(out, _true_);

    spec_dat_set_out(specDat, out);

    out = spec_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    spec_dat_set_flags(specDat, flags);

    flags = print_flags_free(flags);


    return specDat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CALCULATE SPECTRA     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Setup Data Struct   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _setup_dat_append_y(_spec_info_t *info, spec_dat_t *specDat, dat_t *dat, char *label, bool err, size_t index, int *depth);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static dat_t *_spec_dat_setup_dat(_spec_info_t *info, spec_dat_t *specDat)
{
    /*

        Setup the output dat_t struct

    */


    /* Variables of the data */

    sample_shape_t *sampleShape = flss_get_sample_shape(specDat -> id);

    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    /* Data dimensions */
    size_t xDim = 1 + 3 * (sampleShape -> dim - 1) - 1; // 1 for temporal; 3 * (dim - 1) - 1 for spatial (3 * dim for 3 spatial vectors, -3 for delta and -1 for free choice of one variable)
    size_t yDim = 0;
    size_t size = sampleArgZ -> size * sampleShape -> size;

    /* Allocate memory */
    dat_t *dat = dat_new(xDim, yDim, size);

    /* Set the x labels */
    char *label, buf[5];
    size_t position = 0;

    label = misc_scat(3, "'", _idVarZ_, "'");
    dat_set_xlabel(dat, position++, label);
    free(label);

    for (size_t i = 0; i < sampleShape -> dimLength; i++)
      {
        if (sampleShape -> dimLength == 1)
            label = misc_scat(3, "'", _idVarK_, "'");
        else
          {
            snprintf(buf, 5, "%ld", i + 1);
            label = misc_scat(4, "'", _idVarK_, buf, "'");
          }

        dat_set_xlabel(dat, position++, label);
        free(label);
      }

    for (size_t i = 0; i < sampleShape -> dimOrientation; i++)
      {
        if (sampleShape -> dimOrientation == 1)
            label = misc_scat(3, "'", _idVarMu_, "'");
        else
          {
            snprintf(buf, 5, "%ld", i + 1);
            label = misc_scat(4, "'", _idVarMu_, buf, "'");
          }

        dat_set_xlabel(dat, position++, label);
        free(label);
      }


    /* Data columns */

    /* Search for the labels and append new columns to the data struct */
    for (size_t i = 0; i < info -> partsSize; i++)
      {
        if (!(*(info -> partsExist[i])))
            continue;

        bool success;
        size_t index = misc_ssearch_index(specDat -> out -> labels, specDat -> out -> size, info -> partsLabels[i], &success);

        if (info -> partsDerivLog != NULL)
          {
            info -> partsDerivLog[i] = specDat -> out -> log[index];
          }

        if (info -> partsDerivVals != NULL)
          {
            info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);
          }

        if (success)
          {
            int depth = 0;

            _setup_dat_append_y(info, specDat, dat, info -> partsLabels[i], *info -> partsHasErr[i] && specDat -> out -> err[index], i, &depth);
          }
      }

    /* Did not find any matching label */
    if (dat -> yDim == 0)
      {
        printf("Labels for the '%s' function were either not provided or did not fit any of the expected labels.\n", info -> id);
      }

    return dat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _setup_dat_append_y(_spec_info_t *info, spec_dat_t *specDat, dat_t *dat, char *label, bool err, size_t index, int *depth)
{
    /*

        Append y columns to dat

    */

    /* No derivatives */
    if (info -> partsDerivMult == NULL)
        goto appendCol;

    /* Variables */
    sample_shape_t *sampleShape = flss_get_sample_shape(specDat -> id);

    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    sample_raw_t *sampleRawK = sampleShape -> sampleRawLength;
    sample_arg_t *sampleArgK = sampleRawK -> sampleArg;

    sample_raw_t *sampleRawMu = sampleShape -> sampleRawOrientation;
    sample_arg_t *sampleArgMu = sampleRawMu -> sampleArg;

    int currDepth = *depth;
    *depth += 1;

    /* Redshift z */
    if (currDepth == 0)
      {
        /* Full multiplicity */
        if (info -> partsDerivMult[index][0] == 1)
          {
            /* z sample as dat struct */
            dat_t *sampleZDat = dat_new(1, 0, sampleArgZ -> size);

            for (size_t i = 0; i < sampleArgZ -> size; i++)
                dat_set_value(sampleZDat, 0, i, 'x', sampleRawZ -> array[i]);

            /* Derivative values */
            dat_t *derivVals = dat_cat(2, info -> partsDerivVals[index], sampleZDat);
            info -> partsDerivVals[index] = dat_free(info -> partsDerivVals[index]);
            info -> partsDerivVals[index] = dat_cat(1, derivVals);
            derivVals = dat_free(derivVals);
            sampleZDat = dat_free(sampleZDat);

            dat_set_xlabel(info -> partsDerivVals[index], info -> partsDerivVals[index] -> xDim - 1, (char*) _idVarZ_);

            for (size_t i = 0; i < sampleArgZ -> size; i++)
              {
                /* Get the z value as a string */
                char *zToString = misc_dtos(sampleRawZ -> array[i], __XSBUFFER__);

                /* Total label */
                char *totLabel = misc_scat(6, label, " [", _idVarZ_, " = ", zToString, "]");

                /* Append a new column */
                _setup_dat_append_y(info, specDat, dat, totLabel, err, index, depth);
                *depth = currDepth + 1;

                /* Free memory */
                free(zToString);
                free(totLabel);
              }

            return 0;
          }

        /* Partial multiplicity */
        else if (info -> partsDerivMult[index][0] == -1)
          {
            /* Total label */
            char *totLabel = misc_scat(6, label, " [", _idVarZ_, " = ", _idVarZ_, "']");

            /* Append a new column */
            _setup_dat_append_y(info, specDat, dat, totLabel, err, index, depth);

            /* Free memory */
            free(totLabel);

            return 0;
          }

        /* No multiplicity */
        else
          {
            _setup_dat_append_y(info, specDat, dat, label, err, index, depth);

            return 0;
          }
      }


    /* Scale k */
    if (currDepth == 1)
      {
        /* Full multiplicity */
        if (info -> partsDerivMult[index][1] == 1)
          {
            /* k sample as dat struct */
            dat_t *sampleKDat = dat_new(1, 0, sampleArgK -> size);

            for (size_t i = 0; i < sampleArgK -> size; i++)
                dat_set_value(sampleKDat, 0, i, 'x', sampleRawK -> array[i]);

            /* Derivative values */
            dat_t *derivVals = dat_cat(2, info -> partsDerivVals[index], sampleKDat);
            info -> partsDerivVals[index] = dat_free(info -> partsDerivVals[index]);
            info -> partsDerivVals[index] = dat_cat(1, derivVals);
            derivVals = dat_free(derivVals);
            sampleKDat = dat_free(sampleKDat);

            dat_set_xlabel(info -> partsDerivVals[index], info -> partsDerivVals[index] -> xDim - 1, (char*) _idVarK_);

            for (size_t i = 0; i < sampleArgK -> size; i++)
              {
                /* Get the k value as a string */
                char *kToString = misc_dtos(sampleRawK -> array[i], __XSBUFFER__);

                /* Total label */
                char *totLabel = misc_scat(6, label, " [", _idVarK_, " = ", kToString, "]");

                /* Append a new column */
                _setup_dat_append_y(info, specDat, dat, totLabel, err, index, depth);
                *depth = currDepth + 1;

                /* Free memory */
                free(kToString);
                free(totLabel);
              }

            return 0;
          }

        /* Partial multiplicity */
        else if (info -> partsDerivMult[index][1] == -1)
          {
            /* Total label */
            char *totLabel = misc_scat(6, label, " [", _idVarK_, " = ", _idVarK_, "']");

            /* Append a new column */
            _setup_dat_append_y(info, specDat, dat, totLabel, err, index, depth);

            /* Free memory */
            free(totLabel);

            return 0;
          }

        /* No multiplicity */
        else
          {
            _setup_dat_append_y(info, specDat, dat, label, err, index, depth);

            return 0;
          }
      }


    /* Angle mu */
    if (currDepth == 2)
      {
        /* Full multiplicity */
        if (info -> partsDerivMult[index][2] == 1)
          {
            /* mu sample as dat struct */
            dat_t *sampleMuDat = dat_new(1, 0, sampleArgMu -> size);

            for (size_t i = 0; i < sampleArgMu -> size; i++)
                dat_set_value(sampleMuDat, 0, i, 'x', sampleRawMu -> array[i]);

            /* Derivative values */
            dat_t *derivVals = dat_cat(2, info -> partsDerivVals[index], sampleMuDat);
            info -> partsDerivVals[index] = dat_free(info -> partsDerivVals[index]);
            info -> partsDerivVals[index] = dat_cat(1, derivVals);
            derivVals = dat_free(derivVals);
            sampleMuDat = dat_free(sampleMuDat);

            dat_set_xlabel(info -> partsDerivVals[index], info -> partsDerivVals[index] -> xDim - 1, (char*) _idVarMu_);

            for (size_t i = 0; i < sampleArgMu -> size; i++)
              {
                /* Get the mu value as a string */
                char *muToString = misc_dtos(sampleRawMu -> array[i], __XSBUFFER__);

                /* Total label */
                char *totLabel = misc_scat(6, label, " [", _idVarMu_, " = ", muToString, "]");

                /* Append a new column */
                _setup_dat_append_y(info, specDat, dat, totLabel, err, index, depth);
                *depth = currDepth + 1;

                /* Free memory */
                free(muToString);
                free(totLabel);
              }

            return 0;
          }

        /* Partial multiplicity */
        else if (info -> partsDerivMult[index][2] == -1)
          {
            /* Total label */
            char *totLabel = misc_scat(6, label, " [", _idVarMu_, " = ", _idVarMu_, "']");

            /* Append a new column */
            _setup_dat_append_y(info, specDat, dat, totLabel, err, index, depth);

            /* Free memory */
            free(totLabel);

            return 0;
          }

        /* No multiplicity */
        else
          {
            _setup_dat_append_y(info, specDat, dat, label, err, index, depth);

            return 0;
          }
      }

    appendCol:

    /* Append another column */
    dat_append_y(dat, label);

    /* Append another column for errors */
    if (err)
      {
        /* Total label */
        char *totLabel = misc_scat(3, label, " ", _idIntError_);

        /* Append a new column */
        dat_append_y(dat, totLabel);

        /* Free memory */
        free(totLabel);
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------   Integration Bins for Derivatives   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


// TODO: Improve this...


static int _spec_dpnl_setup_bin(void)
{
    /*

        Setup the integration bins for the power spectrum derivatives

    */

    int threadsNum = omp_get_max_threads();

    /* Integral bins */
    _pnl1LoopBin = malloc(sizeof(double*) * (size_t) threadsNum);
    _dpnlK1LoopBin = malloc(sizeof(double*) * (size_t) threadsNum);
    _dpnlMu1LoopBin = malloc(sizeof(double*) * (size_t) threadsNum);

    for (size_t i = 0; i < (size_t) threadsNum; i++)
      {
        _pnl1LoopBin[i] = NULL;
        _dpnlK1LoopBin[i] = NULL;
        _dpnlMu1LoopBin[i] = NULL;
      }

    return 0;
}


static int _spec_dpnl_reset_bin(int thread)
{
    /*

        Reset the integration bin for the power spectrum derivatives for a single thread

    */

    if (_pnl1LoopBin[thread] != NULL) free(_pnl1LoopBin[thread]);
    _pnl1LoopBin[thread] = NULL;

    if (_dpnlK1LoopBin[thread] != NULL) free(_dpnlK1LoopBin[thread]);
    _dpnlK1LoopBin[thread] = NULL;

    if (_dpnlMu1LoopBin[thread] != NULL) free(_dpnlMu1LoopBin[thread]);
    _dpnlMu1LoopBin[thread] = NULL;

    return 0;
}


static int _spec_dpnl_free_bin(void)
{
    /*

        Free the integration bins for the power spectrum derivatives

    */

    int threadsNum = omp_get_max_threads();

    if (_pnl1LoopBin != NULL)
      {
        for (size_t i = 0; i < (size_t) threadsNum; i++)
          {
            if (_pnl1LoopBin[i] != NULL) free(_pnl1LoopBin[i]);
          }

        free(_pnl1LoopBin);

        _pnl1LoopBin = NULL;
      }

    if (_dpnlK1LoopBin != NULL)
      {
        for (size_t i = 0; i < (size_t) threadsNum; i++)
          {
            if (_dpnlK1LoopBin[i] != NULL) free(_dpnlK1LoopBin[i]);
          }

        free(_dpnlK1LoopBin);

        _dpnlK1LoopBin = NULL;
      }

    if (_dpnlMu1LoopBin != NULL)
      {
        for (size_t i = 0; i < (size_t) threadsNum; i++)
          {
            if (_dpnlMu1LoopBin[i] != NULL) free(_dpnlMu1LoopBin[i]);
          }

        free(_dpnlMu1LoopBin);

        _dpnlMu1LoopBin = NULL;
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Spectra Function   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


dat_t *spec_dat_poly(spec_dat_t *specDat)
{
    /*

        Calculate either the power spectrum or bispectrum for several configurations

    */


    /* Exit if nothing should be calculated */

    /* specDat itself is NULL */
    if (specDat == NULL)
      {
        return NULL;
      }

    /* Outsize is 0 */
    if (specDat -> out -> size == 0)
      {
        return NULL;
      }


    /* Get parameters from specDat */

    /* Variables */
    sample_shape_t *sampleShape = flss_get_sample_shape(specDat -> id);

    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    sample_raw_t *sampleRawK = sampleShape -> sampleRawLength;
    sample_arg_t *sampleArgK = sampleRawK -> sampleArg;

    sample_raw_t *sampleRawMu = sampleShape -> sampleRawOrientation;
    sample_arg_t *sampleArgMu = sampleRawMu -> sampleArg;

    size_t sizes[4] = {sampleArgZ -> size, sampleArgK -> size, sampleArgMu -> size, sampleArgZ -> size * sampleShape -> size};

    /* Info on spectrum */
    _spec_info_t *info = _spec_info_get_struct(specDat -> id);

    /* Get the output data struct */
    dat_t *dat = _spec_dat_setup_dat(info, specDat);

    /* Print struct */
    print_t *print = print_new();

    print_set_id(print, specDat -> id);
    print_set_flags(print, specDat -> flags);

    print_set_sizes(print, sizes);

    /* TODO: Improve this */
    if (!strcmp(specDat -> id, _idSpecDPnl_))
      {
        _spec_dpnl_setup_bin();
      }


        /* Integrate */

    /* Parallel threading */

    #pragma omp parallel

      { // Start pragma parallel

        /* Spectra arguments struct */
        spec_arg_t *specArg = spec_arg_new(info -> id);

        spec_arg_set_dz(specArg, sampleArgZ -> step);
        spec_arg_set_dk(specArg, sampleArgK -> step);
        spec_arg_set_dmu(specArg, sampleArgMu -> step);

        /* Kern struct */
        kern_t *kern = spec_arg_get_kern(specArg);

        /* Deriv struct */
        spec_deriv_t *deriv = spec_arg_get_deriv(specArg);

        /* Print stage */
        #pragma omp single
          {
            print_spec(print);
          }


        /* Calculate the spectrum */

        for (size_t n = 0; n < sampleArgZ -> size; n++)

          { // Start temporal for

            /* Redshift + Fiducials */
            kernels_set_z(kern, sampleRawZ -> array[n]);

            /* Print fiducials */
            #pragma omp single
              {
                print_set_kern(print, kern);

                print_spec(print);
              }


            /* Calculate the spectra for the current redshift */

            #pragma omp for schedule(dynamic)

            for (size_t m = 0; m < sampleShape -> size; m++)

              { // Start spatial for

                /* Current row in the data array */
                size_t row = sampleShape -> size * n + m;

                /* Column in the data array */
                size_t col = 0;

                /* Set kernel variables */
                shape_t *shape = sampleShape -> arrayShape[m];

                for (size_t i = 0; i < shape -> dim; i++)
                  {
                    kernels_qset_k(kern, i, shape_get_vertex_length(shape, i));
                    kernels_qset_mu(kern, i, shape_get_vertex_orientation(shape, i));

                    for (size_t j = i + 1; j < shape -> dim; j++)
                        kernels_qset_nu(kern, i, j, shape_get_vertex_angle(shape, i, j));
                  }

                /* Print what every thread is calculating. */

                #pragma omp critical
                  {
                    print_set_thread(print, omp_get_thread_num());
                    print_set_kern(print, kern);

                    print_set_subfinished(print, 0);

                    print_spec(print);
                  }


                /* Store z, k and mu */

                dat_set_value(dat, col++, row, 'x', kern -> z);

                for (size_t i = 0; i < sampleShape -> dimLength; i++)
                    dat_set_value(dat, col++, row, 'x', kernels_qget_k(kern, i));

                for (size_t i = 0; i < sampleShape -> dimOrientation; i++)
                    dat_set_value(dat, col++, row, 'x', kernels_qget_mu(kern, i));


                /* Calculate the spectrum contributions */

                /* Reset column for y data */
                col = 0;

                double result[3];
                double resultFull[2];

                resultFull[0] = 0.;
                resultFull[1] = 0.;

                for (size_t i = 0; i < info -> partsSize; i++)

                  { // Start info -> partsSize for

                    /* Spectrum */
                    if (info -> partsDerivLog == NULL)

                      { // Start NULL specArg -> deriv

                        /* Calculate the parts of the full result */
                        if (!misc_sin(info -> partsLabels[i], _idFull_))

                          { // Start _idFull_ not in info -> partsLabels[i]

                            (spec_poly(specDat -> id, info -> partsLabels[i]))(specArg, result);

                            /* Calculate the full result */
                            resultFull[0] += result[0];
                            resultFull[1] = sqrt(resultFull[1]*resultFull[1] + result[1]*result[1]);

                            /* Add the result to the data struct */
                            if (col < dat -> yDim && misc_sin(dat -> yLabels[col], info -> partsLabels[i]))
                              {
                                dat_set_value(dat, col++, row, 'y', result[0]);

                                /* Add the integral error to the data struct */
                                if (col < dat -> yDim && misc_sin(dat -> yLabels[col], _idIntError_))
                                  {
                                    dat_set_value(dat, col++, row, 'y', (result[0] == 0.) ? result[1] : result[1] / result[0]);
                                  }
                              }

                          } // End _idFull_ not in info -> partsLabels[i]

                        /* Full result */
                        else

                          { // Start _idFull_ in info -> partsLabels[i]

                            /* Add the full result to the data struct */
                            if (col < dat -> yDim && misc_sin(dat -> yLabels[col], info -> partsLabels[i]))
                              {
                                dat_set_value(dat, col++, row, 'y', resultFull[0]);

                                /* Add the full integral error to the data struct */
                                if (col < dat -> yDim && misc_sin(dat -> yLabels[col], _idIntError_))
                                  {
                                    dat_set_value(dat, col++, row, 'y', (resultFull[0] == 0.) ? resultFull[1] : resultFull[1] / resultFull[0]);
                                  }
                              }

                          } // End _idFull_ in info -> partsLabels[i]

                      } // End NULL specArg -> deriv

                    /* Spectrum Derivative */
                    else

                      { // Start not NULL specArg -> deriv

                        /* Logarithmic derivative */
                        spec_deriv_set_log(deriv, info -> partsDerivLog[i]);

                        /* No or partial multiplicity */
                        if (info -> partsDerivVals[i] == NULL)
                          {
                            /* Set z for partial multiplicity (TODO: What about k and mu? Have different options...) */
                            spec_deriv_set_z(deriv, kern -> z);

                            /* Calculate the derivative and add the result to the data struct */
                            if (col < dat -> yDim && misc_sin(dat -> yLabels[col], info -> partsLabels[i]))
                              {
                                /* Calculate the derivative */
                                (spec_poly(info -> id, info -> partsLabels[i]))(specArg, result);

                                dat_set_value(dat, col++, row, 'y', result[0]);

                                /* Add the integral error to the data struct */
                                if (col < dat -> yDim && misc_sin(dat -> yLabels[col], _idIntError_))
                                  {
                                    dat_set_value(dat, col++, row, 'y', (result[0] == 0.) ? result[1] : result[1] / result[0]);
                                  }
                              }

                            continue;
                          }

                        /* Full multiplicity */
                        for (size_t i1 = 0; i1 < info -> partsDerivVals[i] -> size; i1++)
                          {
                            /* Set the deriv variables */
                            for (size_t i2 = 0; i2 < info -> partsDerivVals[i] -> xDim; i2++)
                                spec_deriv_set_var(deriv, info -> partsDerivVals[i] -> xLabels[i2], dat_get_value(info -> partsDerivVals[i], i2, i1, 'x'));

                            /* Calculate the derivative and add the result to the data struct */
                            if (col < dat -> yDim && misc_sin(dat -> yLabels[col], info -> partsLabels[i]))
                              {
                                /* Calculate the derivative */
                                (spec_poly(specDat -> id, info -> partsLabels[i]))(specArg, result);

                                dat_set_value(dat, col++, row, 'y', result[0]);

                                /* Add the integral error to the data struct */
                                if (col < dat -> yDim && misc_sin(dat -> yLabels[col], _idIntError_))
                                  {
                                    dat_set_value(dat, col++, row, 'y', (result[0] == 0.) ? result[1] : fabs(result[1] / result[0]));
                                  }
                              }
                          }

                      } // End not NULL specArg -> deriv

                  } // End info -> partsSize for

                /* Thread has finished */
                #pragma omp critical
                  {
                    print_set_subfinished(print, 1);
                    print_update_progress(print, 1. / ((double) dat -> size));

                    print_spec(print);
                  }

                /* TODO: Improve this */
                if (!strcmp(specDat -> id, _idSpecDPnl_))
                  {
                    _spec_dpnl_reset_bin(omp_get_thread_num());
                  }

              } // End spatial for

          } // End temporal for


        /* Stage has finished */

        #pragma omp single
          {
            print_set_stagefinished(print, 1);
            print_set_progress(print, 1.);

            print_spec(print);
          }


        /* Free memory */

        specArg = spec_arg_free(specArg);

      } // End pragma parallel


    /* Free memory */

    print = print_free(print);


    /* Write the result to file */

    if (specDat -> out -> file != NULL && specDat -> out -> precision != 0)
      {
        char *outFile = misc_scat(3, __FISHERLSSDIR__, _outDir, specDat -> out -> file);

        printf("%s\n", outFile);

        if (specDat -> out -> binary)
            dat_output(outFile, dat, NULL);

        else
            dat_output(outFile, dat, &specDat -> out -> precision);

        free(outFile);
      }

    /* TODO: Improve this... */
    if (!strcmp(specDat -> id, _idSpecDPnl_))
      {
        _spec_dpnl_free_bin();
      }

    return dat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     POWER SPECTRUM     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------   (Direct) Non-Linear Power Spectrum   -------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_pnl_tree(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the tree level power spectrum.

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_get_k(kern, 0);

    /* Smoothing */
    double smooth = kernels_smooth(kern);

    /* Result */
    result[0] = pow(kern -> growth * kernels_z1(kern), 2.) * smooth * _fidPk_(&k, _fidParamsPk_);
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_pnl_p22(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the p22 contribution to the one-loop correction.

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Tree-level reuslt */
    if (_pnlInfo -> loopOrder == 0)
      {
        return 0;
      }

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /* Integrate */

    /* One-loop */
    double result1Loop[3] = {0., 0., 0.};
    intgrt_t *intgrt1Loop = integrate_cp(_pnlInfo -> loopIntgrt[0]);
    integrate_set_params(intgrt1Loop, kern);

    integrate(integrand_spec_pnl_p22_1loop, intgrt1Loop, result1Loop);

    intgrt1Loop = integrate_free(intgrt1Loop);

    /* Get the result */
    double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
    result1Loop[0] *= factor1Loop;
    result1Loop[1] *= factor1Loop;

    result[0] += result1Loop[0];
    result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    return 0;
}


int spec_pnl_p13(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the P13 contribution to the one-loop correction.

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Tree-level result */
    if (_pnlInfo -> loopOrder == 0)
      {
        return 0;
      }

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /* Integrate */

    /* One-loop */
    double result1Loop[3] = {0., 0., 0.};
    intgrt_t *intgrt1Loop = integrate_cp(_pnlInfo -> loopIntgrt[0]);
    integrate_set_params(intgrt1Loop, kern);

    integrate(integrand_spec_pnl_p13_1loop, intgrt1Loop, result1Loop);

    intgrt1Loop = integrate_free(intgrt1Loop);

    /* Get the result */
    double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
    result1Loop[0] *= factor1Loop;
    result1Loop[1] *= factor1Loop;

    result[0] += result1Loop[0];
    result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_pnl_ctr(spec_arg_t *specArg, double *result)
{
    /*

        Counter terms for the power spectrum.

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Tree-level result */
    if (_pnlInfo -> loopOrder == 0)
      {
        return 0;
      }

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Smoothing */
    double smooth = kernels_smooth(kern);

    /* Result */
    result[0] += -2. * (
                        kern -> ctr -> c0
                        + 3./2. * kern -> ctr -> c2 * mu*mu
                        - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                            * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                       )

                        * k*k * pow(kern -> growth, 2.) * _fidPk_(&k, _fidParamsPk_) * smooth;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_pnl_sn(spec_arg_t *specArg, double *result)
{
    /*

        Shot noise for the power spectrum

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Result */
    result[0] = kern -> surv -> sn;
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_pnl_full(spec_arg_t *specArg, double *result)
{
    /*

        Full non-linear power spectrum

    */

    /* Must have specOrder = 2 (must have this in case pnl is calculated for shot noise from bi- or trispectrum, or covpp, ...) */
    kern_t *kern = specArg -> kern;
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Tree-level power spectrum */
    double pnlTree[3];
    spec_pnl_tree(specArg, pnlTree);

    /* One-loop-level power spectrum */
    double pnlLoopP22[3];
    double pnlLoopP13[3];

    spec_pnl_p22(specArg, pnlLoopP22);
    spec_pnl_p13(specArg, pnlLoopP13);

    /* Counter terms */
    double pnlCtr[3];
    spec_pnl_ctr(specArg, pnlCtr);

    /* Shot noise */
    double pnlSn[3];
    spec_pnl_sn(specArg, pnlSn);

    /* Final result */
    result[0] = pnlTree[0] + pnlLoopP22[0] + pnlLoopP13[0] + pnlCtr[0] + pnlSn[0];
    result[1] = sqrt(pnlTree[1]*pnlTree[1] + pnlLoopP22[1]*pnlLoopP22[1] + pnlLoopP13[1]*pnlLoopP13[1] + pnlCtr[1]*pnlCtr[1] + pnlSn[1]*pnlSn[1]);

    /* Reset specOrder */
    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*spec_pnl(const char *label))(spec_arg_t*, double*)
{
    /*

        Return the requested function that calculates part (or full) of the non-linear power spectrum

    */

    /* Tree-Level Power Spectrum */
    if (!strcmp(label, _pnlLabelTree))
      {
        return spec_pnl_tree;
      }

    /* One-Loop-Level Power Spectrum (P22) */
    if (!strcmp(label, _pnlLabelP22))
      {
        return spec_pnl_p22;
      }

    /* One-Loop-Level Power Spectrum (P13) */
    if (!strcmp(label, _pnlLabelP13))
      {
        return spec_pnl_p13;
      }

    /* Power Spectrum Counter Terms */
    if (!strcmp(label, _pnlLabelCtr))
      {
        return spec_pnl_ctr;
      }

    /* Power Spectrum Shot Noise */
    if (!strcmp(label, _pnlLabelSn))
      {
        return spec_pnl_sn;
      }

    /* Full Power Spectrum */
    if (!strcmp(label, _pnlLabelFull))
      {
        return spec_pnl_full;
      }

    /* Label not found */
    return _spec_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------   (Indirect) Non-Linear Power Spectrum   ------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_ptr(void *var, void *params)
{
    /*

        Calculate the tree-level power spectrum (plus shot noise)

    */

    (void) params;

    /* Spec arg struct */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Get the tree level power specturm (+ shot noise) */
    double ptrPart[2];
    double ptr = 0.;

    /* Tree level */
    spec_pnl_tree(specArg, ptrPart);
    ptr += ptrPart[0];

    /* Shot noise */
    spec_pnl_sn(specArg, ptrPart);
    ptr += ptrPart[0];

    return ptr;
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_pnl(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum

    */

    (void) params;

    /* No interpolation function */
    if (_pnlInterpInPnl.interp == NULL)
      {
        double result[3];
        spec_pnl_full(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* kern */
    kern_t *kern = specArg -> kern;

    /* Variables (TODO: What if z only has one value in data?) */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* mu can be set to its absolute value TODO: Improve this... (for the case where mu is negative in interpolation function) */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    return interpolate_interp_eval(_pnlInterpInPnl.interpEval, vars, _pnlInterpInPnl.interp, specArg);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------   Non-Linear Power Spectrum (Analytical) Derivatives   ------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_dpnl_a2ga(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(a^(2)_γ)

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_a2ga_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double a2ga = _fidBTSTxA2Ga_(&deriv -> z, _fidParamsBTSTxA2Ga_);

        if (a2ga != 0.)
          {
            result[0] *= a2ga;
            result[1] *= fabs(a2ga);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_d2ga(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(d^(2)_γ)

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_d2ga_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double d2ga = _fidBTSTxD2Ga_(&deriv -> z, _fidParamsBTSTxD2Ga_);

        if (d2ga != 0.)
          {
            result[0] *= d2ga;
            result[1] *= fabs(d2ga);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_a3gaa(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(a^(3)_γa)

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_a3gaa_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double a3gaa = _fidBTSTxA3GaA_(&deriv -> z, _fidParamsBTSTxA3GaA_);

        if (a3gaa != 0.)
          {
            result[0] *= a3gaa;
            result[1] *= fabs(a3gaa);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_a3gab(spec_arg_t *specArg, double *result)
{
    /*

         Calculate dPnl/d(a^(3)_γb)

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_a3gab_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double a3gab = _fidBTSTxA3GaB_(&deriv -> z, _fidParamsBTSTxA3GaB_);

        if (a3gab != 0.)
          {
            result[0] *= a3gab;
            result[1] *= fabs(a3gab);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_d3gaa(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(d^(3)_γa)

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_d3gaa_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double d3gaa = _fidBTSTxD3GaA_(&deriv -> z, _fidParamsBTSTxD3GaA_);

        if (d3gaa != 0.)
          {
            result[0] *= d3gaa;
            result[1] *= fabs(d3gaa);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_d3gab(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(d^(3)_γb)

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_d3gab_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double d3gab = _fidBTSTxD3GaB_(&deriv -> z, _fidParamsBTSTxD3GaB_);

        if (d3gab != 0.)
          {
            result[0] *= d3gab;
            result[1] *= fabs(d3gab);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_b1(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/db1

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 2. * kernels_z1(kern) * kernels_dz1_b1(kern) * pow(kern -> growth, 2.) * smooth * _fidPk_(&k, _fidParamsPk_);
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 2. * (
                              kern -> ctr -> c4 * pow(kern -> rsd -> f * k*mu, 4.)
                                 * (kern -> bias -> b1 + kern -> rsd -> f * mu*mu)
                             )

                            * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_b1_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double b1 = _fidBiasxB1_(&deriv -> z, _fidParamsBiasxB1_);

        if (b1 != 0.)
          {
            result[0] *= b1;
            result[1] *= fabs(b1);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_b2(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/db2

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_b2_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double b2 = _fidBiasxB2_(&deriv -> z, _fidParamsBiasxB2_);

        if (b2 != 0.)
          {
            result[0] *= b2;
            result[1] *= fabs(b2);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_c2ga(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dc^(2)_γ

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_c2ga_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double c2ga = _fidBiasxC2Ga_(&deriv -> z, _fidParamsBiasxC2Ga_);

        if (c2ga != 0.)
          {
            result[0] *= c2ga;
            result[1] *= fabs(c2ga);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_bgam3(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/db_Γ3

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = 0.;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        integrate(integrand_spec_dpnl_bgam3_1loop, intgrt1Loop, result1Loop);

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double bgam3 = _fidBiasxBGam3_(&deriv -> z, _fidParamsBiasxBGam3_);

        if (bgam3 != 0.)
          {
            result[0] *= bgam3;
            result[1] *= fabs(bgam3);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dpnl_f(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/df

            dPnl(z)/df(z') = (...) + δPnl(z)/δf(z')

        The ellipsis (...) contains the growth factor contribution:

            G(z) = exp(- int f(z') / (1+z') dz')

        the contribution is

                               {   - G(z) / (1+z') dz'   if z0 < z' < z
            dG(z)/d(f(z')) =  {
                               {            0            else

        where z0 is the first redshift bin.

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;
    result[2] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only lower triangular terms */
    if (deriv -> z > z)
      {
        return 0;
      }

    /* Thread */
    int thread = omp_get_thread_num();

    /* Growth contribution (dlogG/df) */
    double dGrowth = - (1. / (1. + deriv -> z)) * specArg -> dz;

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /***   Growth contribution   ***/

      {
        /**  Tree-level contributiom  **/

        double dpTree = 2. * pow(kern -> growth * kernels_z1(kern), 2.) * smooth * _fidPk_(&k, _fidParamsPk_) * dGrowth;
        result[0] += dpTree;


        /**  One-loop contribution  **/

        if (_dpnlInfo -> loopOrder >= 1)
          {
            /* Counter terms contributions */
            double dpCtr = -4. * (
                                   kern -> ctr -> c0
                                   + 3./2. * kern -> ctr -> c2 * mu*mu
                                   - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                       * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                                 )

                                   * k*k * pow(kern -> growth, 2.) * _fidPk_(&k, _fidParamsPk_) * smooth * dGrowth;
            result[0] += dpCtr;

            /* One-loop */
            double result1Loop[3] = {0., 0., 0.};
            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            /* Can use the bin if it was setup beforehand */
            if (_pnl1LoopBin != NULL)
              {
                if (_pnl1LoopBin[thread] == NULL)
                  {
                    _pnl1LoopBin[thread] = malloc(sizeof(double) * 3);

                    integrate(integrand_spec_pnl_1loop, intgrt1Loop, _pnl1LoopBin[thread]);
                  }

                result1Loop[0] = _pnl1LoopBin[thread][0];
                result1Loop[1] = _pnl1LoopBin[thread][1];
              }

            /* Must integrate without bin */
            else
              {
                integrate(integrand_spec_pnl_1loop, intgrt1Loop, result1Loop);
              }

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 8. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth * dGrowth;
            result1Loop[0] *= factor1Loop;
            result1Loop[1] *= factor1Loop;

            result[0] += result1Loop[0];
            result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
            result[2] = result1Loop[2];
          }
      }


    /***   Partial f derivative contribution   ***/

    if (deriv -> z == z)
      {
        /**  Tree-level contribution  **/

        double dpTree = 2. * kernels_z1(kern) * kernels_dz1_f(kern) * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;
        result[0] += dpTree;


        /**  One-loop contribution  **/

        if (_dpnlInfo -> loopOrder >= 1)
          {
            /* Counter terms contributions */
            double dpCtr = 4. * (
                              kern -> ctr -> c4 * pow(kern -> rsd -> f, 3.) * pow(k*mu, 4.)
                                 * ( kern -> bias -> b1 + kern -> rsd -> f * mu*mu )
                                 * ( kern -> rsd -> f * mu*mu
                                        + (kern -> bias -> b1 + kern -> rsd -> f * mu*mu) )
                            )

                             * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;
            result[0] += dpCtr;

            /* One-loop */
            double result1Loop[3] = {0., 0., 0.};
            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            integrate(integrand_spec_dpnl_f_1loop, intgrt1Loop, result1Loop);

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
            result1Loop[0] *= factor1Loop;
            result1Loop[1] *= factor1Loop;

            result[0] += result1Loop[0];
            result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
            result[2] = result1Loop[2];
          }
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double f = _fidRSDxF_(&deriv -> z, _fidParamsRSDxF_);

        if (f != 0.)
          {
            result[0] *= f;
            result[1] *= fabs(f);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_sigv(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dsigv

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Thread */
    int thread = omp_get_thread_num();

    /* Smoothing */
    double dsmooth = kernels_dsmooth_sigv(kern);


    /**  Tree-level contribution  **/

    double dpTree = pow(kernels_z1(kern), 2.) * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * dsmooth;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = -2. * (
                               kern -> ctr -> c0
                               + 3./2. * kern -> ctr -> c2 * mu*mu
                               - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                   * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                             )

                            * k*k * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * dsmooth;
        result[0] += dpCtr;

        /* One-loop */
        double result1Loop[3] = {0., 0., 0.};
        intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
        integrate_set_params(intgrt1Loop, kern);

        /* Can use the bin if it was setup beforehand */
        if (_pnl1LoopBin != NULL)
          {
            if (_pnl1LoopBin[thread] == NULL)
              {
                _pnl1LoopBin[thread] = malloc(sizeof(double) * 3);

                integrate(integrand_spec_pnl_1loop, intgrt1Loop, _pnl1LoopBin[thread]);
              }

            result1Loop[0] = _pnl1LoopBin[thread][0];
            result1Loop[1] = _pnl1LoopBin[thread][1];
          }

        /* Must integrate without bin */
        else
          {
            integrate(integrand_spec_pnl_1loop, intgrt1Loop, result1Loop);
          }

        intgrt1Loop = integrate_free(intgrt1Loop);

        /* Get the result */
        double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * dsmooth;
        result1Loop[0] *= factor1Loop;
        result1Loop[1] *= factor1Loop;

        result[0] += result1Loop[0];
        result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double sigv = _fidRSDxSigv_(&deriv -> z, _fidParamsRSDxSigv_);

        if (sigv != 0.)
          {
            result[0] *= sigv;
            result[1] *= fabs(sigv);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dpnl_d(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(d)

            dPnl/d(d) = -2 Pnl + δPnl/δk δk/δd + δPnl/δmu δmu/δd

        where d = Da / Da' is the angular diameter distance ratio which enters via the AP effect (assuming h = H/H' = 1, d = Da/Da' = 1 -> alpha = 1)
        such that:

            δk/δd = k (mu^2 - 1)
            δmu/δd = mu (1 - mu^2)

        Note: The AP effect does not effect the counter terms.

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Thread */
    int thread = omp_get_thread_num();

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /***   Volume correction contribution   ***/

      {
        /**  Tree-level contributiom  **/

        double dpTree = -2. * pow(kernels_z1(kern), 2.) * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth; // -2 * PTree
        result[0] += dpTree;


        /**  One-loop contribution  **/

        if (_dpnlInfo -> loopOrder >= 1)
          {
            /* Counter terms contributions */
            double dpCtr = 4. * (
                                  kern -> ctr -> c0
                                  + 3./2. * kern -> ctr -> c2 * mu*mu
                                  - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                      * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                                )

                               * k*k * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth; // -2 * Pctr
            result[0] += dpCtr;

            /* One-loop */
            double result1Loop[3] = {0., 0., 0.};
            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            /* Can use the bin if it was setup beforehand */
            if (_pnl1LoopBin != NULL)
              {
                if (_pnl1LoopBin[thread] == NULL)
                  {
                    _pnl1LoopBin[thread] = malloc(sizeof(double) * 3);

                    integrate(integrand_spec_pnl_1loop, intgrt1Loop, _pnl1LoopBin[thread]);
                  }

                result1Loop[0] = _pnl1LoopBin[thread][0];
                result1Loop[1] = _pnl1LoopBin[thread][1];
              }

            /* Must integrate without bin */
            else
              {
                integrate(integrand_spec_pnl_1loop, intgrt1Loop, result1Loop);
              }

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = -4. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;

            result1Loop[0] *= factor1Loop;
            result1Loop[1] *= factor1Loop;

            result[0] += result1Loop[0];
            result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
          }
      }


    /***   Partial D derivative contribution   ***/

    if (fabs(mu) < 1.) // Must have |mu| < 1
      {
        /**  Tree-level contributiom  **/

        double dpTreeK = pow(kernels_z1(kern), 2.) * pow(kern -> growth, 2.) * smooth * _fidDPk_(&k, _fidParamsDPk_);
        double dpTreeMu = 2. * kernels_z1(kern) * kernels_dz1_mu(kern) * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;

        result[0] += (dpTreeMu * mu - dpTreeK * k) * (1. - mu*mu);


        /**  One-loop contribution  **/

        if (_dpnlInfo -> loopOrder >= 1)
          {
            /* Counter terms contributions */
            double dpCtrK = -2. * (
                                  kern -> ctr -> c0
                                  + 3./2. * kern -> ctr -> c2 * mu*mu
                                  - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                      * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                                 )

                               * k*k * _fidDPk_(&k, _fidParamsDPk_) * pow(kern -> growth, 2.) * smooth;

            result[0] += -dpCtrK * k * (1. - mu*mu);

            /* One-loop */
            double result1LoopK[3];
            double result1LoopMu[3];

            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            /* Can only use the bin if it was setup beforehand */
            if (_dpnlK1LoopBin != NULL)
              {
                if (_dpnlK1LoopBin[thread] == NULL)
                  {
                    _dpnlK1LoopBin[thread] = malloc(sizeof(double) * 3);

                    integrate(integrand_spec_dpnl_k_1loop, intgrt1Loop, _dpnlK1LoopBin[thread]);
                  }

                result1LoopK[0] = _dpnlK1LoopBin[thread][0];
                result1LoopK[1] = _dpnlK1LoopBin[thread][1];
              }

            /* Must integrate without bin */
            else
              {
                integrate(integrand_spec_dpnl_k_1loop, intgrt1Loop, result1LoopK);
              }

            /* Can only use the bin if it was setup beforehand */
            if (_dpnlMu1LoopBin != NULL)
              {
                if (_dpnlMu1LoopBin[thread] == NULL)
                  {
                    _dpnlMu1LoopBin[thread] = malloc(sizeof(double) * 3);

                    integrate(integrand_spec_dpnl_mu_1loop, intgrt1Loop, _dpnlMu1LoopBin[thread]);
                  }

                result1LoopMu[0] = _dpnlMu1LoopBin[thread][0];
                result1LoopMu[1] = _dpnlMu1LoopBin[thread][1];
              }

            /* Must integrate without bin */
            else
              {
                integrate(integrand_spec_dpnl_mu_1loop, intgrt1Loop, result1LoopMu);
              }


            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;

            result1LoopK[0] *= factor1Loop;
            result1LoopK[1] *= factor1Loop;

            result1LoopMu[0] *= factor1Loop;
            result1LoopMu[1] *= factor1Loop;

            result[0] += (result1LoopMu[0] * mu - result1LoopK[0] * k) * (1. - mu*mu);
            result[1] = sqrt(result[1]*result[1] + ( pow(result1LoopMu[1] * mu, 2.) + pow(result1LoopK[1] * k, 2.) ) * pow(1. - mu*mu, 2.));
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_h(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/d(h)

            dPnl/d(h) = Pnl + δPnl/δk δk/δh + δPnl/δmu δmu/δh

        where h = H / H' is the Hubble function ratio which enters via the AP effect (assuming h = H/H' = 1, d = Da/Da' = 1 -> alpha = 1).
        such that:

            δk/δh = k mu^2
            δmu/δh = mu (1 - mu^2)

        Note: The AP effect does not effect the counter terms.

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Thread */
    int thread = omp_get_thread_num();

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /***   Volume correction contribution   ***/

      {
        /**  Tree-level contributiom  **/

        double dpTree = pow(kernels_z1(kern), 2.) * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth; // +1 * PTree
        result[0] += dpTree;


        /**  One-loop contribution  **/

        if (_dpnlInfo -> loopOrder >= 1)
          {
            /* Counter terms contributions */
            double dpCtr = -2. * (
                                   kern -> ctr -> c0
                                   + 3./2. * kern -> ctr -> c2 * mu*mu
                                   - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                       * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                                 )

                                * k*k * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth; // +1 * Pctr
            result[0] += dpCtr;

            /* One-loop */
            double result1Loop[3] = {0., 0., 0.};
            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            /* Can use the bin if it was setup beforehand */
            if (_pnl1LoopBin != NULL)
              {
                if (_pnl1LoopBin[thread] == NULL)
                  {
                    _pnl1LoopBin[thread] = malloc(sizeof(double) * 3);

                    integrate(integrand_spec_pnl_1loop, intgrt1Loop, _pnl1LoopBin[thread]);
                  }

                result1Loop[0] = _pnl1LoopBin[thread][0];
                result1Loop[1] = _pnl1LoopBin[thread][1];
              }

            /* Must integrate without bin */
            else
              {
                integrate(integrand_spec_pnl_1loop, intgrt1Loop, result1Loop);
              }

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;

            result1Loop[0] *= factor1Loop;
            result1Loop[1] *= factor1Loop;

            result[0] += result1Loop[0];
            result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
          }
      }


    /***   Partial H derivative contribution   ***/

      {
        /**  Tree-level contributiom  **/

        double dpTreeK = pow(kernels_z1(kern), 2.) * pow(kern -> growth, 2.) * smooth * _fidDPk_(&k, _fidParamsDPk_);
        double dpTreeMu = 2. * kernels_z1(kern) * kernels_dz1_mu(kern) * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;

        result[0] += (dpTreeMu * (1. - mu*mu) + dpTreeK * k*mu) * mu;


        /**  One-loop contribution  **/

        if (_dpnlInfo -> loopOrder >= 1)
          {
            /* Counter terms contributions */
            double dpCtrK = -2. * (
                                  kern -> ctr -> c0
                                  + 3./2. * kern -> ctr -> c2 * mu*mu
                                  - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                      * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                                 )

                               * k*k * _fidDPk_(&k, _fidParamsDPk_) * pow(kern -> growth, 2.) * smooth;

            result[0] += dpCtrK * k * mu*mu;

            /* One-loop */
            double result1LoopK[3] = {0., 0., 0};
            double result1LoopMu[3] = {0., 0., 0.};

            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            /* Can only use the bin if it was setup beforehand */
            if (_dpnlK1LoopBin != NULL)
              {
                if (_dpnlK1LoopBin[thread] == NULL)
                  {
                    _dpnlK1LoopBin[thread] = malloc(sizeof(double) * 3);

                    integrate(integrand_spec_dpnl_k_1loop, intgrt1Loop, _dpnlK1LoopBin[thread]);
                  }

                result1LoopK[0] = _dpnlK1LoopBin[thread][0];
                result1LoopK[1] = _dpnlK1LoopBin[thread][1];
              }

            /* Must integrate without bin */
            else
              {
                integrate(integrand_spec_dpnl_k_1loop, intgrt1Loop, result1LoopK);
              }

            if (fabs(mu) < 1.) // Must have |mu| < 1
              {
                /* Can only use the bin if it was setup beforehand */
                if (_dpnlMu1LoopBin != NULL)
                  {
                    if (_dpnlMu1LoopBin[thread] == NULL)
                      {
                        _dpnlMu1LoopBin[thread] = malloc(sizeof(double) * 3);

                        integrate(integrand_spec_dpnl_mu_1loop, intgrt1Loop, _dpnlMu1LoopBin[thread]);
                      }

                    result1LoopMu[0] = _dpnlMu1LoopBin[thread][0];
                    result1LoopMu[1] = _dpnlMu1LoopBin[thread][1];
                  }

                /* Must integrate without bin */
                else
                  {
                    integrate(integrand_spec_dpnl_mu_1loop, intgrt1Loop, result1LoopMu);
                  }
              }

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;

            result1LoopK[0] *= factor1Loop;
            result1LoopK[1] *= factor1Loop;

            result1LoopMu[0] *= factor1Loop;
            result1LoopMu[1] *= factor1Loop;

            result[0] += (result1LoopMu[0] * (1. - mu*mu) + result1LoopK[0] * k*mu) * mu;
            result[1] = sqrt(result[1]*result[1] + ( pow(result1LoopMu[1] * (1. - mu*mu), 2.) + pow(result1LoopK[1] * k*mu, 2.) ) * mu*mu);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dpnl_c0(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dc0

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = -2. * k*k * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;
        result[0] += dpCtr;
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double c0 = _fidCtrxC0_(&deriv -> z, _fidParamsCtrxC0_);

        if (c0 != 0.)
          {
            result[0] *= c0;
            result[1] *= fabs(c0);
          }
      }

    /* Reset variables */
    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_c2(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dc2

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = -3. * mu*mu * k*k * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;
        result[0] += dpCtr;
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double c2 = _fidCtrxC2_(&deriv -> z, _fidParamsCtrxC2_);

        if (c2 != 0.)
          {
            result[0] *= c2;
            result[1] *= fabs(c2);
          }
      }

    /* Reset variables */
    kern -> specOrder = specOrder;

    return 0;
}


int spec_dpnl_c4(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dc4

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        return 0;
      }

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = 0.;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = pow(kern -> rsd -> f * mu*k, 4.)
                         * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2.)
                         * _fidPk_(&k, _fidParamsPk_) * pow(kern -> growth, 2.) * smooth;
        result[0] += dpCtr;
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double c4 = _fidCtrxC4_(&deriv -> z, _fidParamsCtrxC4_);

        if (c4 != 0.)
          {
            result[0] *= c4;
            result[1] *= fabs(c4);
          }
      }

    /* Reset variables */
    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dpnl_psn(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dPsn

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal term */
    if (fabs(deriv -> z - z) > __ABSTOL__)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Result */
    result[0] = kern -> surv -> sn;
    result[1] = 0.;

    return 0;
}


int spec_dpnl_bsn1(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dBsn1

    */

    (void) specArg;

    /* Result */
    result[0] = 0.;
    result[1] = 0.;

    return 0;
}


int spec_dpnl_bsn2(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl/dBsn2

    */

    (void) specArg;

    /* Result */
    result[0] = 0.;
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dpnl_pk(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dPnl(k_)/dP^(0)(k')

    */

    /* Initialise result */
    result[0] = 0.;
    result[1] = 0.;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Must have specOrder = 2 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 2;

    /* Smoothing */
    double smooth = kernels_smooth(kern);


    /**  Tree-level contribution  **/

    double dpTree = (deriv -> k != k) ? 0.
                                      : pow(kernels_z1(kern) * kern -> growth, 2.) * smooth;
    result[0] += dpTree;


    /**  One-loop contribution  **/

    if (_dpnlInfo -> loopOrder >= 1)
      {
        /* Counter terms contributions */
        double dpCtr = (deriv -> k != k) ? 0.
                                         : -2. * (
                                                    kern -> ctr -> c0
                                                  + 3./2. * kern -> ctr -> c2 * mu*mu
                                                  - 0.5 * kern -> ctr -> c4 * pow(kern -> rsd -> f * mu, 4.)
                                                      * k*k * pow(kern -> bias -> b1 + kern -> rsd -> f * mu*mu, 2)
                                                 )

                                                   * k*k * pow(kern -> growth, 2.) * smooth;

        result[0] += dpCtr;

        /* One-loop */

        /* Volume part */
        if (deriv -> k == k)
          {
            double result1Loop[3] = {0., 0., 0.};
            intgrt_t *intgrt1Loop = integrate_cp(_dpnlInfo -> loopIntgrt[0]);
            integrate_set_params(intgrt1Loop, kern);

            integrate(integrand_spec_dpnl_pk_vol_1loop, intgrt1Loop, result1Loop);

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 2. / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
            result1Loop[0] *= factor1Loop;
            result1Loop[1] *= factor1Loop;

            result[0] += result1Loop[0];
            result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
          }

        /* Angular part */
          {
            kernels_qset_k(kern, 1, deriv -> k); // Must know k' in integrand function

            double result1Loop[3] = {0., 0., 0.};
            intgrt_t *intgrt1Loop = integrate_new();

            size_t dim = 2;
            double upperBounds[2] = {1., 2. * M_PI};
            double lowerBounds[2] = {-1., 0.};

            integrate_set_dim(intgrt1Loop, dim);
            integrate_set_bounds_upper(intgrt1Loop, upperBounds);
            integrate_set_bounds_lower(intgrt1Loop, lowerBounds);

            integrate_set_vegas(intgrt1Loop, _dpnlInfo -> loopIntgrt[0] -> vegas);

            integrate_set_params(intgrt1Loop, kern);

            integrate(integrand_spec_dpnl_pk_ang_1loop, intgrt1Loop, result1Loop);

            intgrt1Loop = integrate_free(intgrt1Loop);

            /* Get the result */
            double factor1Loop = 2. * pow(deriv -> k, 2.) * specArg -> dk / pow(2. * M_PI, 3.) * pow(kern -> growth, 4.) * smooth;
            result1Loop[0] *= factor1Loop;
            result1Loop[1] *= factor1Loop;

            result[0] += result1Loop[0];
            result[1] = sqrt(result[1]*result[1] + result1Loop[1]*result1Loop[1]);
          }
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double pk = _fidPk_(&deriv -> k, _fidParamsPk_);

        if (pk != 0.)
          {
            result[0] *= pk;
            result[1] *= fabs(pk);
          }
      }

    /* Reset variables */
    kernels_qset_k(kern, 1, k);
    kernels_qset_mu(kern, 1, -mu);
    kernels_qset_nu(kern, 0, 1, -1.);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*spec_dpnl(const char *label))(spec_arg_t*, double*)
{
    /*

        Get the non-linear power spectrum derivative function w.r.t. to parameter with given label.

        TODO: Might dereference NULL pointer (_dpnlLabel...) if ini functions are not called...

    */


    /**  Bootstrap parameters  **/

    /* a2ga */
    if (!strcmp(label, _idParamsA2Ga_) || !strcmp(label, _dpnlLabelA2Ga))
      {
        return spec_dpnl_a2ga;
      }

    /* d2ga */
    if (!strcmp(label, _idParamsD2Ga_) || !strcmp(label, _dpnlLabelD2Ga))
      {
        return spec_dpnl_d2ga;
      }

    /* a3gaa */
    if (!strcmp(label, _idParamsA3GaA_) || !strcmp(label, _dpnlLabelA3GaA))
      {
        return spec_dpnl_a3gaa;
      }

    /* a3gab */
    if (!strcmp(label, _idParamsA3GaB_) || !strcmp(label, _dpnlLabelA3GaB))
      {
        return spec_dpnl_a3gab;
      }

    /* d3gaa */
    if (!strcmp(label, _idParamsD3GaA_) || !strcmp(label, _dpnlLabelD3GaA))
      {
        return spec_dpnl_d3gaa;
      }

    /* d3gab */
    if (!strcmp(label, _idParamsD3GaB_) || !strcmp(label, _dpnlLabelD3GaB))
      {
        return spec_dpnl_d3gab;
      }


    /**  Bias parameters  **/

    /* b1 */
    if (!strcmp(label, _idParamsB1_) || !strcmp(label, _dpnlLabelB1))
      {
        return spec_dpnl_b1;
      }

    /* b2 */
    if (!strcmp(label, _idParamsB2_) || !strcmp(label, _dpnlLabelB2))
      {
        return spec_dpnl_b2;
      }

    /* c2ga */
    if (!strcmp(label, _idParamsC2Ga_) || !strcmp(label, _dpnlLabelC2Ga))
      {
        return spec_dpnl_c2ga;
      }

    /* bgam3 */
    if (!strcmp(label, _idParamsBGam3_) || !strcmp(label, _dpnlLabelBGam3))
      {
        return spec_dpnl_bgam3;
      }


    /**  Redshift space distortion parameters  **/

    /* f */
    if (!strcmp(label, _idParamsF_) || !strcmp(label, _dpnlLabelF))
      {
        return spec_dpnl_f;
      }

    /* sigv */
    if (!strcmp(label, _idParamsSigv_) || !strcmp(label, _dpnlLabelSigv))
      {
        return spec_dpnl_sigv;
      }


    /**  AP parameters  **/

    /* d */
    if (!strcmp(label, _idParamsD_) || !strcmp(label, _dpnlLabelD))
      {
        return spec_dpnl_d;
      }

    /* h */
    if (!strcmp(label, _idParamsH_) || !strcmp(label, _dpnlLabelH))
      {
        return spec_dpnl_h;
      }


    /**  Counter term parameters  **/

    /* c0 */
    if (!strcmp(label, _idParamsC0_) || !strcmp(label, _dpnlLabelC0))
      {
        return spec_dpnl_c0;
      }

    /* c2 */
    if (!strcmp(label, _idParamsC2_) || !strcmp(label, _dpnlLabelC2))
      {
        return spec_dpnl_c2;
      }

    /* c4 */
    if (!strcmp(label, _idParamsC4_) || !strcmp(label, _dpnlLabelC4))
      {
        return spec_dpnl_c4;
      }


    /**  Shot noise parameters  **/

    /* Psn */
    if (!strcmp(label, _idParamsPsn_) || !strcmp(label, _dpnlLabelPsn))
      {
        return spec_dpnl_psn;
      }

    /* Bsn1 */
    if (!strcmp(label, _idParamsBsn1_) || !strcmp(label, _dpnlLabelBsn1))
      {
        return spec_dpnl_bsn1;
      }

    /* Bsn2 */
    if (!strcmp(label, _idParamsBsn2_) || !strcmp(label, _dpnlLabelBsn2))
      {
        return spec_dpnl_bsn2;
      }


    /**  Linear power spectrum  **/

    /* Pk */
    if (!strcmp(label, _idParamsPk_) || !strcmp(label, _dpnlLabelPk))
      {
        return spec_dpnl_pk;
      }


    /**  Label not found  **/

    return _spec_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------   (Indirect) Non-Linear Power Spectrum (Analytical) Derivatives   -------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_a2ga(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t a^(2)_γ

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInA2Ga.interp == NULL)
      {
        double result[3];
        spec_dpnl_a2ga(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInA2Ga.interpEval, vars, _dpnlInterpInA2Ga.interp, kern);

    if (deriv -> log)
      {
        double a2ga = _fidBTSTxA2Ga_(&deriv -> z, _fidParamsBTSTxA2Ga_);

        if (a2ga != 0.)
            result *= a2ga;
      }

    return result;
}


static double _spec_dpnl_d2ga(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. d^(2)_γ

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInD2Ga.interp == NULL)
      {
        double result[3];
        spec_dpnl_d2ga(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInD2Ga.interpEval, vars, _dpnlInterpInD2Ga.interp, kern);

    if (deriv -> log)
      {
        double d2ga = _fidBTSTxD2Ga_(&deriv -> z, _fidParamsBTSTxD2Ga_);

        if (d2ga != 0.)
            result *= d2ga;
      }

    return result;
}


static double _spec_dpnl_a3gaa(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. a^(3)_γa

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInA3GaA.interp == NULL)
      {
        double result[3];
        spec_dpnl_a3gaa(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInA3GaA.interpEval, vars, _dpnlInterpInA3GaA.interp, kern);

    if (deriv -> log)
      {
        double a3gaa = _fidBTSTxA3GaA_(&deriv -> z, _fidParamsBTSTxA3GaA_);

        if (a3gaa != 0.)
            result *= a3gaa;
      }

    return result;
}


static double _spec_dpnl_a3gab(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. a^(3)_γb

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInA3GaB.interp == NULL)
      {
        double result[3];
        spec_dpnl_a3gab(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInA3GaB.interpEval, vars, _dpnlInterpInA3GaB.interp, kern);

    if (deriv -> log)
      {
        double a3gab = _fidBTSTxA3GaB_(&deriv -> z, _fidParamsBTSTxA3GaB_);

        if (a3gab != 0.)
            result *= a3gab;
      }

    return result;
}


static double _spec_dpnl_d3gaa(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. d^(3)_γa

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInD3GaA.interp == NULL)
      {
        double result[3];
        spec_dpnl_d3gaa(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInD3GaA.interpEval, vars, _dpnlInterpInD3GaA.interp, kern);

    if (deriv -> log)
      {
        double d3gaa = _fidBTSTxD3GaA_(&deriv -> z, _fidParamsBTSTxD3GaA_);

        if (d3gaa != 0.)
            result *= d3gaa;
      }

    return result;
}


static double _spec_dpnl_d3gab(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. d^(3)_γb

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInD3GaB.interp == NULL)
      {
        double result[3];
        spec_dpnl_d3gab(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInD3GaB.interpEval, vars, _dpnlInterpInD3GaB.interp, kern);

    if (deriv -> log)
      {
        double d3gab = _fidBTSTxD3GaB_(&deriv -> z, _fidParamsBTSTxD3GaB_);

        if (d3gab != 0.)
            result *= d3gab;
      }

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_b1(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. b1

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInB1.interp == NULL)
      {
        double result[3];
        spec_dpnl_b1(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInB1.interpEval, vars, _dpnlInterpInB1.interp, kern);

    if (deriv -> log)
      {
        double b1 = _fidBiasxB1_(&deriv -> z, _fidParamsBiasxB1_);

        if (b1 != 0.)
            result *= b1;
      }

    return result;
}


static double _spec_dpnl_b2(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. b2

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInB2.interp == NULL)
      {
        double result[3];
        spec_dpnl_b2(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInB2.interpEval, vars, _dpnlInterpInB2.interp, kern);

    if (deriv -> log)
      {
        double b2 = _fidBiasxB2_(&deriv -> z, _fidParamsBiasxB2_);

        if (b2 != 0.)
            result *= b2;
      }

    return result;
}


static double _spec_dpnl_c2ga(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. c^(2)_γ

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInC2Ga.interp == NULL)
      {
        double result[3];
        spec_dpnl_c2ga(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInC2Ga.interpEval, vars, _dpnlInterpInC2Ga.interp, kern);

    if (deriv -> log)
      {
        double c2ga = _fidBiasxC2Ga_(&deriv -> z, _fidParamsBiasxC2Ga_);

        if (c2ga != 0.)
            result *= c2ga;
      }

    return result;
}


static double _spec_dpnl_bgam3(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. b_Γ_3

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInBGam3.interp == NULL)
      {
        double result[3];
        spec_dpnl_bgam3(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInBGam3.interpEval, vars, _dpnlInterpInBGam3.interp, kern);

    if (deriv -> log)
      {
        double bgam3 = _fidBiasxBGam3_(&deriv -> z, _fidParamsBiasxBGam3_);

        if (bgam3 != 0.)
            result *= bgam3;
      }

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_f(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivatives

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInF.interp == NULL)
      {
        double result[3];
        spec_dpnl_f(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[4] = {z, k, mu, deriv -> z};

    /* Only "lower diagonal" terms survive */
    if (deriv -> z > z)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInF.interpEval, vars, _dpnlInterpInF.interp, kern);

    if (deriv -> log)
      {
        double f = _fidRSDxF_(&deriv -> z, _fidParamsRSDxF_);

        if (f != 0.)
            result *= f;
      }

    return result;
}


static double _spec_dpnl_sigv(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. sigv

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInSigv.interp == NULL)
      {
        double result[3];
        spec_dpnl_sigv(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInSigv.interpEval, vars, _dpnlInterpInSigv.interp, kern);

    if (deriv -> log)
      {
        double sigv = _fidRSDxSigv_(&deriv -> z, _fidParamsRSDxSigv_);

        if (sigv != 0.)
            result *= sigv;
      }

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_d(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. d

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInD.interp == NULL)
      {
        double result[3];
        spec_dpnl_d(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    return interpolate_interp_eval(_dpnlInterpInD.interpEval, vars, _dpnlInterpInD.interp, kern);
}


static double _spec_dpnl_h(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. h

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInH.interp == NULL)
      {
        double result[3];
        spec_dpnl_h(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    return interpolate_interp_eval(_dpnlInterpInH.interpEval, vars, _dpnlInterpInH.interp, kern);
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_c0(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. c0

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInC0.interp == NULL)
      {
        double result[3];
        spec_dpnl_c0(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInC0.interpEval, vars, _dpnlInterpInC0.interp, kern);

    if (deriv -> log)
      {
        double c0 = _fidCtrxC0_(&deriv -> z, _fidParamsCtrxC0_);

        if (c0 != 0.)
            result *= c0;
      }

    return result;
}


static double _spec_dpnl_c2(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. c2

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInC2.interp == NULL)
      {
        double result[3];
        spec_dpnl_c2(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInC2.interpEval, vars, _dpnlInterpInC2.interp, kern);

    if (deriv -> log)
      {
        double c2 = _fidCtrxC2_(&deriv -> z, _fidParamsCtrxC2_);

        if (c2 != 0.)
            result *= c2;
      }

    return result;
}


static double _spec_dpnl_c4(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. c4

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInC4.interp == NULL)
      {
        double result[3];
        spec_dpnl_c4(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInC4.interpEval, vars, _dpnlInterpInC4.interp, kern);

    if (deriv -> log)
      {
        double c4 = _fidCtrxC4_(&deriv -> z, _fidParamsCtrxC4_);

        if (c4 != 0.)
            result *= c4;
      }

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_psn(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. psn

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInPsn.interp == NULL)
      {
        double result[3];
        spec_dpnl_psn(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    return interpolate_interp_eval(_dpnlInterpInPsn.interpEval, vars, _dpnlInterpInPsn.interp, kern);
}


static double _spec_dpnl_bsn1(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. bsn1

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInBsn1.interp == NULL)
      {
        double result[3];
        spec_dpnl_bsn1(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    return interpolate_interp_eval(_dpnlInterpInBsn1.interpEval, vars, _dpnlInterpInBsn1.interp, kern);
}


static double _spec_dpnl_bsn2(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. bsn2

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInBsn1.interp == NULL)
      {
        double result[3];
        spec_dpnl_bsn2(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[3] = {z, k, mu};

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms survive */
    if (fabs(deriv -> z - z) > __ABSTOL__)
        return 0.;

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    return interpolate_interp_eval(_dpnlInterpInBsn2.interpEval, vars, _dpnlInterpInBsn2.interp, kern);
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dpnl_pk(void *var, void *params)
{
    /*

        Evaluate the interpolation function for the non-linear power spectrum derivative w.r.t. Pk

    */

    (void) params;

    /* No interpolation function */
    if (_dpnlInterpInPk.interp == NULL)
      {
        double result[3];
        spec_dpnl_pk(var, result);

        return result[0];
      }

    /* specArg */
    spec_arg_t *specArg = (spec_arg_t*) var;

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Derivative struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Variables */
    double z = kernels_get_z(kern);
    double k = kernels_get_k(kern, 0);
    double mu = kernels_get_mu(kern, 0);

    double vars[4] = {z, k, mu, deriv -> k};

    /* mu can be set to its positive value */
    vars[2] = fabs(vars[2]);

    /* mu must be less than one */
    if (vars[2] > 1.)
        vars[2] = 1.;

    /* Result */
    double result = interpolate_interp_eval(_dpnlInterpInPk.interpEval, vars, _dpnlInterpInPk.interp, kern);

    if (deriv -> log)
      {
        double pk = _fidPk_(&deriv -> k, _fidParamsPk_);

        if (pk != 0.)
            result *= pk;
      }

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double (*_spec_dpnl(const char *label))(void*, void*)
{
    /*

        Get the non-linear power spectrum derivative function for the given label

    */


    /**  Bootstrap parameters  **/

    /* a2ga */
    if (!strcmp(label, _idParamsA2Ga_) || !strcmp(label, _dpnlLabelA2Ga))
      {
        return _specDPnlA2Ga_;
      }

    /* d2ga */
    if (!strcmp(label, _idParamsD2Ga_) || !strcmp(label, _dpnlLabelD2Ga))
      {
        return _specDPnlD2Ga_;
      }

    /* a3gaa */
    if (!strcmp(label, _idParamsA3GaA_) || !strcmp(label, _dpnlLabelA3GaA))
      {
        return _specDPnlA3GaA_;
      }

    /* a3gab */
    if (!strcmp(label, _idParamsA3GaB_) || !strcmp(label, _dpnlLabelA3GaB))
      {
        return _specDPnlA3GaB_;
      }

    /* d3gaa */
    if (!strcmp(label, _idParamsD3GaA_) || !strcmp(label, _dpnlLabelD3GaA))
      {
        return _specDPnlD3GaA_;
      }

    /* d3gab */
    if (!strcmp(label, _idParamsD3GaB_) || !strcmp(label, _dpnlLabelD3GaB))
      {
        return _specDPnlD3GaB_;
      }


    /**  Bias parameters  **/

    /* b1 */
    if (!strcmp(label, _idParamsB1_) || !strcmp(label, _dpnlLabelB1))
      {
        return _specDPnlB1_;
      }

    /* b2 */
    if (!strcmp(label, _idParamsB2_) || !strcmp(label, _dpnlLabelB2))
      {
        return _specDPnlB2_;
      }

    /* c2ga */
    if (!strcmp(label, _idParamsC2Ga_) || !strcmp(label, _dpnlLabelC2Ga))
      {
        return _specDPnlC2Ga_;
      }

    /* bgam3 */
    if (!strcmp(label, _idParamsBGam3_) || !strcmp(label, _dpnlLabelBGam3))
      {
        return _specDPnlBGam3_;
      }


    /**  Redshift space distortion parameters  **/

    /* f */
    if (!strcmp(label, _idParamsF_) || !strcmp(label, _dpnlLabelF))
      {
        return _specDPnlF_;
      }

    /* sigv */
    if (!strcmp(label, _idParamsSigv_) || !strcmp(label, _dpnlLabelF))
      {
        return _specDPnlSigv_;
      }


    /**  AP parameters  **/

    /* d */
    if (!strcmp(label, _idParamsD_) || !strcmp(label, _dpnlLabelD))
      {
        return _specDPnlD_;
      }

    /* h */
    if (!strcmp(label, _idParamsH_) || !strcmp(label, _dpnlLabelH))
      {
        return _specDPnlH_;
      }


    /**  Counter term parameters  **/

    /* c0 */
    if (!strcmp(label, _idParamsC0_) || !strcmp(label, _dpnlLabelC0))
      {
        return _specDPnlC0_;
      }

    /* c2 */
    if (!strcmp(label, _idParamsC2_) || !strcmp(label, _dpnlLabelC2))
      {
        return _specDPnlC2_;
      }

    /* c4 */
    if (!strcmp(label, _idParamsC4_) || !strcmp(label, _dpnlLabelC4))
      {
        return _specDPnlC4_;
      }


    /**  Shot noise parameters  **/

    /* Psn */
    if (!strcmp(label, _idParamsPsn_) || !strcmp(label, _dpnlLabelPsn))
      {
        return _specDPnlPsn_;
      }

    /* Bsn1 */
    if (!strcmp(label, _idParamsBsn1_) || !strcmp(label, _dpnlLabelBsn1))
      {
        return _specDPnlBsn1_;
      }

    /* Bsn2 */
    if (!strcmp(label, _idParamsBsn2_) || !strcmp(label, _dpnlLabelBsn2))
      {
        return _specDPnlBsn2_;
      }


    /**  Linear power spectrum  **/

    /* Pk */
    if (!strcmp(label, _idParamsPk_) || !strcmp(label, _dpnlLabelPk))
      {
        return _specDPnlPk_;
      }


    /**  Label not found  **/

    return _zeroFunc_;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Bispectrum     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   (Direct) Tree-Level Bispectrum   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_btr_tree(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the tree level bispectrum.

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * z2Kernel1 * pk1 * pk2
                                                       + z1Kernel1 * z1Kernel3 * z2Kernel2 * pk1 * pk3
                                                       + z1Kernel2 * z1Kernel3 * z2Kernel3 * pk2 * pk3);

    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_btr_sn(spec_arg_t *specArg, double *result)
{
    /*

        Shot noise for the bispectrum

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Interpolated power spectrum */

    /* Pnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl1 = (k1 != 0.) ? _specPnl_(specArg, NULL) : kern -> surv -> sn;

    /* Pnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl2 = (k2 != 0.) ? _specPnl_(specArg, NULL) : kern -> surv -> sn;

    /* Pnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl3 = (k3 != 0.) ? _specPnl_(specArg, NULL) : kern -> surv -> sn;

    /* Result (power spectrum has shot noise) */
    result[0] = kern -> surv -> sn * (pnl1 + pnl2 + pnl3) - 2. * kern -> surv -> sn*kern -> surv -> sn;
    result[1] = 0.;

    /* Reset variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_btr_full(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the full tree-level bispectrum plus shot noise

    */

    /* Must have specOrder = 3 */
    kern_t *kern = specArg -> kern;
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Tree level */
    double resultBtrTree[2];
    spec_btr_tree(specArg, resultBtrTree);

    /* Shot noise */
    double resultBtrSN[2];
    spec_btr_sn(specArg, resultBtrSN);

    /* Result */
    result[0] = resultBtrTree[0] + resultBtrSN[0];
    result[1] = sqrt(resultBtrTree[1]*resultBtrTree[1] + resultBtrSN[1]*resultBtrSN[1]);

    /* Reset specOrder */
    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*spec_btr(const char *label))(spec_arg_t*, double*)
{
    /*

        Return the requested function that calculates part (or full) of the tree level bispectrum

    */

    /* Tree-Level Bispectrum (excluding shot noise) */
    if (!strcmp(label, _btrLabelTree))
      {
        return spec_btr_tree;
      }

    /* Bispectrum Shot Noise */
    if (!strcmp(label, _btrLabelSn))
      {
        return spec_btr_sn;
      }

    /* Full Bispectrum */
    if (!strcmp(label, _btrLabelFull))
      {
        return spec_btr_full;
      }

    /* Label not found */
    return _spec_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------   (Indirect) Tree-Level Bispectrum   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_btr(void *var, void *params)
{
    /*

        Calculate the tree-level bispectrum (plus shot noise)

        Note: This indirect function is used as the initial function for 'double _specBtr_(void*, void*)'.
              The direct function 'int spec_btr(spec_arg_t*, double*)' is of different shape (but here they calculate the same thing in the same way)!

    */

    (void) params;

    double result[3];
    spec_btr_full(var, result);

    return result[0];
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------   (Direct) Tree-Level Bispectrum (Analytical) Derivatives   ---------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_a2ga(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(a^(2)_γ)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double dz2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_a2ga(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double dz2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_dz2_a2ga(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double dz2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_a2ga(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * dz2Kernel1 * pk1 * pk2
                                                       + z1Kernel1 * z1Kernel3 * dz2Kernel2 * pk1 * pk3
                                                       + z1Kernel2 * z1Kernel3 * dz2Kernel3 * pk2 * pk3);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double a2ga = _fidBTSTxA2Ga_(&deriv -> z, _fidParamsBTSTxA2Ga_);

        if (a2ga != 0.)
          {
            result[0] *= a2ga;
            result[1] *= fabs(a2ga);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsA2Ga_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsA2Ga_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsA2Ga_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_d2ga(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(d^(2)_γ)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double dz2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_d2ga(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double dz2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_dz2_d2ga(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double dz2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_d2ga(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * dz2Kernel1 * pk1 * pk2
                                                       + z1Kernel1 * z1Kernel3 * dz2Kernel2 * pk1 * pk3
                                                       + z1Kernel2 * z1Kernel3 * dz2Kernel3 * pk2 * pk3);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double d2ga = _fidBTSTxD2Ga_(&deriv -> z, _fidParamsBTSTxD2Ga_);

        if (d2ga != 0.)
          {
            result[0] *= d2ga;
            result[1] *= fabs(d2ga);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsD2Ga_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsD2Ga_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsD2Ga_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_a3gaa(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(a^(3)_γa)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsA3GaA_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsA3GaA_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsA3GaA_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


int spec_dbtr_a3gab(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(a^(3)_γb)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsA3GaB_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsA3GaB_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsA3GaB_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


int spec_dbtr_d3gaa(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(d^(3)_γa)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsD3GaA_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsD3GaA_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsD3GaA_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


int spec_dbtr_d3gab(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(d^(3)_γb)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsD3GaB_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsD3GaB_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsD3GaB_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_b1(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(b1)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);
    double dz1Kernel1 = kernels_dz1_b1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);
    double dz1Kernel2 = kernels_dz1_b1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);
    double dz1Kernel3 = kernels_dz1_b1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_b1(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_dz2_b1(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_b1(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * ((dz1Kernel1 * z1Kernel2 * z2Kernel1 + z1Kernel1 * dz1Kernel2 * z2Kernel1 + z1Kernel1 * z1Kernel2 * dz2Kernel1) * pk1 * pk2
                                                       + (dz1Kernel1 * z1Kernel3 * z2Kernel2 + z1Kernel1 * dz1Kernel3 * z2Kernel2 + z1Kernel1 * z1Kernel3 * dz2Kernel2) * pk1 * pk3
                                                       + (dz1Kernel2 * z1Kernel3 * z2Kernel3 + z1Kernel2 * dz1Kernel3 * z2Kernel3 + z1Kernel2 * z1Kernel3 * dz2Kernel3) * pk2 * pk3);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double b1 = _fidBiasxB1_(&deriv -> z, _fidParamsBiasxB1_);

        if (b1 != 0.)
          {
            result[0] *= b1;
            result[1] *= fabs(b1);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsB1_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsB1_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsB1_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_b2(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(b2)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double dz2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_b2(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double dz2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_dz2_b2(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double dz2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_b2(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * dz2Kernel1 * pk1 * pk2
                                                       + z1Kernel1 * z1Kernel3 * dz2Kernel2 * pk1 * pk3
                                                       + z1Kernel2 * z1Kernel3 * dz2Kernel3 * pk2 * pk3);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double b2 = _fidBiasxB2_(&deriv -> z, _fidParamsBiasxB2_);

        if (b2 != 0.)
          {
            result[0] *= b2;
            result[1] *= fabs(b2);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsB2_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsB2_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsB2_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_c2ga(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(c^(2)_γ)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double dz2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_c2ga(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double dz2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_dz2_c2ga(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double dz2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_c2ga(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * dz2Kernel1 * pk1 * pk2
                                                       + z1Kernel1 * z1Kernel3 * dz2Kernel2 * pk1 * pk3
                                                       + z1Kernel2 * z1Kernel3 * dz2Kernel3 * pk2 * pk3);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double c2ga = _fidBiasxC2Ga_(&deriv -> z, _fidParamsBiasxC2Ga_);

        if (c2ga != 0.)
          {
            result[0] *= c2ga;
            result[1] *= fabs(c2ga);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsC2Ga_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsC2Ga_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsC2Ga_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_bgam3(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(b_Γ3)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsBGam3_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsBGam3_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsBGam3_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_f(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(f)

            dBtr(z)/df(z') = (...) + δBtr(z)/δf(z')

        The ellipsis (...) contains the growth factor contribution:

            G(z) = exp(- int f(z') / (1+z') dz')

        the contribution is

                               {   - G(z) / (1+z') dz'   if z' <= z
            dG(z)/d(f(z')) =  {
                               {            0            else.

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* No contribution */
    if (kern -> z < deriv -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

    /* Growth contribution (dlogG/df) */
    double dGrowth = - (1. / (1. + deriv -> z)) * specArg -> dz;

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);
    double dz1Kernel1 = kernels_dz1_f(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);
    double dz1Kernel2 = kernels_dz1_f(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);
    double dz1Kernel3 = kernels_dz1_f(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_f(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_dz2_f(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_f(kern) : 0.;


    /* Result */

    result[0] = 0.;
    result[1] = 0.;

    if (kern -> z == specArg -> deriv -> z)
      {
        result[0] += 2. * pow(kern -> growth, 4.) * smooth * ((dz1Kernel1 * z1Kernel2 * z2Kernel1 + z1Kernel1 * dz1Kernel2 * z2Kernel1 + z1Kernel1 * z1Kernel2 * dz2Kernel1) * pk1 * pk2
                                                            + (dz1Kernel1 * z1Kernel3 * z2Kernel2 + z1Kernel1 * dz1Kernel3 * z2Kernel2 + z1Kernel1 * z1Kernel3 * dz2Kernel2) * pk1 * pk3
                                                            + (dz1Kernel2 * z1Kernel3 * z2Kernel3 + z1Kernel2 * dz1Kernel3 * z2Kernel3 + z1Kernel2 * z1Kernel3 * dz2Kernel3) * pk2 * pk3);
      }

    if (kern -> z <= specArg -> deriv -> z)
      {
        result[0] += 8. * pow(kern -> growth, 4.) * dGrowth * smooth * (z1Kernel1 * z1Kernel2 * z2Kernel1 * pk1 * pk2
                                                                       + z1Kernel1 * z1Kernel3 * z2Kernel2 * pk1 * pk3
                                                                       + z1Kernel2 * z1Kernel3 * z2Kernel3 * pk2 * pk3);
      }

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double f = _fidRSDxF_(&deriv -> z, _fidParamsRSDxF_);

        if (f != 0.)
          {
            result[0] *= f;
            result[1] *= fabs(f);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsF_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsF_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsF_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_sigv(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(sig_v)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double dsmooth = kernels_dsmooth_sigv(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * dsmooth * (z1Kernel1 * z1Kernel2 * z2Kernel1 * pk1 * pk2
                                                        + z1Kernel1 * z1Kernel3 * z2Kernel2 * pk1 * pk3
                                                        + z1Kernel2 * z1Kernel3 * z2Kernel3 * pk2 * pk3);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double sigv = _fidRSDxSigv_(&deriv -> z, _fidParamsRSDxSigv_);

        if (sigv != 0.)
          {
            result[0] *= sigv;
            result[1] *= fabs(sigv);
          }
      }


    /* Shote noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsSigv_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsSigv_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsSigv_))(specArg, NULL) : 0.;


    /* Result */

    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);


    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_d(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(d)

            dBtr/d(d) = -4 Btr + δBtr/δk_i δk_i/δd + δBtr/δmu_i δmu_i/δd + δBtr/δnu_ij δnu_ij/δd

        where d = Da / Da' is the angular diameter distance ratio which enters via the AP effect (assuming h = H/H' = 1, d = Da/Da' = 1 -> alpha = 1)
        such that:

            δk_i/δd = k_i (mu_i^2 - 1)
            δmu_i/δd = mu_i (1 - mu_i^2)
            δnu_ij/δd = mu_i (1 - mu_i^2)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double dpk1 = _fidDPk_(&k1, _fidParamsDPk_);

    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double dpk2 = _fidDPk_(&k2, _fidParamsDPk_);

    double pk3 = _fidPk_(&k3, _fidParamsPk_);
    double dpk3 = _fidDPk_(&k3, _fidParamsDPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);
    double dz1Kernel1 = kernels_dz1_mu(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);
    double dz1Kernel2 = kernels_dz1_mu(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);
    double dz1Kernel3 = kernels_dz1_mu(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel1Nu = (k1 != 0. && k2 != 0.) ? kernels_dz2_nu(kern) : 0.;

    double dz2Kernel1K1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel1Mu1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_mu(kern) : 0.;

    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_nu(kern, 0, 1, nu12);

    double dz2Kernel1K2 = (k1 != 0. && k2 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel1Mu2 = (k1 != 0. && k2 != 0.) ? kernels_dz2_mu(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel2Nu = (k1 != 0. && k3 != 0.) ? kernels_dz2_nu(kern) : 0.;

    double dz2Kernel2K1 = (k1 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel2Mu1 = (k1 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;

    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_nu(kern, 0, 1, nu13);

    double dz2Kernel2K3 = (k1 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel2Mu3 = (k1 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel3Nu = (k2 != 0. && k3 != 0.) ? kernels_dz2_nu(kern) : 0.;

    double dz2Kernel3K2 = (k2 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel3Mu2 = (k2 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;

    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu23);

    double dz2Kernel3K3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel3Mu3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;


    /* Result */

    /* Volume correction */
    result[0] = -8. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * z2Kernel1 * pk1 * pk2
                                                        + z1Kernel1 * z1Kernel3 * z2Kernel2 * pk1 * pk3
                                                        + z1Kernel2 * z1Kernel3 * z2Kernel3 * pk2 * pk3);


    /* k and mu derivatives */
    result[0] += 2. * pow(kern -> growth, 4.) * smooth *
                  (
                    (
                      z1Kernel1 * z1Kernel2 * (dz2Kernel1K1 * pk1 + z2Kernel1 * dpk1) * pk2
                    + z1Kernel1 * z1Kernel3 * (dz2Kernel2K1 * pk1 + z2Kernel2 * dpk1) * pk3
                    ) * k1 * (mu1 * mu1 - 1.)
                  +
                    (
                      z1Kernel1 * z1Kernel2 * (dz2Kernel1K2 * pk2 + z2Kernel2 * dpk2) * pk1
                    + z1Kernel2 * z1Kernel3 * (dz2Kernel3K2 * pk2 + z2Kernel3 * dpk2) * pk3
                    ) * k2 * (mu2 * mu2 - 1.)
                  +
                    (
                      z1Kernel1 * z1Kernel3 * (dz2Kernel2K3 * pk3 + z2Kernel2 * dpk3) * pk1
                    + z1Kernel2 * z1Kernel3 * (dz2Kernel3K3 * pk3 + z2Kernel3 * dpk3) * pk2
                    ) * k3 * (mu3 * mu3 - 1.)

                  +

                    (
                      z1Kernel1 * z1Kernel2 * dz2Kernel1Nu * pk1 * pk2
                    ) * (-nu12) * (mu1 * mu1 + mu2 * mu2) + 2. * mu1 * mu2
                  +
                    (
                      z1Kernel1 * z1Kernel3 * dz2Kernel2Nu * pk1 * pk3
                    ) * (-nu13) * (mu1 * mu1 + mu3 * mu3) + 2. * mu1 * mu3
                  +
                    (
                      z1Kernel2 * z1Kernel3 * dz2Kernel3Nu * pk2 * pk3
                    ) * (-nu23) * (mu2 * mu2 + mu3 * mu3) + 2. * mu2 * mu3

                  +

                    (
                      z1Kernel2 * (dz1Kernel1 * z2Kernel1 + z1Kernel1 * dz2Kernel1Mu1) * pk1 * pk2
                    + z1Kernel3 * (dz1Kernel1 * z2Kernel2 + z1Kernel1 * dz2Kernel2Mu1) * pk1 * pk3
                    ) * mu1 * (1. - mu1 * mu1)
                  +
                    (
                      z1Kernel1 * (dz1Kernel2 * z2Kernel1 + z1Kernel2 * dz2Kernel1Mu2) * pk1 * pk2
                    + z1Kernel3 * (dz1Kernel2 * z2Kernel3 + z1Kernel2 * dz2Kernel3Mu2) * pk2 * pk3
                    ) * mu2 * (1. - mu2 * mu2)
                  +
                    (
                      z1Kernel1 * (dz1Kernel3 * z2Kernel2 + z1Kernel3 * dz2Kernel2Mu3) * pk1 * pk3
                    + z1Kernel2 * (dz1Kernel3 * z2Kernel3 + z1Kernel3 * dz2Kernel3Mu3) * pk2 * pk3
                    ) * mu3 * (1. - mu3 * mu3)
                  );

    result[1] = 0.;


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsD_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsD_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsD_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


int spec_dbtr_h(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(h)

            dBtr/d(h) = 2 Btr + δBtr/δk_i δk_i/δh + δBtr/δmu_i δmu_i/δh

        where h = H / H' is the Hubble function ratio which enters via the AP effect (assuming h = H/H' = 1, d = Da/Da' = 1 -> alpha = 1)
        such that:

            δk_i/δh = k_i mu_i^2
            δmu_i/δh = mu_i (1 - mu_i^2)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double dpk1 = _fidDPk_(&k1, _fidParamsDPk_);

    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double dpk2 = _fidDPk_(&k2, _fidParamsDPk_);

    double pk3 = _fidPk_(&k3, _fidParamsPk_);
    double dpk3 = _fidDPk_(&k3, _fidParamsDPk_);


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);
    double dz1Kernel1 = kernels_dz1_mu(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);
    double dz1Kernel2 = kernels_dz1_mu(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);
    double dz1Kernel3 = kernels_dz1_mu(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel1Nu = (k1 != 0. && k2 != 0.) ? kernels_dz2_nu(kern) : 0.;

    double dz2Kernel1K1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel1Mu1 = (k1 != 0. && k2 != 0.) ? kernels_dz2_mu(kern) : 0.;

    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_nu(kern, 0, 1, nu12);

    double dz2Kernel1K2 = (k1 != 0. && k2 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel1Mu2 = (k1 != 0. && k2 != 0.) ? kernels_dz2_mu(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel2Nu = (k1 != 0. && k3 != 0.) ? kernels_dz2_nu(kern) : 0.;

    double dz2Kernel2K1 = (k1 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel2Mu1 = (k1 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;

    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_nu(kern, 0, 1, nu13);

    double dz2Kernel2K3 = (k1 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel2Mu3 = (k1 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;
    double dz2Kernel3Nu = (k2 != 0. && k3 != 0.) ? kernels_dz2_nu(kern) : 0.;

    double dz2Kernel3K2 = (k2 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel3Mu2 = (k2 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;

    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu23);

    double dz2Kernel3K3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_k(kern) : 0.;
    double dz2Kernel3Mu3 = (k2 != 0. && k3 != 0.) ? kernels_dz2_mu(kern) : 0.;


    /* Result */

    /* Volume correction */
    result[0] = 4. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * z2Kernel1 * pk1 * pk2
                                                       + z1Kernel1 * z1Kernel3 * z2Kernel2 * pk1 * pk3
                                                       + z1Kernel2 * z1Kernel3 * z2Kernel3 * pk2 * pk3);


    /* k and mu derivatives */
    result[0] += 2. * pow(kern -> growth, 4.) * smooth *
                  (
                    (
                      z1Kernel1 * z1Kernel2 * (dz2Kernel1K1 * pk1 + z2Kernel1 * dpk1) * pk2
                    + z1Kernel1 * z1Kernel3 * (dz2Kernel2K1 * pk1 + z2Kernel2 * dpk1) * pk3
                    ) * k1 * mu1 * mu1
                  +
                    (
                      z1Kernel1 * z1Kernel2 * (dz2Kernel1K2 * pk2 + z2Kernel2 * dpk2) * pk1
                    + z1Kernel2 * z1Kernel3 * (dz2Kernel3K2 * pk2 + z2Kernel3 * dpk2) * pk3
                    ) * k2 * mu2 * mu2
                  +
                    (
                      z1Kernel1 * z1Kernel3 * (dz2Kernel2K3 * pk3 + z2Kernel2 * dpk3) * pk1
                    + z1Kernel2 * z1Kernel3 * (dz2Kernel3K3 * pk3 + z2Kernel3 * dpk3) * pk2
                    ) * k3 * mu3 * mu3

                  +

                    (
                      z1Kernel1 * z1Kernel2 * dz2Kernel1Nu * pk1 * pk2
                    ) * (-nu12) * (mu1 * mu1 + mu2 * mu2) + 2. * mu1 * mu2
                  +
                    (
                      z1Kernel1 * z1Kernel3 * dz2Kernel2Nu * pk1 * pk3
                    ) * (-nu13) * (mu1 * mu1 + mu3 * mu3) + 2. * mu1 * mu3
                  +
                    (
                      z1Kernel2 * z1Kernel3 * dz2Kernel3Nu * pk2 * pk3
                    ) * (-nu23) * (mu2 * mu2 + mu3 * mu3) + 2. * mu2 * mu3

                  +

                    (
                      z1Kernel2 * (dz1Kernel1 * z2Kernel1 + z1Kernel1 * dz2Kernel1Mu1) * pk1 * pk2
                    + z1Kernel3 * (dz1Kernel1 * z2Kernel2 + z1Kernel1 * dz2Kernel2Mu1) * pk1 * pk3
                    ) * mu1 * (1. - mu1 * mu1)
                  +
                    (
                      z1Kernel1 * (dz1Kernel2 * z2Kernel1 + z1Kernel2 * dz2Kernel1Mu2) * pk1 * pk2
                    + z1Kernel3 * (dz1Kernel2 * z2Kernel3 + z1Kernel2 * dz2Kernel3Mu2) * pk2 * pk3
                    ) * mu2 * (1. - mu2 * mu2)
                  +
                    (
                      z1Kernel1 * (dz1Kernel3 * z2Kernel2 + z1Kernel3 * dz2Kernel2Mu3) * pk1 * pk3
                    + z1Kernel2 * (dz1Kernel3 * z2Kernel3 + z1Kernel3 * dz2Kernel3Mu3) * pk2 * pk3
                    ) * mu3 * (1. - mu3 * mu3)
                  );

    result[1] = 0.;


    /* Result */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsH_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsH_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsH_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_c0(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(c0)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsC0_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsC0_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsC0_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


int spec_dbtr_c2(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(c2)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsC2_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsC2_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsC2_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


int spec_dbtr_c4(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(c4)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }


    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsC4_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsC4_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsC4_))(specArg, NULL) : 0.;

    /* Result */
    result[0] = kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_psn(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(Psn)

    */

    (void) specArg;

    /* Result (no Psn only Bsn for bispectrum) */
    result[0] = 0.;
    result[1] = 0.;

    return 0;
}


int spec_dbtr_bsn1(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(Bsn1)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);

    double nu12 = kernels_qget_nu(kern, 0, 1);


    /* Shot noise */

    /* Pnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl1 = (k1 != 0.) ? _specPnl_(specArg, NULL) : 0.;

    /* Pnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl2 = (k2 != 0.) ? _specPnl_(specArg, NULL) : 0.;

    /* Pnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl3 = (k3 != 0.) ? _specPnl_(specArg, NULL) : 0.;

    /* Result (must remove the shot noise from pnl) */
    result[0] = kern -> surv -> sn * (pnl1 + pnl2 + pnl3) - 3. * pow(kern -> surv -> sn, 2.);
    result[1] = 0.;

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    return 0;
}


int spec_dbtr_bsn2(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(Bsn2)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Only diagonal terms */
    if (deriv -> z != kern -> z)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Result */
    result[0] = pow(kern -> surv -> sn, 2.);
    result[1] = 0.;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_dbtr_pk(spec_arg_t *specArg, double *result)
{
    /*

        Calculate dBtr/d(Pk)

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Deriv struct */
    spec_deriv_t *deriv = specArg -> deriv;

    /* Must have specOrder = 3 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 3;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

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

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);

    double dpk12 = 0.;
    double dpk13 = 0.;
    double dpk23 = 0.;

    if (fabs(specArg -> deriv -> k - k1) < __ABSTOL__)
      {
        dpk12 += pk2;
        dpk13 += pk3;
      }

    if (fabs(specArg -> deriv -> k - k2) < __ABSTOL__)
      {
        dpk12 += pk1;
        dpk23 += pk3;
      }

    if (fabs(specArg -> deriv -> k - k3) < __ABSTOL__)
      {
        dpk13 += pk1;
        dpk23 += pk2;
      }


    /* Kernels */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z2(k1_, k2_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    double z2Kernel1 = (k1 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k1_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu13);

    double z2Kernel2 = (k1 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k2_, k3_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, nu23);

    double z2Kernel3 = (k2 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;


    /* Result */

    result[0] = 2. * pow(kern -> growth, 4.) * smooth * (z1Kernel1 * z1Kernel2 * z2Kernel1 * dpk12
                                                       + z1Kernel1 * z1Kernel3 * z2Kernel2 * dpk13
                                                       + z1Kernel2 * z1Kernel3 * z2Kernel3 * dpk23);

    result[1] = 0.;

    /* Logarithmic derivative */
    if (deriv -> log)
      {
        double pk = _fidPk_(&deriv -> k, _fidParamsPk_);

        if (pk != 0.)
          {
            result[0] *= pk;
            result[1] *= fabs(pk);
          }
      }


    /* Shot noise */

    /* DPnl(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl1 = (k1 != 0.) ? (_specDPnl_(_idParamsPk_))(specArg, NULL) : 0.;

    /* DPnl(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl2 = (k2 != 0.) ? (_specDPnl_(_idParamsPk_))(specArg, NULL) : 0.;

    /* DPnl(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double dpnl3 = (k3 != 0.) ? (_specDPnl_(_idParamsPk_))(specArg, NULL) : 0.;

    /* Result */
    result[0] += kern -> surv -> sn * (dpnl1 + dpnl2 + dpnl3);

    /* Reset the variables */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, nu12);

    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*spec_dbtr(const char *label))(spec_arg_t*, double*)
{
    /*

        Get the tree-level bispectrum derivative function w.r.t. to parameter with given label.

        TODO: Might dereference NULL pointer (_dbtrLabel...) if ini functions are not called...

    */


    /**  Bootstrap parameters  **/

    /* a2ga */
    if (!strcmp(label, _idParamsA2Ga_) || !strcmp(label, _dbtrLabelA2Ga))
      {
        return spec_dbtr_a2ga;
      }

    /* d2ga */
    if (!strcmp(label, _idParamsD2Ga_) || !strcmp(label, _dbtrLabelD2Ga))
      {
        return spec_dbtr_d2ga;
      }

    /* a3gaa */
    if (!strcmp(label, _idParamsA3GaA_) || !strcmp(label, _dbtrLabelA3GaA))
      {
        return spec_dbtr_a3gaa;
      }

    /* a3gab */
    if (!strcmp(label, _idParamsA3GaB_) || !strcmp(label, _dbtrLabelA3GaB))
      {
        return spec_dbtr_a3gab;
      }

    /* d3gaa */
    if (!strcmp(label, _idParamsD3GaA_) || !strcmp(label, _dbtrLabelD3GaA))
      {
        return spec_dbtr_d3gaa;
      }

    /* d3gab */
    if (!strcmp(label, _idParamsD3GaB_) || !strcmp(label, _dbtrLabelD3GaB))
      {
        return spec_dbtr_d3gab;
      }


    /**  Bias parameters  **/

    /* b1 */
    if (!strcmp(label, _idParamsB1_) || !strcmp(label, _dbtrLabelB1))
      {
        return spec_dbtr_b1;
      }

    /* b2 */
    if (!strcmp(label, _idParamsB2_) || !strcmp(label, _dbtrLabelB2))
      {
        return spec_dbtr_b2;
      }

    /* c2ga */
    if (!strcmp(label, _idParamsC2Ga_) || !strcmp(label, _dbtrLabelC2Ga))
      {
        return spec_dbtr_c2ga;
      }

    /* bgam3 */
    if (!strcmp(label, _idParamsBGam3_) || !strcmp(label, _dbtrLabelBGam3))
      {
        return spec_dbtr_bgam3;
      }


    /**  Redshift space distortion parameters  **/

    /* f */
    if (!strcmp(label, _idParamsF_) || !strcmp(label, _dbtrLabelF))
      {
        return spec_dbtr_f;
      }

    /* sigv */
    if (!strcmp(label, _idParamsSigv_) || !strcmp(label, _dbtrLabelSigv))
      {
        return spec_dbtr_sigv;
      }


    /**  AP parameters  **/

    /* d */
    if (!strcmp(label, _idParamsD_) || !strcmp(label, _dbtrLabelD))
      {
        return spec_dbtr_d;
      }

    /* h */
    if (!strcmp(label, _idParamsH_) || !strcmp(label, _dbtrLabelH))
      {
        return spec_dbtr_h;
      }


    /**  Counter term parameters  **/

    /* c0 */
    if (!strcmp(label, _idParamsC0_) || !strcmp(label, _dbtrLabelC0))
      {
        return spec_dbtr_c0;
      }

    /* c2 */
    if (!strcmp(label, _idParamsC2_) || !strcmp(label, _dbtrLabelC2))
      {
        return spec_dbtr_c2;
      }

    /* c4 */
    if (!strcmp(label, _idParamsC4_) || !strcmp(label, _dbtrLabelC4))
      {
        return spec_dbtr_c4;
      }


    /**  Shot noise parameters  **/

    /* Psn */
    if (!strcmp(label, _idParamsPsn_) || !strcmp(label, _dbtrLabelPsn))
      {
        return spec_dbtr_psn;
      }

    /* Bsn1 */
    if (!strcmp(label, _idParamsBsn1_) || !strcmp(label, _dbtrLabelBsn1))
      {
        return spec_dbtr_bsn1;
      }

    /* Bsn2 */
    if (!strcmp(label, _idParamsBsn2_) || !strcmp(label, _dbtrLabelBsn2))
      {
        return spec_dbtr_bsn2;
      }


    /**  Linear power spectrum  **/

    /* Pk */
    if (!strcmp(label, _idParamsPk_) || !strcmp(label, _dbtrLabelPk))
      {
        return spec_dbtr_pk;
      }


    /**  Label not found  **/

    return _spec_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------   (Indirect) Tree-Level Bispectrum (Analytical) Derivatives   ---------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_a2ga(void *var, void *params)
{
    /*

        Calculate dBtr/d(a^(2)_γ)

    */

    (void) params;

    double result[3];
    spec_dbtr_a2ga((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_d2ga(void *var, void *params)
{
    /*

        Calculate dBtr/d(d^(2)_γ)

    */

    (void) params;

    double result[3];
    spec_dbtr_d2ga((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_a3gaa(void *var, void *params)
{
    /*

        Calculate dBtr/d(a^(3)_γa)

    */

    (void) params;

    double result[3];
    spec_dbtr_a3gaa((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_a3gab(void *var, void *params)
{
    /*

        Calculate dBtr/d(a^(3)_γb)

    */

    (void) params;

    double result[3];
    spec_dbtr_a3gab((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_d3gaa(void *var, void *params)
{
    /*

        Calculate dBtr/d(d^(3)_γa)

    */

    (void) params;

    double result[3];
    spec_dbtr_d3gaa((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_d3gab(void *var, void *params)
{
    /*

        Calculate dBtr/d(d^(3)_γb)

    */

    (void) params;

    double result[3];
    spec_dbtr_d3gab((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_b1(void *var, void *params)
{
    /*

        Calculate dBtr/d(b1)

    */

    (void) params;

    double result[3];
    spec_dbtr_b1((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_b2(void *var, void *params)
{
    /*

        Calculate dBtr/d(b2)

    */

    (void) params;

    double result[3];
    spec_dbtr_b2((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_c2ga(void *var, void *params)
{
    /*

        Calculate dBtr/d(c^(2)_γ)

    */

    (void) params;

    double result[3];
    spec_dbtr_c2ga((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_bgam3(void *var, void *params)
{
    /*

        Calculate dBtr/d(b_Γ3)

    */

    (void) params;

    double result[3];
    spec_dbtr_bgam3((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_f(void *var, void *params)
{
    /*

        Calculate dBtr/d(f)

    */

    (void) params;

    double result[3];
    spec_dbtr_f((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_sigv(void *var, void *params)
{
    /*

        Calculate dBtr/d(sig_v)

    */

    (void) params;

    double result[3];
    spec_dbtr_sigv((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_d(void *var, void *params)
{
    /*

        Calculate dBtr/d(d)

    */

    (void) params;

    double result[3];
    spec_dbtr_d((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_h(void *var, void *params)
{
    /*

        Calculate dBtr/d(h)

    */

    (void) params;

    double result[3];
    spec_dbtr_h((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_c0(void *var, void *params)
{
    /*

        Calculate dBtr/d(c0)

    */

    (void) params;

    double result[3];
    spec_dbtr_c0((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_c2(void *var, void *params)
{
    /*

        Calculate dBtr/d(c2)

    */

    (void) params;

    double result[3];
    spec_dbtr_c2((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_c4(void *var, void *params)
{
    /*

        Calculate dBtr/d(c4)

    */

    (void) params;

    double result[3];
    spec_dbtr_c4((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_psn(void *var, void *params)
{
    /*

        Calculate dBtr/d(Psn)

    */

    (void) params;

    double result[3];
    spec_dbtr_psn((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_bsn1(void *var, void *params)
{
    /*

        Calculate dBtr/d(Bsn1)

    */

    (void) params;

    double result[3];
    spec_dbtr_bsn1((spec_arg_t*) var, result);

    return result[0];
}


static double _spec_dbtr_bsn2(void *var, void *params)
{
    /*

        Calculate dBtr/d(Bsn2)

    */

    (void) params;

    double result[3];
    spec_dbtr_bsn2((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_dbtr_pk(void *var, void *params)
{
    /*

        Calculate dBtr/d(Pk)

    */

    (void) params;

    double result[3];
    spec_dbtr_pk((spec_arg_t*) var, result);

    return result[0];
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double (*_spec_dbtr(const char *label))(void*, void*)
{
    /*

        Get the tree-level bispectrum derivative function w.r.t. to parameter with given label.

        TODO: Might dereference NULL pointer (_dbtrLabel...) if ini functions are not called...

    */


    /**  Bootstrap parameters  **/

    /* a2ga */
    if (!strcmp(label, _idParamsA2Ga_) || !strcmp(label, _dbtrLabelA2Ga))
      {
        return _specDBtrA2Ga_;
      }

    /* d2ga */
    if (!strcmp(label, _idParamsD2Ga_) || !strcmp(label, _dbtrLabelD2Ga))
      {
        return _specDBtrD2Ga_;
      }

    /* a3gaa */
    if (!strcmp(label, _idParamsA3GaA_) || !strcmp(label, _dbtrLabelA3GaA))
      {
        return _specDBtrA3GaA_;
      }

    /* a3gab */
    if (!strcmp(label, _idParamsA3GaB_) || !strcmp(label, _dbtrLabelA3GaB))
      {
        return _specDBtrA3GaB_;
      }

    /* d3gaa */
    if (!strcmp(label, _idParamsD3GaA_) || !strcmp(label, _dbtrLabelD3GaA))
      {
        return _specDBtrD3GaA_;
      }

    /* d3gab */
    if (!strcmp(label, _idParamsD3GaB_) || !strcmp(label, _dbtrLabelD3GaB))
      {
        return _specDBtrD3GaB_;
      }


    /**  Bias parameters  **/

    /* b1 */
    if (!strcmp(label, _idParamsB1_) || !strcmp(label, _dbtrLabelB1))
      {
        return _specDBtrB1_;
      }

    /* b2 */
    if (!strcmp(label, _idParamsB2_) || !strcmp(label, _dbtrLabelB2))
      {
        return _specDBtrB2_;
      }

    /* c2ga */
    if (!strcmp(label, _idParamsC2Ga_) || !strcmp(label, _dbtrLabelC2Ga))
      {
        return _specDBtrC2Ga_;
      }

    /* bgam3 */
    if (!strcmp(label, _idParamsBGam3_) || !strcmp(label, _dbtrLabelBGam3))
      {
        return _specDBtrBGam3_;
      }


    /**  Redshift space distortion parameters  **/

    /* f */
    if (!strcmp(label, _idParamsF_) || !strcmp(label, _dbtrLabelF))
      {
        return _specDBtrF_;
      }

    /* sigv */
    if (!strcmp(label, _idParamsSigv_) || !strcmp(label, _dbtrLabelF))
      {
        return _specDBtrSigv_;
      }


    /**  AP parameters  **/

    /* d */
    if (!strcmp(label, _idParamsD_) || !strcmp(label, _dbtrLabelD))
      {
        return _specDBtrD_;
      }

    /* h */
    if (!strcmp(label, _idParamsH_) || !strcmp(label, _dbtrLabelH))
      {
        return _specDBtrH_;
      }


    /**  Counter term parameters  **/

    /* c0 */
    if (!strcmp(label, _idParamsC0_) || !strcmp(label, _dbtrLabelC0))
      {
        return _specDBtrC0_;
      }

    /* c2 */
    if (!strcmp(label, _idParamsC2_) || !strcmp(label, _dbtrLabelC2))
      {
        return _specDBtrC2_;
      }

    /* c4 */
    if (!strcmp(label, _idParamsC4_) || !strcmp(label, _dbtrLabelC4))
      {
        return _specDBtrC4_;
      }


    /**  Shot noise parameters  **/

    /* Psn */
    if (!strcmp(label, _idParamsPsn_) || !strcmp(label, _dbtrLabelPsn))
      {
        return _specDBtrPsn_;
      }

    /* Bsn1 */
    if (!strcmp(label, _idParamsBsn1_) || !strcmp(label, _dbtrLabelBsn1))
      {
        return _specDBtrBsn1_;
      }

    /* Bsn2 */
    if (!strcmp(label, _idParamsBsn2_) || !strcmp(label, _dbtrLabelBsn2))
      {
        return _specDBtrBsn2_;
      }


    /**  Linear power spectrum  **/

    /* Pk */
    if (!strcmp(label, _idParamsPk_) || !strcmp(label, _dbtrLabelPk))
      {
        return _specDBtrPk_;
      }


    /**  Label not found  **/

    return _zeroFunc_;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Trispectrum     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   (Direct) Tree-Level Trispectrum   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int spec_ttr_tree(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the tree level trispectrum.

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Must have specOrder = 4 */
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 4;

    /* Smoothing */
    double smooth = kernels_smooth(kern);

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);
    double k4 = kernels_qget_k(kern, 3);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);
    double mu4 = kernels_qget_mu(kern, 3);

    double nu12 = kernels_qget_nu(kern, 0, 1);
    double nu13 = kernels_qget_nu(kern, 0, 2);
    double nu14 = kernels_qget_nu(kern, 0, 3);
    double nu23 = kernels_qget_nu(kern, 1, 2);
    double nu24 = kernels_qget_nu(kern, 1, 3);
    double nu34 = kernels_qget_nu(kern, 2, 3);

    /* kij = |ki_ + kj_| = sqrt( ki^2 + kj^2 + 2 ki kj nuij ) */
    double k12 = sqrt(k1*k1 + k2*k2 + 2*k1*k2*nu12);
    double k13 = sqrt(k1*k1 + k3*k3 + 2*k1*k3*nu13);
    double k14 = sqrt(k1*k1 + k4*k4 + 2*k1*k4*nu14);
//    double k23 = k14;
//    double k24 = k13;
//    double k34 = k12;

    /* muij = (ki_ + kj_).s_ / kij = (mui ki + muj kj) / kij */
    double mu12 = (k12 != 0.) ? (mu1*k1 + mu2*k2) / k12 : 0.;
    double mu13 = (k13 != 0.) ? (mu1*k1 + mu3*k3) / k13 : 0.;
    double mu14 = (k14 != 0.) ? (mu1*k1 + mu4*k4) / k14 : 0.;
//    double mu23 = - mu14;
//    double mu24 = - mu13;
//    double mu34 = - mu12;

    /* nuij_k = (ki_ + kj_).kk_ / (kij kk) = (nuik ki + nujk kj) / kij */
    double nu12_1 = (k12 != 0.) ? (nu12 * k2 + k1) / k12 : 0.;
    double nu12_2 = (k12 != 0.) ? (nu12 * k1 + k2) / k12 : 0.;

    double nu13_1 = (k13 != 0.) ? (nu13 * k3 + k1) / k13 : 0.;
    double nu13_3 = (k13 != 0.) ? (nu13 * k1 + k3) / k13 : 0.;

    double nu14_1 = (k14 != 0.) ? (nu14 * k4 + k1) / k14 : 0.;
    double nu14_4 = (k14 != 0.) ? (nu14 * k1 + k4) / k14 : 0.;

    double nu23_2 = (k14 != 0.) ? (nu23 * k3 + k2) / k14 : 0.;
    double nu23_3 = (k14 != 0.) ? (nu23 * k2 + k3) / k14 : 0.;

    double nu24_2 = (k13 != 0.) ? (nu24 * k4 + k2) / k13 : 0.;
    double nu24_4 = (k13 != 0.) ? (nu24 * k2 + k4) / k13 : 0.;

    double nu34_3 = (k12 != 0.) ? (nu34 * k4 + k3) / k12 : 0.;
    double nu34_4 = (k12 != 0.) ? (nu34 * k3 + k4) / k12 : 0.;

    /* Linear power spectra */
    double pk1 = _fidPk_(&k1, _fidParamsPk_);
    double pk2 = _fidPk_(&k2, _fidParamsPk_);
    double pk3 = _fidPk_(&k3, _fidParamsPk_);
    double pk4 = _fidPk_(&k4, _fidParamsPk_);

    double pk12 = _fidPk_(&k12, _fidParamsPk_);
    double pk13 = _fidPk_(&k13, _fidParamsPk_);
    double pk14 = _fidPk_(&k14, _fidParamsPk_);


    /**  Kernels  **/

    /* Z1 */

    /* Z1(k1_) */
    kernels_qset_mu(kern, 0, mu1);
    double z1Kernel1 = kernels_z1(kern);

    /* Z1(k2_) */
    kernels_qset_mu(kern, 0, mu2);
    double z1Kernel2 = kernels_z1(kern);

    /* Z1(k3_) */
    kernels_qset_mu(kern, 0, mu3);
    double z1Kernel3 = kernels_z1(kern);

    /* Z1(k4_) */
    kernels_qset_mu(kern, 0, mu4);
    double z1Kernel4 = kernels_z1(kern);


    /* Z2 */

    /* Z2(k12_, -k2_) */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu12);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -nu12_2);

    double z2Kernel12_2 = (k12 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k12_, -k1_) */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu12);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -nu12_1);

    double z2Kernel12_1 = (k12 != 0. && k1 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k13_, -k3_) */
    kernels_qset_k(kern, 0, k13);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu13);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -nu13_3);

    double z2Kernel13_3 = (k13 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k13_, -k1_) */
    kernels_qset_k(kern, 0, k13);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu13);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -nu13_1);

    double z2Kernel13_1 = (k13 != 0. && k1 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k14_, -k4_) */
    kernels_qset_k(kern, 0, k14);
    kernels_qset_k(kern, 1, k4);
    kernels_qset_mu(kern, 0, mu14);
    kernels_qset_mu(kern, 1, -mu4);
    kernels_qset_nu(kern, 0, 1, -nu14_4);

    double z2Kernel14_4 = (k14 != 0. && k4 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k14_, -k1_) */
    kernels_qset_k(kern, 0, k14);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu14);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -nu14_1);

    double z2Kernel14_1 = (k14 != 0. && k1 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k23_, -k3_) = Z2(k14_, k3_) */
    kernels_qset_k(kern, 0, k14);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu14);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, -nu23_3);

    double z2Kernel23_3 = (k14 != 0. && k3 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k23_, -k2_) = Z2(k14_, k2_) */
    kernels_qset_k(kern, 0, k14);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu14);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, -nu23_2);

    double z2Kernel23_2 = (k14 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k24_, -k4_) = Z2(k13_, k4_) */
    kernels_qset_k(kern, 0, k13);
    kernels_qset_k(kern, 1, k4);
    kernels_qset_mu(kern, 0, mu13);
    kernels_qset_mu(kern, 1, mu4);
    kernels_qset_nu(kern, 0, 1, -nu24_4);

    double z2Kernel24_4 = (k13 != 0. && k4 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k24_, -k2_) = Z2(k13_, k2_) */
    kernels_qset_k(kern, 0, k13);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu13);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_nu(kern, 0, 1, -nu24_2);

    double z2Kernel24_2 = (k13 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k34_, -k4_) = Z2(k12_, k4_) */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k4);
    kernels_qset_mu(kern, 0, mu12);
    kernels_qset_mu(kern, 1, mu4);
    kernels_qset_nu(kern, 0, 1, -nu34_4);

    double z2Kernel34_4 = (k12 != 0. && k4 != 0.) ? kernels_z2(kern) : 0.;

    /* Z2(k34_, -k3_) = Z2(k12_, k3_) */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu12);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_nu(kern, 0, 1, -nu34_3);

    double z2Kernel34_3 = (k12 != 0. && k2 != 0.) ? kernels_z2(kern) : 0.;


    /* Z3 */

    /* Z3(k1_, k2_, k3_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k3);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu3);

    kernels_qset_nu(kern, 0, 1, nu12);
    kernels_qset_nu(kern, 0, 2, nu13);
    kernels_qset_nu(kern, 1, 2, nu23);

    double z3Kernel123 = (k1 != 0. && k2 != 0. && k3 != 0.) ? kernels_z3(kern) : 0.;

    /* Z3(k1_, k2_, k4_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k4);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu4);

    kernels_qset_nu(kern, 0, 1, nu12);
    kernels_qset_nu(kern, 0, 2, nu14);
    kernels_qset_nu(kern, 1, 2, nu24);

    double z3Kernel124 = (k1 != 0. && k2 != 0. && k4 != 0.) ? kernels_z3(kern) : 0.;

    /* Z3(k1_, k3_, k4_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_k(kern, 2, k4);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_mu(kern, 2, mu4);

    kernels_qset_nu(kern, 0, 1, nu13);
    kernels_qset_nu(kern, 0, 2, nu14);
    kernels_qset_nu(kern, 1, 2, nu34);

    double z3Kernel134 = (k1 != 0. && k3 != 0. && k4 != 0.) ? kernels_z3(kern) : 0.;

    /* Z3(k2_, k3_, k4_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_k(kern, 2, k4);

    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_mu(kern, 2, mu4);

    kernels_qset_nu(kern, 0, 1, nu23);
    kernels_qset_nu(kern, 0, 2, nu24);
    kernels_qset_nu(kern, 1, 2, nu34);

    double z3Kernel234 = (k2 != 0. && k3 != 0. && k4 != 0.) ? kernels_z3(kern) : 0.;


    /**  Result  **/

    /* T2211 */
    double tri2211 = 4. * (z2Kernel13_3 * z2Kernel24_4 * pk13 + z2Kernel14_4 * z2Kernel23_3 * pk14) * z1Kernel3 * z1Kernel4 *  pk3 * pk4
                   + 4. * (z2Kernel12_2 * z2Kernel34_4 * pk12 + z2Kernel14_4 * z2Kernel23_2 * pk14) * z1Kernel2 * z1Kernel4 *  pk2 * pk4
                   + 4. * (z2Kernel13_3 * z2Kernel24_2 * pk13 + z2Kernel12_2 * z2Kernel34_3 * pk12) * z1Kernel2 * z1Kernel3 *  pk2 * pk3
                   + 4. * (z2Kernel13_1 * z2Kernel24_4 * pk13 + z2Kernel34_4 * z2Kernel12_1 * pk12) * z1Kernel1 * z1Kernel4 *  pk1 * pk4
                   + 4. * (z2Kernel34_3 * z2Kernel12_1 * pk12 + z2Kernel14_1 * z2Kernel23_3 * pk14) * z1Kernel1 * z1Kernel3 *  pk1 * pk3
                   + 4. * (z2Kernel13_1 * z2Kernel24_2 * pk13 + z2Kernel14_1 * z2Kernel23_2 * pk14) * z1Kernel1 * z1Kernel2 *  pk1 * pk2;

    /* T3111 */
    double tri3111 = 6. * z3Kernel123 * z1Kernel1 * z1Kernel2 * z1Kernel3 * pk1 * pk2 * pk3
                   + 6. * z3Kernel124 * z1Kernel1 * z1Kernel2 * z1Kernel4 * pk1 * pk2 * pk4
                   + 6. * z3Kernel134 * z1Kernel1 * z1Kernel3 * z1Kernel4 * pk1 * pk3 * pk4
                   + 6. * z3Kernel234 * z1Kernel2 * z1Kernel3 * z1Kernel4 * pk2 * pk3 * pk4;

    /* Combine the results */
    result[0] = pow(kern -> growth, 6.) * smooth * (tri2211 + tri3111);
    result[1] = 0.;

    /* Reset the variables */
    kern -> specOrder = specOrder;

    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k3);

    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu3);

    kernels_qset_nu(kern, 0, 1, nu12);
    kernels_qset_nu(kern, 0, 2, nu13);
    kernels_qset_nu(kern, 1, 2, nu23);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_ttr_sn(spec_arg_t *specArg, double *result)
{
    /*

        Shot noise for the trispectrum

    */

    /* Kern struct */
    kern_t *kern = specArg -> kern;

    /* Variables */
    double k1 = kernels_qget_k(kern, 0);
    double k2 = kernels_qget_k(kern, 1);
    double k3 = kernels_qget_k(kern, 2);
    double k4 = kernels_qget_k(kern, 3);

    double mu1 = kernels_qget_mu(kern, 0);
    double mu2 = kernels_qget_mu(kern, 1);
    double mu3 = kernels_qget_mu(kern, 2);
    double mu4 = kernels_qget_mu(kern, 3);

    double nu12 = kernels_qget_nu(kern, 0, 1);
    double nu13 = kernels_qget_nu(kern, 0, 2);
    double nu14 = kernels_qget_nu(kern, 0, 3);
    double nu23 = kernels_qget_nu(kern, 1, 2);
    double nu24 = kernels_qget_nu(kern, 1, 3);
    double nu34 = kernels_qget_nu(kern, 2, 3);

    /* kij = |ki_ + kj_| = sqrt( ki^2 + kj^2 + 2 ki kj nuij ) */
    double k12 = sqrt(k1*k1 + k2*k2 + 2*k1*k2*nu12);
    double k13 = sqrt(k1*k1 + k3*k3 + 2*k1*k3*nu13);
    double k14 = sqrt(k1*k1 + k4*k4 + 2*k1*k4*nu14);

//    double k23 = k14;
//    double k24 = k13;
//    double k34 = k12;

    /* muij = (ki_ + kj_).s_ / kij = (mui ki + muj kj) / kij */
    double mu12 = (k12 != 0.) ? (mu1*k1 + mu2*k2) / k12 : 0.;
    double mu13 = (k13 != 0.) ? (mu1*k1 + mu3*k3) / k13 : 0.;
    double mu14 = (k14 != 0.) ? (mu1*k1 + mu4*k4) / k14 : 0.;
//    double mu23 = -mu14;
//    double mu24 = -mu13;
//    double mu34 = -mu12;

    /* nuij_k = (ki_ + kj_).kk_ / (kij kk) = (nuik ki + nujk kj) / kij */
    double nu12_3 = (k12 != 0.) ? (nu13 * k1 + nu23 * k2) / k12 : 0.;
    double nu12_4 = (k12 != 0.) ? (nu14 * k1 + nu24 * k2) / k12 : 0.;
    double nu13_2 = (k13 != 0.) ? (nu12 * k1 + nu23 * k3) / k13 : 0.;
    double nu13_4 = (k13 != 0.) ? (nu14 * k1 + nu34 * k3) / k13 : 0.;
    double nu14_2 = (k14 != 0.) ? (nu12 * k1 + nu24 * k4) / k14 : 0.;
    double nu14_3 = (k14 != 0.) ? (nu13 * k1 + nu34 * k4) / k14 : 0.;

    double nu23_1 = (k14 != 0.) ? (nu12 * k2 + nu13 * k3) / k14 : 0.;
    double nu23_4 = (k14 != 0.) ? (nu24 * k2 + nu34 * k3) / k14 : 0.;
    double nu24_1 = (k13 != 0.) ? (nu12 * k2 + nu14 * k4) / k13 : 0.;
    double nu24_3 = (k13 != 0.) ? (nu23 * k2 + nu34 * k4) / k13 : 0.;

    double nu34_1 = (k12 != 0.) ? (nu13 * k3 + nu14 * k4) / k12 : 0.;
    double nu34_2 = (k12 != 0.) ? (nu23 * k3 + nu24 * k4) / k12 : 0.;


    /**  Spectra  **/

    /* Bispectrum */

    /* B(k12_, k3_, k4_) */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_k(kern, 2, k4);

    kernels_qset_mu(kern, 0, mu12);
    kernels_qset_mu(kern, 1, mu3);
    kernels_qset_mu(kern, 2, mu4);

    kernels_qset_nu(kern, 0, 1, nu12_3);
    kernels_qset_nu(kern, 0, 2, nu12_4);
    kernels_qset_nu(kern, 1, 2, nu34);

    double btr12 = _specBtr_(specArg, NULL);

    /* B(k13_, k2_, k4_) */
    kernels_qset_k(kern, 0, k13);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k4);

    kernels_qset_mu(kern, 0, mu13);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu4);

    kernels_qset_nu(kern, 0, 1, nu13_2);
    kernels_qset_nu(kern, 0, 2, nu13_4);
    kernels_qset_nu(kern, 1, 2, nu24);

    double btr13 = _specBtr_(specArg, NULL);

    /* B(k14_, k2_, k3_) */
    kernels_qset_k(kern, 0, k14);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_k(kern, 2, k3);

    kernels_qset_mu(kern, 0, mu14);
    kernels_qset_mu(kern, 1, mu2);
    kernels_qset_mu(kern, 2, mu3);

    kernels_qset_nu(kern, 0, 1, nu14_2);
    kernels_qset_nu(kern, 0, 2, nu14_3);
    kernels_qset_nu(kern, 1, 2, nu23);

    double btr14 = _specBtr_(specArg, NULL);

    /* B(k23_, k1_, k4_) = B(-k14_, k1_, k4_) */
    kernels_qset_k(kern, 0, k14);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_k(kern, 2, k4);

    kernels_qset_mu(kern, 0, -mu14);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_mu(kern, 2, mu4);

    kernels_qset_nu(kern, 0, 1, nu23_1);
    kernels_qset_nu(kern, 0, 2, nu23_4);
    kernels_qset_nu(kern, 1, 2, nu14);

    double btr23 = _specBtr_(specArg, NULL);

    /* B(k24_, k1_, k3_) = B(-k13_, k1_, k3_) */
    kernels_qset_k(kern, 0, k13);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_k(kern, 2, k3);

    kernels_qset_mu(kern, 0, -mu13);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_mu(kern, 2, mu3);

    kernels_qset_nu(kern, 0, 1, nu24_1);
    kernels_qset_nu(kern, 0, 2, nu24_3);
    kernels_qset_nu(kern, 1, 2, nu13);

    double btr24 = _specBtr_(specArg, NULL);

    /* B(k34_, k1_, k2_) = B(-k12_, k1_, k2_) */
    kernels_qset_k(kern, 0, k12);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_k(kern, 2, k2);

    kernels_qset_mu(kern, 0, -mu12);
    kernels_qset_mu(kern, 1, mu1);
    kernels_qset_mu(kern, 2, mu2);

    kernels_qset_nu(kern, 0, 1, nu34_1);
    kernels_qset_nu(kern, 0, 2, nu34_2);
    kernels_qset_nu(kern, 1, 2, nu12);

    double btr34 = _specBtr_(specArg, NULL);


    /* Power spectrum */

    /* P(k1_) */
    kernels_qset_k(kern, 0, k1);
    kernels_qset_k(kern, 1, k1);
    kernels_qset_mu(kern, 0, mu1);
    kernels_qset_mu(kern, 1, -mu1);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl1 = _specPnl_(specArg, NULL);

    /* P(k2_) */
    kernels_qset_k(kern, 0, k2);
    kernels_qset_k(kern, 1, k2);
    kernels_qset_mu(kern, 0, mu2);
    kernels_qset_mu(kern, 1, -mu2);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl2 = _specPnl_(specArg, NULL);

    /* P(k3_) */
    kernels_qset_k(kern, 0, k3);
    kernels_qset_k(kern, 1, k3);
    kernels_qset_mu(kern, 0, mu3);
    kernels_qset_mu(kern, 1, -mu3);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl3 = _specPnl_(specArg, NULL);

    /* P(k4_) */
    kernels_qset_k(kern, 0, k4);
    kernels_qset_k(kern, 1, k4);
    kernels_qset_mu(kern, 0, mu4);
    kernels_qset_mu(kern, 1, -mu4);
    kernels_qset_nu(kern, 0, 1, -1.);

    double pnl4 = _specPnl_(specArg, NULL);


    /* Result */

    result[0] = kern -> surv -> sn * (btr12 + btr13 + btr14 + btr23 + btr24 + btr34)
                - 2. * kern -> surv -> sn*kern -> surv -> sn * (pnl1 + pnl2 + pnl3 + pnl4)
                + 3. * kern -> surv -> sn*kern -> surv -> sn*kern -> surv -> sn;

    result[1] = 0.;

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

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int spec_ttr_full(spec_arg_t *specArg, double *result)
{
    /*

        Calculate the full tree-level trispectrum (ie tree-level + shot noise)

    */

    /* Must have specOrder = 4 */
    kern_t *kern = specArg -> kern;
    size_t specOrder = kern -> specOrder;
    kern -> specOrder = 4;

    /* Tree level */
    double resultTtrTree[2];
    spec_ttr_tree(specArg, resultTtrTree);

    /* Shot noise */
    double resultTtrSn[2];
    spec_ttr_sn(specArg, resultTtrSn);

    /* Result */
    result[0] = resultTtrTree[0] + resultTtrSn[0];
    result[1] = sqrt(resultTtrTree[1]*resultTtrTree[1] + resultTtrSn[1]*resultTtrSn[1]);

    /* Reset specOrder */
    kern -> specOrder = specOrder;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*spec_ttr(const char *label))(spec_arg_t*, double*)
{
    /*

        Return the requested function that calculates part (or full) of the tree level trispectrum

    */

    /* Tree-Level Trispectrum (excluding shot noise) */
    if (!strcmp(label, _ttrLabelTree))
      {
        return spec_ttr_tree;
      }

    /* Trispectrum Shot Noise */
    if (!strcmp(label, _ttrLabelSn))
      {
        return spec_ttr_sn;
      }

    /* Full Trispectrum */
    if (!strcmp(label, _ttrLabelFull))
      {
        return spec_ttr_full;
      }

    /* Label not found */
    return _spec_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------   (Indirect) Tree-Level Trispectrum   --------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _spec_ttr(void *var, void *params)
{
    /*

        Calculate the tree-level trispectrum (plus shot noise)

        Note: This indirect function is used as the initial function for 'double _specTtr_(void*, void*)'.
              The direct function 'int spec_ttr(spec_arg_t*, double*)' is of different shape (but here they calculate the same thing in the same way)!

    */

    (void) params;

    double result[3];
    spec_ttr_full(var, result);

    return result[0];
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     POLY SPECTRUM     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   (Direct) Poly Spectrum   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*spec_poly(const char *specLabel, const char *label))(spec_arg_t*, double*)
{
    /*

        Get the 'label' part of the poly spectrum with label 'specLabel'.

    */

    /**  Power Spectrum  **/

    /* Non-Linear Power Spectrum */
    if (!strcmp(specLabel, _idSpecPnl_))
      {
        return spec_pnl(label);
      }

    /* Non-Linear Power Spectrum Derivatives */
    if (!strcmp(specLabel, _idSpecDPnl_))
      {
        return spec_dpnl(label);
      }


    /**  Bispectrum  **/

    /* Tree-Level Bispectrum */
    if (!strcmp(specLabel, _idSpecBtr_))
      {
        return spec_btr(label);
      }

    /* Tree-Level Bispectrum Derivatives */
    if (!strcmp(specLabel, _idSpecDBtr_))
      {
        return spec_dbtr(label);
      }


    /**  Trispectrum  **/

    /* Tree-Level Trispectrum */
    if (!strcmp(specLabel, _idSpecTtr_))
      {
        return spec_ttr(label);
      }


    /**  Label not found  **/

    return _spec_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   (Indirect) Poly Spectrum   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double (*_spec_poly(const char *specLabel, const char *label))(void*, void*)
{
    /*

        Get the 'label' part of the poly spectrum with label 'specLabel'.

    */

    /**  Power Spectrum  **/

    /* Non-Linear Power Spectrum */
    if (!strcmp(specLabel, _idSpecPnl_))
      {
        return _specPnl_;
      }

    /* Non-Linear Power Spectrum Derivatives */
    if (!strcmp(specLabel, _idSpecDPnl_))
      {
        return _specDPnl_(label);
      }


    /**  Bispectrum  **/

    /* Tree-Level Bispectrum */
    if (!strcmp(specLabel, _idSpecBtr_))
      {
        return _specBtr_;
      }

    /* Tree-Level Bispectrum Derivatives */
    if (!strcmp(specLabel, _idSpecDBtr_))
      {
        return _specDBtr_(label);
      }


    /**  Trispectrum  **/

    /* Tree-Level Trispectrum */
    if (!strcmp(specLabel, _idSpecTtr_))
      {
        return _specTtr_;
      }


    /**  Label not found  **/

    return _zeroFunc_;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

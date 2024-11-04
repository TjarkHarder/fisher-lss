/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     FIDUCIALS.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "fiducials.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE / SETUP / FREE VARIABLES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Local Variables   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Fiducial file names  **/

static const char *_fidDir = "/fiducials/";

static char *_fidFile;
static char *_fidFilePk;


/**  Interpolation structures for the fiducials  **/

typedef struct
{
    /*

        Interpolate a fiducial parameter from the file

    */

    interp_t *interp;

    interp_t *(*interpInit)(dat_t*, double (*)(void*, void*, void*));
    double (*interpEval)(void*, interp_t*, void*);
    interp_t *(*interpFree)(interp_t*);

} _fid_interp_single_t;


typedef struct
{
    /*

        Interpolate fiducial parameters from the file

    */

    size_t size;
    char **labels;
    _fid_interp_single_t **fidInterpSingle;

} _fid_interp_t;



/* Number of fiducials */
static const size_t _fidSize = 22;

/* Fiducials interpolation struct */
static _fid_interp_t *_fidInterp = NULL;


/* LCDM */

/* omegaM0 */
static char *_fidLabelOmegaM0 = "omegaM0";
static _fid_interp_single_t _fidInterpOmegaM0;

/* ns */
static char *_fidLabelNs = "ns";
static _fid_interp_single_t _fidInterpNs;

/* growthIndex */
static char *_fidLabelGrowthIndex = "growthIndex";
static _fid_interp_single_t _fidInterpGrowthIndex;


/* Bootstrap */

/* a2Ga */
static char *_fidLabelA2Ga = "a2Ga";
static _fid_interp_single_t _fidInterpA2Ga;

/* d2Ga */
static char *_fidLabelD2Ga = "d2Ga";
static _fid_interp_single_t _fidInterpD2Ga;

/* a3GaA */
static char *_fidLabelA3GaA = "a3GaA";
static _fid_interp_single_t _fidInterpA3GaA;

/* a3GaB */
static char *_fidLabelA3GaB = "a3GaB";
static _fid_interp_single_t _fidInterpA3GaB;

/* d3GaA */
static char *_fidLabelD3GaA = "d3GaA";
static _fid_interp_single_t _fidInterpD3GaA;

/* d3GaB */
static char *_fidLabelD3GaB = "d3GaB";
static _fid_interp_single_t _fidInterpD3GaB;


/* Bias */

/* b1 */
static char *_fidLabelB1 = "b1";
static _fid_interp_single_t _fidInterpB1;

/* b2 */
static char *_fidLabelB2 = "b2";
static _fid_interp_single_t _fidInterpB2;

/* bG2 */
static char *_fidLabelBG2 = "bG2";
static _fid_interp_single_t _fidInterpBG2;

/* bGam3 */
static char *_fidLabelBGam3 = "bGam3";
static _fid_interp_single_t _fidInterpBGam3;


/* RSD */

/* sigv */
static char *_fidLabelSigv = "sigv";
static _fid_interp_single_t _fidInterpSigv;

/* sigs */
static char *_fidLabelSigs = "sigs";
static _fid_interp_single_t _fidInterpSigs;


/* Conter terms */

/* c0 */
static char *_fidLabelC0 = "c0";
static _fid_interp_single_t _fidInterpC0;

/* c2 */
static char *_fidLabelC2 = "c2";
static _fid_interp_single_t _fidInterpC2;

/* c4 */
static char *_fidLabelC4 = "c4";
static _fid_interp_single_t _fidInterpC4;


/* Survey */

/* V */
static char *_fidLabelV = "V";
static _fid_interp_single_t _fidInterpV;

/* n */
static char *_fidLabelN = "n";
static _fid_interp_single_t _fidInterpN;


/* Linear Power Spectrum */

/* Pk */
static char *_fidLabelPk = "Pk";
static _fid_interp_single_t _fidInterpPk;

/* dPk */
static char *_fidLabelDPk = "dPk";
static _fid_interp_single_t _fidInterpDPk;



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Initialise Local Variables   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_file(void);
static int _ini_interp(void);
static int _ini_func(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int fiducials_ini(void)
{
    /*

        Initialise local variables

    */

    /* Fiducial file names */
    _ini_file();

    /* Assign interpolation functions */
    _ini_interp();

    /* Assign evaluation functions */
    _ini_func();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_file(void)
{
    /*

        Default file name of the fiducials

    */

    _fidFile = NULL;
    _fidFilePk = NULL;

    misc_scp(&_fidFile, "fiducials.dat");
    misc_scp(&_fidFilePk, "pk.dat");

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_interp(void)
{
    /*

        Default interpolation functions (init, eval, free) for the fiducials

    */

    _fidInterp = malloc(sizeof(_fid_interp_t));

    _fidInterp -> size = _fidSize;
    _fidInterp -> labels = malloc(sizeof(char*) * _fidInterp -> size);
    _fidInterp -> fidInterpSingle = malloc(sizeof(_fid_interp_single_t*) * _fidInterp -> size);


    /* LCDM */

    /* omegaM0 */
    _fidInterp -> labels[0] = _fidLabelOmegaM0;
    _fidInterp -> fidInterpSingle[0] = &_fidInterpOmegaM0;

    /* ns */
    _fidInterp -> labels[1] = _fidLabelNs;
    _fidInterp -> fidInterpSingle[1] = &_fidInterpNs;

    /* growthIndex */
    _fidInterp -> labels[2] = _fidLabelGrowthIndex;
    _fidInterp -> fidInterpSingle[2] = &_fidInterpGrowthIndex;


    /* Bootstrap */

    /* a2Ga */
    _fidInterp -> labels[3] = _fidLabelA2Ga;
    _fidInterp -> fidInterpSingle[3] = &_fidInterpA2Ga;

    /* d2Ga */
    _fidInterp -> labels[4] = _fidLabelD2Ga;
    _fidInterp -> fidInterpSingle[4] = &_fidInterpD2Ga;

    /* a3GaA */
    _fidInterp -> labels[5] = _fidLabelA3GaA;
    _fidInterp -> fidInterpSingle[5] = &_fidInterpA3GaA;

    /* a3GaB */
    _fidInterp -> labels[6] = _fidLabelA3GaB;
    _fidInterp -> fidInterpSingle[6] = &_fidInterpA3GaB;

    /* d3GaA */
    _fidInterp -> labels[7] = _fidLabelD3GaA;
    _fidInterp -> fidInterpSingle[7] = &_fidInterpD3GaA;

    /* d3GaB */
    _fidInterp -> labels[8] = _fidLabelD3GaB;
    _fidInterp -> fidInterpSingle[8] = &_fidInterpD3GaB;


    /* Bias */

    /* b1 */
    _fidInterp -> labels[9] = _fidLabelB1;
    _fidInterp -> fidInterpSingle[9] = &_fidInterpB1;

    /* b2 */
    _fidInterp -> labels[10] = _fidLabelB2;
    _fidInterp -> fidInterpSingle[10] = &_fidInterpB2;

    /* bG2 */
    _fidInterp -> labels[11] = _fidLabelBG2;
    _fidInterp -> fidInterpSingle[11] = &_fidInterpBG2;

    /* bGam3 */
    _fidInterp -> labels[12] = _fidLabelBGam3;
    _fidInterp -> fidInterpSingle[12] = &_fidInterpBGam3;


    /* RSD */

    /* sigv */
    _fidInterp -> labels[13] = _fidLabelSigv;
    _fidInterp -> fidInterpSingle[13] = &_fidInterpSigv;

    /* sigs */
    _fidInterp -> labels[14] = _fidLabelSigs;
    _fidInterp -> fidInterpSingle[14] = &_fidInterpSigs;


    /* Ctr */

    /* c0 */
    _fidInterp -> labels[15] = _fidLabelC0;
    _fidInterp -> fidInterpSingle[15] = &_fidInterpC0;

    /* c2 */
    _fidInterp -> labels[16] = _fidLabelC2;
    _fidInterp -> fidInterpSingle[16] = &_fidInterpC2;

    /* c4 */
    _fidInterp -> labels[17] = _fidLabelC4;
    _fidInterp -> fidInterpSingle[17] = &_fidInterpC4;


    /* Survey */

    /* V */
    _fidInterp -> labels[18] = _fidLabelV;
    _fidInterp -> fidInterpSingle[18] = &_fidInterpV;

    /* n */
    _fidInterp -> labels[19] = _fidLabelN;
    _fidInterp -> fidInterpSingle[19] = &_fidInterpN;


    /* Linear Power Spectrum */
    _fidInterp -> labels[20] = _fidLabelPk;
    _fidInterp -> fidInterpSingle[20] = &_fidInterpPk;


    /* Linear Power Spectrum Derivative */
    _fidInterp -> labels[21] = _fidLabelDPk;
    _fidInterp -> fidInterpSingle[21] = &_fidInterpDPk;


    /* Set the interpolation functions */

    for (size_t i = 0; i < _fidInterp -> size; i++)
      {
        _fidInterp -> fidInterpSingle[i] -> interp = NULL;

        if (!strcmp(_fidInterp -> labels[i], "Pk"))
          {
            _fidInterp -> fidInterpSingle[i] -> interpInit = interpolate_interp_init_gsl;
            _fidInterp -> fidInterpSingle[i] -> interpEval = interpolate_interp_eval_gsl;
            _fidInterp -> fidInterpSingle[i] -> interpFree = interpolate_interp_free_gsl;
          }

        else if (!strcmp(_fidInterp -> labels[i], "dPk"))
          {
            _fidInterp -> fidInterpSingle[i] -> interpInit = interpolate_interp_init_gsl;
            _fidInterp -> fidInterpSingle[i] -> interpEval = interpolate_interp_eval_deriv_gsl;
            _fidInterp -> fidInterpSingle[i] -> interpFree = interpolate_interp_free_gsl;
          }

        else
          {
            _fidInterp -> fidInterpSingle[i] -> interpInit = interpolate_interp_init_splinter;
            _fidInterp -> fidInterpSingle[i] -> interpEval = interpolate_interp_eval_splinter;
            _fidInterp -> fidInterpSingle[i] -> interpFree = interpolate_interp_free_splinter;
          }
      }


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

/* LCDM */

static fid_lcdm_t *_fid_lcdm(void *var, void *params);

static double _fid_lcdm_omegam0(void *var, void *params);
static double _fid_lcdm_ns(void *var, void *params);
static double _fid_lcdm_growthindex(void *var, void *params);


/* Bootstrap */

static fid_btst_t *_fid_btst(void *var, void *params);

static double _fid_btst_a2ga(void *var, void *params);
static double _fid_btst_d2ga(void *var, void *params);

static double _fid_btst_a3gaa(void *var, void *params);
static double _fid_btst_a3gab(void *var, void *params);

static double _fid_btst_d3gaa(void *var, void *params);
static double _fid_btst_d3gab(void *var, void *params);

static double _fid_btst_h(void *var, void *params);


/* Bias */

static fid_bias_t *_fid_bias(void *var, void *params);

static double _fid_bias_b1(void *var, void *params);
static double _fid_bias_b2(void *var, void *params);

static double _fid_bias_bg2(void *var, void *params);
static double _fid_bias_c2ga(void *var, void *params);

static double _fid_bias_bgam3(void *var, void *params);


/* RSD */

static fid_rsd_t *_fid_rsd(void *var, void *params);

static double _fid_rsd_f(void *var, void *params);

static double _fid_rsd_sigv(void *var, void *params);
static double _fid_rsd_sigs(void *var, void *params);


/* Ctr */

static fid_ctr_t *_fid_ctr(void *var, void *params);

static double _fid_ctr_c0(void *var, void *params);
static double _fid_ctr_c2(void *var, void *params);
static double _fid_ctr_c4(void *var, void *params);


/* Survey */

static fid_surv_t *_fid_surv(void *var, void *params);

static double _fid_surv_v(void *var, void *params);
static double _fid_surv_n(void *var, void *params);
static double _fid_surv_sn(void *var, void *params);


/* Growth */

static double _fid_growth(void *var, void *params);


/* Linear Power Spectrum */

static double _fid_pk(void *var, void *params);
static double _fid_dpk(void *var, void *params);

static double _fid_pk_extrapolate(void *var, void *interp, void *params);
static double _fid_dpk_extrapolate(void *var, void *interp, void *params);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func(void)
{
    /*

        Initialise the fiducial functions that will be called in spectra and fisher modules

    */


    /* LCDM */

    _fidLCDM_ = _fid_lcdm;
    _fidParamsLCDM_ = NULL;


    _fidLCDMxOmegaM0_ = _fid_lcdm_omegam0;
    _fidParamsLCDMxOmegaM0_ = NULL;

    _fidLCDMxNs_ = _fid_lcdm_ns;
    _fidParamsLCDMxNs_ = NULL;

    _fidLCDMxGrowthIndex_ = _fid_lcdm_growthindex;
    _fidParamsLCDMxGrowthIndex_ = NULL;


    /* Bootstrap */

    _fidBTST_ = _fid_btst;
    _fidParamsBTST_ = NULL;


    _fidBTSTxA2Ga_ = _fid_btst_a2ga;
    _fidParamsBTSTxA2Ga_ = NULL;

    _fidBTSTxD2Ga_ = _fid_btst_d2ga;
    _fidParamsBTSTxD2Ga_ = NULL;


    _fidBTSTxA3GaA_ = _fid_btst_a3gaa;
    _fidParamsBTSTxA3GaA_ = NULL;

    _fidBTSTxA3GaB_ = _fid_btst_a3gab;
    _fidParamsBTSTxA3GaB_ = NULL;

    _fidBTSTxD3GaA_ = _fid_btst_d3gaa;
    _fidParamsBTSTxD3GaA_ = NULL;

    _fidBTSTxD3GaB_ = _fid_btst_d3gab;
    _fidParamsBTSTxD3GaB_ = NULL;


    _fidBTSTxH_ = _fid_btst_h;
    _fidParamsBTSTxH_ = NULL;


    /* Bias */

    _fidBias_ = _fid_bias;
    _fidParamsBias_ = NULL;


    _fidBiasxB1_ = _fid_bias_b1;
    _fidParamsBiasxB1_ = NULL;

    _fidBiasxB2_ = _fid_bias_b2;
    _fidParamsBiasxB2_ = NULL;


    _fidBiasxBG2_ = _fid_bias_bg2;
    _fidParamsBiasxBG2_ = NULL;

    _fidBiasxC2Ga_ = _fid_bias_c2ga;
    _fidParamsBiasxC2Ga_ = NULL;


    _fidBiasxBGam3_ = _fid_bias_bgam3;
    _fidParamsBiasxBGam3_ = NULL;


    /* RSD */

    _fidRSD_ = _fid_rsd;
    _fidParamsRSD_ = NULL;


    _fidRSDxF_ = _fid_rsd_f;
    _fidParamsRSDxF_ = NULL;


    _fidRSDxSigv_ = _fid_rsd_sigv;
    _fidParamsRSDxSigv_ = NULL;

    _fidRSDxSigs_ = _fid_rsd_sigs;
    _fidParamsRSDxSigs_ = NULL;


    /* Ctr */

    _fidCtr_ = _fid_ctr;
    _fidParamsCtr_ = NULL;


    _fidCtrxC0_ = _fid_ctr_c0;
    _fidParamsCtrxC0_ = NULL;

    _fidCtrxC2_ = _fid_ctr_c2;
    _fidParamsCtrxC2_ = NULL;

    _fidCtrxC4_ = _fid_ctr_c4;
    _fidParamsCtrxC4_ = NULL;


    /* Survey */

    _fidSurv_ = _fid_surv;
    _fidParamsSurv_ = NULL;


    _fidSurvxV_ = _fid_surv_v;
    _fidParamsSurvxV_ = NULL;

    _fidSurvxN_ = _fid_surv_n;
    _fidParamsSurvxN_ = NULL;

    _fidSurvxSn_ = _fid_surv_sn;
    _fidParamsSurvxSn_ = NULL;



    /* Growth */

    _fidGrowth_ = _fid_growth;
    _fidParamsGrowth_ = NULL;


    /* Pk */

    _fidPk_ = _fid_pk;
    _fidExtrapPk_ = _fid_pk_extrapolate;
    _fidParamsPk_ = NULL;

    _fidDPk_ = _fid_dpk;
    _fidExtrapDPk_ = _fid_dpk_extrapolate;
    _fidParamsDPk_ = NULL;


    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Setup Local Variables   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _setup_interp(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int fiducials_setup(void)
{
    /*

        Setup the fiducials by reading the files and interpolating the fiducials

    */


    /* Interpolate the fiducials */
    _setup_interp();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _setup_interp(void)
{
    /*

        Interpolate the fiducials after reading the in from file

    */


    dat_t *data;

    char *file = NULL;
    char *label;

    double (*extrapolate)(void*, void*, void*);

    for (size_t i = 0; i < _fidInterp -> size; i++)
      {
        /* Set the file name, label of the fiducial and the extrapolation function */
        if (!strcmp(_fidInterp -> labels[i], "Pk"))
          {
            file = misc_scat(3, __FISHERLSSDIR__, _fidDir, _fidFilePk);

            label = NULL;

            extrapolate = _fidExtrapPk_;
          }

        else if (!strcmp(_fidInterp -> labels[i], "dPk"))
          {
            file = misc_scat(3, __FISHERLSSDIR__, _fidDir, _fidFilePk);

            label = NULL;

            extrapolate = _fidExtrapDPk_;
          }

        else
          {
            file = misc_scat(3, __FISHERLSSDIR__, _fidDir, _fidFile);

            label = _fidInterp -> labels[i];

            extrapolate = NULL;
          }

        /* Get the fiducial and its interpolation function */
        data = dat_input(file, label, true);

        _fid_interp_single_t *fidInterp = _fidInterp -> fidInterpSingle[i];

        fidInterp -> interp = interpolate_interp_free(fidInterp -> interpFree, fidInterp -> interp);
        fidInterp -> interp = interpolate_interp_init(fidInterp -> interpInit, data, extrapolate);

        data = dat_free(data);
        free(file);
      }


    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Set Local Variables   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int fiducials_set_file(char type, const char *file)
{
    /*

        Set the file name for the fiducials

    */

    /* Power spectrum fiducial file */
    if (type == 'P')
      {
        misc_scp(&_fidFilePk, file);
      }

    /* Redshift dependent fiducials */
    else if (type == 'F')
      {
        misc_scp(&_fidFile, file);
      }

    else
      {
        printf("\nCan only set the fiducilal file name for 'F' or 'P' type.\n\n");

        exit(1);

        return 1;
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Free Local Variables   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_file(void);
static int _free_interp(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int fiducials_free(void)
{
    /*

        Free local fiducial variables

    */

    /* Free file name */
    _free_file();

    /* Free interpolation structs */
    _free_interp();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_file(void)
{
    /*

        Free the file names for the fiducials

    */

    free(_fidFile);
    free(_fidFilePk);

    _fidFile = NULL;
    _fidFilePk = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _free_interp(void)
{
    /*

        Free the fiducial interpolation structs

    */

    if (_fidInterp == NULL)
        return 0;

    for (size_t i = 0; i < _fidInterp -> size; i++)
      {
        _fidInterp -> fidInterpSingle[i] -> interp = interpolate_interp_free(_fidInterp -> fidInterpSingle[i] -> interpFree, _fidInterp -> fidInterpSingle[i] -> interp);
      }

    free(_fidInterp -> labels);
    free(_fidInterp -> fidInterpSingle);
    free(_fidInterp);

    _fidInterp = NULL;

    return 0;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     NEW / FREE / CP STRUCT FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------------   Lambda CDM   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_lcdm_t *fid_lcdm_new(void)
{
    /*

        Create a new fid_lcdm_t struct

    */

    /* Allocate memory */
    fid_lcdm_t *lcdm = malloc(sizeof(fid_lcdm_t));

    lcdm -> omegaM0 = 0.;
    lcdm -> ns = 0.;
    lcdm -> growthIndex = 0.;

    return lcdm;
}


fid_lcdm_t *fid_lcdm_free(fid_lcdm_t *lcdm)
{
    /*

        Free lcdm

    */

    /* Free lcdm */
    free(lcdm);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_lcdm_t *fid_lcdm_cp(fid_lcdm_t *lcdm)
{
    /*

        Copy a fid_lcdm_t struct

    */

    fid_lcdm_t *lcdmCp = fid_lcdm_new();

    lcdmCp -> omegaM0 = lcdm -> omegaM0;
    lcdmCp -> ns = lcdm -> ns;
    lcdmCp -> growthIndex = lcdm -> growthIndex;

    return lcdmCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------------   Beyond EdS   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_btst_t *fid_btst_new(void)
{
    /*

        Create a new fid_btst_t struct

    */

    /* Allocate memory */
    fid_btst_t *btst = malloc(sizeof(fid_btst_t));

    btst -> a2Ga = 0.;
    btst -> d2Ga = 0.;

    btst -> a3GaA = 0.;
    btst -> a3GaB = 0.;
    btst -> d3GaA = 0.;
    btst -> d3GaB = 0.;

    btst -> h = 0.;

    return btst;
}


fid_btst_t *fid_btst_free(fid_btst_t *btst)
{
    /*

        Free btst

    */

    /* Free btst */
    free(btst);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_btst_t *fid_btst_cp(fid_btst_t *btst)
{
    /*

        Copy a fid_btst_t struct

    */

    fid_btst_t *btstCp = fid_btst_new();

    btstCp -> a2Ga = btst -> a2Ga;
    btstCp -> d2Ga = btst -> d2Ga;

    btstCp -> a3GaA = btst -> a3GaA;
    btstCp -> a3GaB = btst -> a3GaB;

    btstCp -> d3GaA = btst -> d3GaA;
    btstCp -> d3GaB = btst -> d3GaB;

    btstCp -> h = btst -> h;

    return btstCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Bias   ----------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_bias_t *fid_bias_new(void)
{
    /*

        Create a new fid_bias_t struct

    */

    /* Allocate memory */
    fid_bias_t *bias = malloc(sizeof(fid_bias_t));

    /* Initialise the bias parameters */
    bias -> b1 = 1.;
    bias -> b2 = 0.;

    bias -> bG2 = 0.;
    bias -> c2Ga = 0.;

    bias -> bGam3 = 0.;

    return bias;
}


fid_bias_t *fid_bias_free(fid_bias_t *bias)
{
    /*

        Free bias

    */

    free(bias);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_bias_t *fid_bias_cp(fid_bias_t *bias)
{
    /*

        Copy a fid_bias_t struct

    */

    fid_bias_t *biasCp = fid_bias_new();

    biasCp -> b1 = bias -> b1;
    biasCp -> b2 = bias -> b2;

    biasCp -> c2Ga = bias -> c2Ga;
    biasCp -> bG2 = bias -> bG2;

    biasCp -> bGam3 = bias -> bGam3;

    return biasCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Redshift Space Distortions   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_rsd_t *fid_rsd_new(void)
{
    /*

        Create a new fid_rsd_t struct

    */

    /* Allocate memory */
    fid_rsd_t *rsd = malloc(sizeof(fid_rsd_t));

    return rsd;
}


fid_rsd_t *fid_rsd_free(fid_rsd_t *rsd)
{
    /*

        Free rsd

    */

    /* Free rsd */
    free(rsd);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_rsd_t *fid_rsd_cp(fid_rsd_t *rsd)
{
    /*

        Copy a fid_rsd_t struct

    */

    fid_rsd_t *rsdCp = fid_rsd_new();

    rsdCp -> f = rsd -> f;

    rsdCp -> sigv = rsd -> sigv;
    rsdCp -> sigs = rsd -> sigs;

    return rsdCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Counter Terms   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_ctr_t *fid_ctr_new(void)
{
    /*

        Create a new fid_ctr_t struct

    */

    /* Allocate memory */
    fid_ctr_t *ctr = malloc(sizeof(fid_ctr_t));

    return ctr;
}


fid_ctr_t *fid_ctr_free(fid_ctr_t *ctr)
{
    /*

        Free ctr

    */

    free(ctr);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_ctr_t *fid_ctr_cp(fid_ctr_t *ctr)
{
    /*

        Copy a fid_ctr_t struct

    */

    fid_ctr_t *ctrCp = fid_ctr_new();

    ctrCp -> c0 = ctr -> c0;
    ctrCp -> c2 = ctr -> c2;
    ctrCp -> c4 = ctr -> c4;

    return ctrCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------------   Survey   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_surv_t *fid_surv_new(void)
{
    /*

        Create a new fid_surv_t struct

    */

    /* Allocate memory */
    fid_surv_t *surv = malloc(sizeof(fid_surv_t));

    return surv;
}


fid_surv_t *fid_surv_free(fid_surv_t *surv)
{
    /*

        Free surv

    */

    free(surv);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fid_surv_t *fid_surv_cp(fid_surv_t *surv)
{
    /*

        Copy a fid_surv_t struct

    */

    fid_surv_t *survCp = fid_surv_new();

    survCp -> V = surv -> V;
    survCp -> n = surv -> n;
    survCp -> sn = surv -> sn;

    return survCp;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     FIDUCIAL INTERPOLATION FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   LCDM Fiducials   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static fid_lcdm_t *_fid_lcdm(void *var, void *params)
{
    /*

        Calculate the LCDM fiducials.

    */

    (void) params;

    fid_lcdm_t *lcdm = fid_lcdm_new();

    lcdm -> omegaM0 = _fidLCDMxOmegaM0_(var, _fidParamsLCDMxOmegaM0_);
    lcdm -> ns = _fidLCDMxNs_(var, _fidParamsLCDMxNs_);
    lcdm -> growthIndex = _fidLCDMxGrowthIndex_(var, _fidParamsLCDMxGrowthIndex_);

    return lcdm;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_lcdm_omegam0(void *var, void *params)
{
    /*

        Calculate the OmegaM0 fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpOmegaM0.interpEval, var, _fidInterpOmegaM0.interp, params);
}


static double _fid_lcdm_ns(void *var, void *params)
{
    /*

        Calculate the ns fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpNs.interpEval, var, _fidInterpNs.interp, params);
}


static double _fid_lcdm_growthindex(void *var, void *params)
{
    /*

        Calculate the GrowthIndex fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpGrowthIndex.interpEval, var, _fidInterpGrowthIndex.interp, params);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Bootstrap Fiducials   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static fid_btst_t *_fid_btst(void *var, void *params)
{
    /*

        Calculate the bootstrap fiducials.

    */

    (void) params;

    fid_btst_t *btst = fid_btst_new();

    btst -> a2Ga = _fidBTSTxA2Ga_(var, _fidParamsBTSTxA2Ga_);
    btst -> d2Ga = _fidBTSTxD2Ga_(var, _fidParamsBTSTxD2Ga_);

    btst -> a3GaA = _fidBTSTxA3GaA_(var, _fidParamsBTSTxA3GaA_);
    btst -> a3GaB = _fidBTSTxA3GaB_(var, _fidParamsBTSTxA3GaB_);

    btst -> d3GaA = _fidBTSTxD3GaA_(var, _fidParamsBTSTxD3GaA_);
    btst -> d3GaB = _fidBTSTxD3GaB_(var, _fidParamsBTSTxD3GaB_);

    btst -> h = _fidBTSTxH_(var, _fidParamsBTSTxH_);

    return btst;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_btst_a2ga(void *var, void *params)
{
    /*

        Calculate the a2Ga fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpA2Ga.interpEval, var, _fidInterpA2Ga.interp, params);
}


static double _fid_btst_d2ga(void *var, void *params)
{
    /*

        Calculate the d2Ga fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpD2Ga.interpEval, var, _fidInterpD2Ga.interp, params);
}


static double _fid_btst_a3gaa(void *var, void *params)
{
    /*

        Calculate the a3GaA fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpA3GaA.interpEval, var, _fidInterpA3GaA.interp, params);
}


static double _fid_btst_a3gab(void *var, void *params)
{
    /*

        Calculate the a3GaB fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpA3GaB.interpEval, var, _fidInterpA3GaB.interp, params);
}


static double _fid_btst_d3gaa(void *var, void *params)
{
    /*

        Calculate the d3GaA fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpD3GaA.interpEval, var, _fidInterpD3GaA.interp, params);
}


static double _fid_btst_d3gab(void *var, void *params)
{
    /*

        Calculate the d3GaB fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpD3GaB.interpEval, var, _fidInterpD3GaB.interp, params);
}


static double _fid_btst_h(void *var, void *params)
{
    /*

        Calculate the h fiducial via h = a2Ga - 1 (following Amendola et al. 2023 instead of D'Amico et al. 2021)

    */

    (void) params;

    return _fidBTSTxA2Ga_(var, _fidParamsBTSTxA2Ga_) - 1.;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Bias Fiducials   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static fid_bias_t *_fid_bias(void *var, void *params)
{
    /*

        Calculate the bias fiducials.

    */

    (void) params;

    fid_bias_t *bias = fid_bias_new();

    bias -> b1 = _fidBiasxB1_(var, _fidParamsBiasxB1_);
    bias -> b2 = _fidBiasxB2_(var, _fidParamsBiasxB2_);

    bias -> bG2 = _fidBiasxBG2_(var, _fidParamsBiasxBG2_);
    bias -> c2Ga = _fidBiasxC2Ga_(var, _fidParamsBiasxC2Ga_);

    bias -> bGam3 = _fidBiasxBGam3_(var, _fidParamsBiasxBGam3_);

    return bias;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_bias_b1(void *var, void *params)
{
    /*

        Calculate the b1 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclBias_)
      {
        return 1.;
      }

    return interpolate_interp_eval(_fidInterpB1.interpEval, var, _fidInterpB1.interp, params);
}


static double _fid_bias_b2(void *var, void *params)
{
    /*

        Calculate the b2 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclBias_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpB2.interpEval, var, _fidInterpB2.interp, params);
}


static double _fid_bias_bg2(void *var, void *params)
{
    /*

        Calculate the bG2 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclBias_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpBG2.interpEval, var, _fidInterpBG2.interp, params);
}


static double _fid_bias_c2ga(void *var, void *params)
{
    /*

        Calculate the c2Ga fiducial via c2Ga = b1 a2Ga - 2 bG2.

    */

    (void) params;

    return _fidBiasxB1_(var, _fidParamsBiasxB1_) * _fidBTSTxA2Ga_(var, _fidParamsBTSTxA2Ga_) - 2. * _fidBiasxBG2_(var, _fidParamsBiasxBG2_);
}


static double _fid_bias_bgam3(void *var, void *params)
{
    /*

        Calculate the bGam3 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclBias_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpBGam3.interpEval, var, _fidInterpBGam3.interp, params);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   RSD Fiducials   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static fid_rsd_t *_fid_rsd(void *var, void *params)
{
    /*

        Calculate the RSD fiducials.

    */

    (void) params;

    fid_rsd_t *rsd = fid_rsd_new();

    rsd -> f = _fidRSDxF_(var, _fidParamsRSDxF_);

    rsd -> sigv = _fidRSDxSigv_(var, _fidParamsRSDxSigv_);
    rsd -> sigs = _fidRSDxSigs_(var, _fidParamsRSDxSigs_);

    return rsd;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_rsd_f(void *var, void *params)
{
    /*

        Calculate the f fiducial via fid_growth_rate(z, lcdm).

    */

    (void) params;

    if (!_fidInclRSD_)
      {
        return 0.;
      }

    /* Need LCDM */
    fid_lcdm_t *lcdm = _fidLCDM_(var, _fidParamsLCDM_);

    double f = fid_growth_rate(*((double*) var), lcdm);

    /* Free memory */
    lcdm = fid_lcdm_free(lcdm);

    return f;
}


static double _fid_rsd_sigv(void *var, void *params)
{
    /*

        Calculate the sigv fiducial via sigv / sqrt(2) via the corresponding interpolation function.

    */

    if (!_fidInclRSD_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpSigv.interpEval, var, _fidInterpSigv.interp, params) / sqrt(2.);
}


static double _fid_rsd_sigs(void *var, void *params)
{
    /*

        Calculate the sigv fiducial via sigs (1 + z) / fid_hubble(z, lcdm) via the corresponding interpolation function.

    */

    if (!_fidInclRSD_)
      {
        return 0.;
      }

    /* Need LCDM */
    fid_lcdm_t *lcdm = _fidLCDM_(var, _fidParamsLCDM_);

    double sigs = interpolate_interp_eval(_fidInterpSigs.interpEval, var, _fidInterpSigs.interp, params) * (1 + *((double*) var)) / fid_hubble(*((double*) var), lcdm);

    /* Free memory */
    lcdm = fid_lcdm_free(lcdm);

    return sigs;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Ctr Fiducials   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static fid_ctr_t *_fid_ctr(void *var, void *params)
{
    /*

        Calculate the Ctr fiducials.

    */

    (void) params;

    fid_ctr_t *ctr = fid_ctr_new();

    ctr -> c0 = _fidCtrxC0_(var, _fidParamsCtrxC0_);
    ctr -> c2 = _fidCtrxC2_(var, _fidParamsCtrxC2_);
    ctr -> c4 = _fidCtrxC4_(var, _fidParamsCtrxC4_);

    return ctr;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_ctr_c0(void *var, void *params)
{
    /*

        Calculate the c0 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclCtr_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpC0.interpEval, var, _fidInterpC0.interp, params);
}


static double _fid_ctr_c2(void *var, void *params)
{
    /*

        Calculate the c2 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclCtr_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpC2.interpEval, var, _fidInterpC2.interp, params);
}


static double _fid_ctr_c4(void *var, void *params)
{
    /*

        Calculate the c4 fiducial via the corresponding interpolation function.

    */

    if (!_fidInclCtr_)
      {
        return 0.;
      }

    return interpolate_interp_eval(_fidInterpC4.interpEval, var, _fidInterpC4.interp, params);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Survey Fiducials   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static fid_surv_t *_fid_surv(void *var, void *params)
{
    /*

        Calculate the survey fiducials.

    */

    (void) params;

    fid_surv_t *surv = fid_surv_new();

    surv -> V = _fidSurvxV_(var, _fidParamsSurvxV_);
    surv -> n = _fidSurvxN_(var, _fidParamsSurvxN_);
    surv -> sn = _fidSurvxSn_(var, _fidParamsSurvxSn_);

    return surv;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_surv_v(void *var, void *params)
{
    /*

        Calculate the V fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpV.interpEval, var, _fidInterpV.interp, params);
}


static double _fid_surv_n(void *var, void *params)
{
    /*

        Calculate the n fiducial via the corresponding interpolation function.

    */

    return interpolate_interp_eval(_fidInterpN.interpEval, var, _fidInterpN.interp, params);
}


static double _fid_surv_sn(void *var, void *params)
{
    /*

        Calculate the sn fiducial via sn = 1 / n

    */

    (void) params;

    if (!_fidInclSn_)
      {
        return 0.;
      }

    return 1. / _fidSurvxN_(var, _fidParamsSurvxN_);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Growth Fiducial   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_growth(void *var, void *params)
{
    /*

        Calculate the growth fiducial via exp(-int f(z) / (1 + z)).

    */

    (void) params;

    /* Need LCDM */
    fid_lcdm_t *lcdm = _fidLCDM_(var, _fidParamsLCDM_);

    double growth = fid_growth_fct(*((double*) var), lcdm);

    /* Free memory */
    lcdm = fid_lcdm_free(lcdm);

    return growth;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   Linear Power Spectrum Fiducial   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_pk(void *var, void *params)
{
    /*

        Calculate the linear power spectrum fiducial via the corresponding interpolation fiducial.

    */

    return interpolate_interp_eval(_fidInterpPk.interpEval, var, _fidInterpPk.interp, params);
}


static double _fid_dpk(void *var, void *params)
{
    /*

        Calculate the linear power spectrum fiducial via the corresponding interpolation fiducial.

    */

    return interpolate_interp_eval(_fidInterpDPk.interpEval, var, _fidInterpDPk.interp, params);
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_pk_extrapolate(void *var, void *interp, void *params)
{
    /*

        TODO: (logarithmic dependence for large k...)

        Extrapolate the linear power spectrum:

                    P(k) ~ k^ns for k < kmin

        and

                    P(k) ~ k^(ns - 4) for k > kmax

        with kmin and kmax the boundaries of the power spectrum data.

    */

    (void) params;

    interp_t *interpPk = (interp_t*) interp;

    double k = *((double*) var);
    double pk;

    /* k smaller than fidInterpolated region (Note: ns is a constant so it does not matter what is fed to _fidLCDMxNs_) */
    if (k < interpPk -> xBounds[0][0])
        pk = interpPk -> yBounds[0] * pow(k / interpPk -> xBounds[0][0], _fidLCDMxNs_(NULL, _fidParamsLCDMxNs_));

    /* k larger than fidInterpolated region */
    else
        pk = interpPk -> yBounds[1] * pow(k / interpPk -> xBounds[0][1], _fidLCDMxNs_(NULL, _fidParamsLCDMxNs_) - 4.0);

    return pk;
}


static double _fid_dpk_extrapolate(void *var, void *interp, void *params)
{
    /*

        Extrapolate the derivative of the linear power spectrum

    */

    (void) params;

    interp_t *interpDPk = (interp_t*) interp;

    double k = *((double*) var);
    double dpk;
    double xBound[1];

    /* k smaller than fidInterpolated region */
    if (k < interpDPk -> xBounds[0][0])
      {
        xBound[0] = interpDPk -> xBounds[0][0];
        dpk = _fidDPk_(xBound, _fidParamsDPk_) * pow(k / interpDPk -> xBounds[0][0], _fidLCDMxNs_(NULL, _fidParamsLCDMxNs_) - 1.);
      }

    /* k larger than fidInterpolated region */
    else
      {
        xBound[0] = interpDPk -> xBounds[0][1];
        dpk = _fidDPk_(xBound, _fidParamsDPk_) * pow(k / interpDPk -> xBounds[0][1], _fidLCDMxNs_(NULL, _fidParamsLCDMxNs_) - 5.);
      }

    return dpk;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COSMOLOGY     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



double fid_omegam(double z, fid_lcdm_t *lcdm)
{
    /*

        Redshift dependent matter density.

    */

    return lcdm -> omegaM0 * pow(1. + z, 3.) / (lcdm -> omegaM0 * pow(1. + z, 3.) + 1. - lcdm -> omegaM0);
}


double fid_hubble(double z, fid_lcdm_t *lcdm)
{
    /*

        Hubble function in units of h * km / s / Mpc / c -> h / Mpc:

            H = 10^5 / 2.997 * 10^8 * h / Mpc = 1 / 2997 * h / Mpc

    */

    return 1. / (__C_LIGHT__ * 1.0e-5) * pow(lcdm -> omegaM0 * pow(1. + z, 3.) + 1. - lcdm -> omegaM0, 0.5);
}


/*  ------------------------------------------------------------------------------------------------------  */


double fid_growth_rate(double z, fid_lcdm_t *lcdm)
{
    /*

        Linear growth rate in LCDM:

            f(z) = dlog D / dlog a = pow(Omega_M(z), growthIndex)

    */

    return pow(fid_omegam(z, lcdm), lcdm -> growthIndex);
}


/*  ------------------------------------------------------------------------------------------------------  */


static double _fid_growth_fct_integrand(double *var, size_t dim, void *params)
{
    /*

        Integrand for the linear growth function:

            f(z) / (1 + z)

        where f is the (linear) growth rate.

    */

    /* Not used */
    (void) (dim);

    /* Redshift is integration variable */
    double z = var[0];

    /* LCDM parameters */
    fid_lcdm_t *lcdm = (fid_lcdm_t*) params;

    return fid_growth_rate(z, lcdm) / (1 + z);
}


double fid_growth_fct(double z, fid_lcdm_t *lcdm)
{
    /*

        Linear growth function in LCDM:

            D1(z) = exp(- integral(f(z') / (1 + z'), z', 0, z)).

    */

    intgrt_vegas_t *vegas = integrate_vegas_new();

    integrate_vegas_set_wucalls(vegas, 100);
    integrate_vegas_set_mcalls(vegas, 1000);
    integrate_vegas_set_chisqerr(vegas, 0.5);
    integrate_vegas_set_mloop(vegas, 10);
    integrate_vegas_set_verbose(vegas, 0);

    intgrt_t *intgrt = integrate_new();

    double lowerBound = 0.;
    double upperBound = z;

    integrate_set_bounds(intgrt, 1, &upperBound, &lowerBound);
    integrate_set_params(intgrt, lcdm);
    integrate_set_vegas(intgrt, vegas);

    double growthResult[3];

    integrate_vegas(_fid_growth_fct_integrand, intgrt, growthResult);

    intgrt = integrate_free(intgrt);
    vegas = integrate_vegas_free(vegas);

    return exp(-growthResult[0]);
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

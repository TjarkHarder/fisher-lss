/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     FISHER.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "fisher.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE / SETUP / FREE VARIABLES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Local Variables   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Output directory and file extension  **/

static const char *_outDir = "/output/fish/";
static const char *_outExt = ".mat";


/**  Information about the fisher matrix  **/

typedef struct
{
    /*

        Information about the fisher matrix in condensed form

    */

    const char *id;

    size_t specSize;
    const char **specLabels;
    size_t *specOrders;

    size_t covSize;
    const char **covLabels;

    size_t partsSize;

    char **partsLabels;
    bool **partsExist;
    int **partsMult;
    bool **partsDerivLog;

    dat_t **partsDerivVals;

} _fish_info_t;


/*  ------------------------------------------------------------------------------------------------------  */


static _fish_info_t *_fishppInfo = NULL;
static _fish_info_t *_fishbbInfo = NULL;
static _fish_info_t *_fishpbInfo = NULL;


static _fish_info_t *_fish_info_get_struct(const char *id)
{
    /*

        Get the info struct with given id

    */

    /* PP Fisher Matrix */
    if (!strcmp(id, _idFishPP_))
        return _fishppInfo;

    /* BB Fisher Matrix */
    if (!strcmp(id, _idFishBB_))
        return _fishbbInfo;

    /* PB Fisher Matrix */
    if (!strcmp(id, _idFishPB_))
        return _fishpbInfo;

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


static _fish_info_t *_fish_info_free(_fish_info_t *info)
{
    /*

        Free a _fish_info_t struct

    */

    if (info == NULL)
        return 0;

    free(info -> specLabels);
    free(info -> specOrders);

    free(info -> covLabels);

    free(info -> partsLabels);
    free(info -> partsExist);
    free(info -> partsMult);
    free(info -> partsDerivLog);

    for (size_t i = 0; i < info -> partsSize; i++)
        info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);

    free(info -> partsDerivVals);

    free(info);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


size_t fish_info_get_kern_order(const char *id)
{
    /*

        Get the max kern order of all spectra

    */

    _fish_info_t *info = _fish_info_get_struct(id);

    size_t kernOrder = 0;

    for (size_t i = 0; i < info -> specSize; i++)
      {
        size_t specKernOrder = spec_info_get_kern_order(info -> specLabels[i]);

        if (kernOrder < specKernOrder)
            kernOrder = specKernOrder;
      }

    return kernOrder;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Initialise Local Variables   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_info(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int fish_ini(void)
{
    /*


        Initialise fihser matrix modules

    */

    /* Initialise info */
    _ini_info();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_info_pp(void);
static int _ini_info_bb(void);
static int _ini_info_pb(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info(void)
{
    /*

        Initialise the info structs

    */

    /* PP Fisher Matrix */
    _ini_info_pp();

    /* BB Fisher Matrix */
    _ini_info_bb();

    /* PB Fisher Matrix */
    _ini_info_pb();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info_pp(void)
{
    /*

        Initialise the pp fisher matrix info struct

    */


    _fishppInfo = malloc(sizeof(_fish_info_t));

    _fishppInfo -> id = _idFishPP_;


    /**  Spectra Components  **/

    _fishppInfo -> specSize = 1;

    _fishppInfo -> specLabels = malloc(sizeof(const char*) * _fishppInfo -> specSize);
    _fishppInfo -> specLabels[0] = _idSpecDPnl_;

    _fishppInfo -> specOrders = malloc(sizeof(size_t) * _fishppInfo -> specSize);
    _fishppInfo -> specOrders[0] = spec_info_get_order(_idSpecDPnl_);


    /**  Covariance Matrix Components  **/

    _fishppInfo -> covSize = 1;

    _fishppInfo -> covLabels = malloc(sizeof(const char*) * _fishppInfo -> covSize);
    _fishppInfo -> covLabels[0] = _idCovPP_;


    /**  Parts of the Fisher Matrix  **/

    _fishppInfo -> partsSize = __ID_PARAMS_SIZE__;

    _fishppInfo -> partsLabels = malloc(sizeof(char*) * _fishppInfo -> partsSize);
    _fishppInfo -> partsExist = malloc(sizeof(bool*) * _fishppInfo -> partsSize);
    _fishppInfo -> partsMult = malloc(sizeof(int*) * _fishppInfo -> partsSize);
    _fishppInfo -> partsDerivLog = malloc(sizeof(bool*) * _fishppInfo -> partsSize);
    _fishppInfo -> partsDerivVals = malloc(sizeof(dat_t*) * _fishppInfo -> partsSize);

    /* Bootstrap */
    _fishppInfo -> partsLabels[0] = (char*) _idParamsA2Ga_;
    _fishppInfo -> partsExist[0] = (bool*) &_true_;
    _fishppInfo -> partsMult[0] = (int*) _multParamsA2Ga_;
    _fishppInfo -> partsDerivLog[0] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[0] = NULL;

    _fishppInfo -> partsLabels[1] = (char*) _idParamsD2Ga_;
    _fishppInfo -> partsExist[1] = (bool*) &_true_;
    _fishppInfo -> partsMult[1] = (int*) _multParamsD2Ga_;
    _fishppInfo -> partsDerivLog[1] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[1] = NULL;

    _fishppInfo -> partsLabels[2] = (char*) _idParamsA3GaA_;
    _fishppInfo -> partsExist[2] = (bool*) &_true_;
    _fishppInfo -> partsMult[2] = (int*) _multParamsA3GaA_;
    _fishppInfo -> partsDerivLog[2] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[2] = NULL;

    _fishppInfo -> partsLabels[3] = (char*) _idParamsA3GaB_;
    _fishppInfo -> partsExist[3] = (bool*) &_true_;
    _fishppInfo -> partsMult[3] = (int*) _multParamsA3GaB_;
    _fishppInfo -> partsDerivLog[3] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[3] = NULL;

    _fishppInfo -> partsLabels[4] = (char*) _idParamsD3GaA_;
    _fishppInfo -> partsExist[4] = (bool*) &_true_;
    _fishppInfo -> partsMult[4] = (int*) _multParamsD3GaA_;
    _fishppInfo -> partsDerivLog[4] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[4] = NULL;

    _fishppInfo -> partsLabels[5] = (char*) _idParamsD3GaB_;
    _fishppInfo -> partsExist[5] = (bool*) &_true_;
    _fishppInfo -> partsMult[5] = (int*) _multParamsD3GaB_;
    _fishppInfo -> partsDerivLog[5] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[5] = NULL;

    /* Bias */
    _fishppInfo -> partsLabels[6] = (char*) _idParamsB1_;
    _fishppInfo -> partsExist[6] = (bool*) &_true_;
    _fishppInfo -> partsMult[6] = (int*) _multParamsB1_;
    _fishppInfo -> partsDerivLog[6] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[6] = NULL;

    _fishppInfo -> partsLabels[7] = (char*) _idParamsB2_;
    _fishppInfo -> partsExist[7] = (bool*) &_true_;
    _fishppInfo -> partsMult[7] = (int*) _multParamsB2_;
    _fishppInfo -> partsDerivLog[7] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[7] = NULL;

    _fishppInfo -> partsLabels[8] = (char*) _idParamsC2Ga_;
    _fishppInfo -> partsExist[8] = (bool*) &_true_;
    _fishppInfo -> partsMult[8] = (int*) _multParamsC2Ga_;
    _fishppInfo -> partsDerivLog[8] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[8] = NULL;

    _fishppInfo -> partsLabels[9] = (char*) _idParamsBGam3_;
    _fishppInfo -> partsExist[9] = (bool*) &_true_;
    _fishppInfo -> partsMult[9] = (int*) _multParamsBGam3_;
    _fishppInfo -> partsDerivLog[9] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[9] = NULL;

    /* RSD */
    _fishppInfo -> partsLabels[10] = (char*) _idParamsF_;
    _fishppInfo -> partsExist[10] = (bool*) &_true_;
    _fishppInfo -> partsMult[10] = (int*) _multParamsF_;
    _fishppInfo -> partsDerivLog[10] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[10] = NULL;

    _fishppInfo -> partsLabels[11] = (char*) _idParamsSigv_;
    _fishppInfo -> partsExist[11] = (bool*) &_true_;
    _fishppInfo -> partsMult[11] = (int*) _multParamsSigv_;
    _fishppInfo -> partsDerivLog[11] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[11] = NULL;

    /* AP */
    _fishppInfo -> partsLabels[12] = (char*) _idParamsD_;
    _fishppInfo -> partsExist[12] = (bool*) &_true_;
    _fishppInfo -> partsMult[12] = (int*) _multParamsD_;
    _fishppInfo -> partsDerivLog[12] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[12] = NULL;

    _fishppInfo -> partsLabels[13] = (char*) _idParamsH_;
    _fishppInfo -> partsExist[13] = (bool*) &_true_;
    _fishppInfo -> partsMult[13] = (int*) _multParamsH_;
    _fishppInfo -> partsDerivLog[13] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[13] = NULL;

    /* Ctr */
    _fishppInfo -> partsLabels[14] = (char*) _idParamsC0_;
    _fishppInfo -> partsExist[14] = (bool*) &_true_;
    _fishppInfo -> partsMult[14] = (int*) _multParamsC0_;
    _fishppInfo -> partsDerivLog[14] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[14] = NULL;

    _fishppInfo -> partsLabels[15] = (char*) _idParamsC2_;
    _fishppInfo -> partsExist[15] = (bool*) &_true_;
    _fishppInfo -> partsMult[15] = (int*) _multParamsC2_;
    _fishppInfo -> partsDerivLog[15] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[15] = NULL;

    _fishppInfo -> partsLabels[16] = (char*) _idParamsC4_;
    _fishppInfo -> partsExist[16] = (bool*) &_true_;
    _fishppInfo -> partsMult[16] = (int*) _multParamsC4_;
    _fishppInfo -> partsDerivLog[16] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[16] = NULL;

    /* Shot noise */
    _fishppInfo -> partsLabels[17] = (char*) _idParamsPsn_;
    _fishppInfo -> partsExist[17] = (bool*) &_true_;
    _fishppInfo -> partsMult[17] = (int*) _multParamsPsn_;
    _fishppInfo -> partsDerivLog[17] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[17] = NULL;

    _fishppInfo -> partsLabels[18] = (char*) _idParamsBsn1_;
    _fishppInfo -> partsExist[18] = (bool*) &_true_;
    _fishppInfo -> partsMult[18] = (int*) _multParamsBsn1_;
    _fishppInfo -> partsDerivLog[18] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[18] = NULL;

    _fishppInfo -> partsLabels[19] = (char*) _idParamsBsn2_;
    _fishppInfo -> partsExist[19] = (bool*) &_true_;
    _fishppInfo -> partsMult[19] = (int*) _multParamsBsn2_;
    _fishppInfo -> partsDerivLog[19] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[19] = NULL;

    /* Linear power spectrum */
    _fishppInfo -> partsLabels[20] = (char*) _idParamsPk_;
    _fishppInfo -> partsExist[20] = (bool*) &_true_;
    _fishppInfo -> partsMult[20] = (int*) _multParamsPk_;
    _fishppInfo -> partsDerivLog[20] = (bool*) &_false_;
    _fishppInfo -> partsDerivVals[20] = NULL;


    return 0;
}


static int _ini_info_bb(void)
{
    /*

        Initialise the bb fisher matrix info struct

    */


    _fishbbInfo = malloc(sizeof(_fish_info_t));

    _fishbbInfo -> id = _idFishBB_;

    /**  Spectra Components  **/

    _fishbbInfo -> specSize = 1;

    _fishbbInfo -> specLabels = malloc(sizeof(const char*) * _fishbbInfo -> specSize);
    _fishbbInfo -> specLabels[0] = _idSpecDBtr_;

    _fishbbInfo -> specOrders = malloc(sizeof(size_t) * _fishbbInfo -> specSize);
    _fishbbInfo -> specOrders[0] = spec_info_get_order(_idSpecDBtr_);


    /**  Covariance Matrix Components  **/

    _fishbbInfo -> covSize = 1;

    _fishbbInfo -> covLabels = malloc(sizeof(const char*) * _fishbbInfo -> covSize);
    _fishbbInfo -> covLabels[0] = _idCovBB_;


    /**  Parts of the Fisher Matrix  **/

    _fishbbInfo -> partsSize = __ID_PARAMS_SIZE__;

    _fishbbInfo -> partsLabels = malloc(sizeof(char*) * _fishbbInfo -> partsSize);
    _fishbbInfo -> partsExist = malloc(sizeof(bool*) * _fishbbInfo -> partsSize);
    _fishbbInfo -> partsMult = malloc(sizeof(int*) * _fishbbInfo -> partsSize);
    _fishbbInfo -> partsDerivLog = malloc(sizeof(bool*) * _fishbbInfo -> partsSize);
    _fishbbInfo -> partsDerivVals = malloc(sizeof(dat_t*) * _fishbbInfo -> partsSize);

    /* Bootstrap */
    _fishbbInfo -> partsLabels[0] = (char*) _idParamsA2Ga_;
    _fishbbInfo -> partsExist[0] = (bool*) &_true_;
    _fishbbInfo -> partsMult[0] = (int*) _multParamsA2Ga_;
    _fishbbInfo -> partsDerivLog[0] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[0] = NULL;

    _fishbbInfo -> partsLabels[1] = (char*) _idParamsD2Ga_;
    _fishbbInfo -> partsExist[1] = (bool*) &_true_;
    _fishbbInfo -> partsMult[1] = (int*) _multParamsD2Ga_;
    _fishbbInfo -> partsDerivLog[1] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[1] = NULL;

    _fishbbInfo -> partsLabels[2] = (char*) _idParamsA3GaA_;
    _fishbbInfo -> partsExist[2] = (bool*) &_true_;
    _fishbbInfo -> partsMult[2] = (int*) _multParamsA3GaA_;
    _fishbbInfo -> partsDerivLog[2] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[2] = NULL;

    _fishbbInfo -> partsLabels[3] = (char*) _idParamsA3GaB_;
    _fishbbInfo -> partsExist[3] = (bool*) &_true_;
    _fishbbInfo -> partsMult[3] = (int*) _multParamsA3GaB_;
    _fishbbInfo -> partsDerivLog[3] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[3] = NULL;

    _fishbbInfo -> partsLabels[4] = (char*) _idParamsD3GaA_;
    _fishbbInfo -> partsExist[4] = (bool*) &_true_;
    _fishbbInfo -> partsMult[4] = (int*) _multParamsD3GaA_;
    _fishbbInfo -> partsDerivLog[4] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[4] = NULL;

    _fishbbInfo -> partsLabels[5] = (char*) _idParamsD3GaB_;
    _fishbbInfo -> partsExist[5] = (bool*) &_true_;
    _fishbbInfo -> partsMult[5] = (int*) _multParamsD3GaB_;
    _fishbbInfo -> partsDerivLog[5] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[5] = NULL;

    /* Bias */
    _fishbbInfo -> partsLabels[6] = (char*) _idParamsB1_;
    _fishbbInfo -> partsExist[6] = (bool*) &_true_;
    _fishbbInfo -> partsMult[6] = (int*) _multParamsB1_;
    _fishbbInfo -> partsDerivLog[6] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[6] = NULL;

    _fishbbInfo -> partsLabels[7] = (char*) _idParamsB2_;
    _fishbbInfo -> partsExist[7] = (bool*) &_true_;
    _fishbbInfo -> partsMult[7] = (int*) _multParamsB2_;
    _fishbbInfo -> partsDerivLog[7] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[7] = NULL;

    _fishbbInfo -> partsLabels[8] = (char*) _idParamsC2Ga_;
    _fishbbInfo -> partsExist[8] = (bool*) &_true_;
    _fishbbInfo -> partsMult[8] = (int*) _multParamsC2Ga_;
    _fishbbInfo -> partsDerivLog[8] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[8] = NULL;

    _fishbbInfo -> partsLabels[9] = (char*) _idParamsBGam3_;
    _fishbbInfo -> partsExist[9] = (bool*) &_true_;
    _fishbbInfo -> partsMult[9] = (int*) _multParamsBGam3_;
    _fishbbInfo -> partsDerivLog[9] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[9] = NULL;

    /* RSD */
    _fishbbInfo -> partsLabels[10] = (char*) _idParamsF_;
    _fishbbInfo -> partsExist[10] = (bool*) &_true_;
    _fishbbInfo -> partsMult[10] = (int*) _multParamsF_;
    _fishbbInfo -> partsDerivLog[10] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[10] = NULL;

    _fishbbInfo -> partsLabels[11] = (char*) _idParamsSigv_;
    _fishbbInfo -> partsExist[11] = (bool*) &_true_;
    _fishbbInfo -> partsMult[11] = (int*) _multParamsSigv_;
    _fishbbInfo -> partsDerivLog[11] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[11] = NULL;

    /* AP */
    _fishbbInfo -> partsLabels[12] = (char*) _idParamsD_;
    _fishbbInfo -> partsExist[12] = (bool*) &_true_;
    _fishbbInfo -> partsMult[12] = (int*) _multParamsD_;
    _fishbbInfo -> partsDerivLog[12] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[12] = NULL;

    _fishbbInfo -> partsLabels[13] = (char*) _idParamsH_;
    _fishbbInfo -> partsExist[13] = (bool*) &_true_;
    _fishbbInfo -> partsMult[13] = (int*) _multParamsH_;
    _fishbbInfo -> partsDerivLog[13] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[13] = NULL;

    /* Ctr */
    _fishbbInfo -> partsLabels[14] = (char*) _idParamsC0_;
    _fishbbInfo -> partsExist[14] = (bool*) &_true_;
    _fishbbInfo -> partsMult[14] = (int*) _multParamsC0_;
    _fishbbInfo -> partsDerivLog[14] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[14] = NULL;

    _fishbbInfo -> partsLabels[15] = (char*) _idParamsC2_;
    _fishbbInfo -> partsExist[15] = (bool*) &_true_;
    _fishbbInfo -> partsMult[15] = (int*) _multParamsC2_;
    _fishbbInfo -> partsDerivLog[15] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[15] = NULL;

    _fishbbInfo -> partsLabels[16] = (char*) _idParamsC4_;
    _fishbbInfo -> partsExist[16] = (bool*) &_true_;
    _fishbbInfo -> partsMult[16] = (int*) _multParamsC4_;
    _fishbbInfo -> partsDerivLog[16] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[16] = NULL;

    /* Shot noise */
    _fishbbInfo -> partsLabels[17] = (char*) _idParamsPsn_;
    _fishbbInfo -> partsExist[17] = (bool*) &_true_;
    _fishbbInfo -> partsMult[17] = (int*) _multParamsPsn_;
    _fishbbInfo -> partsDerivLog[17] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[17] = NULL;

    _fishbbInfo -> partsLabels[18] = (char*) _idParamsBsn1_;
    _fishbbInfo -> partsExist[18] = (bool*) &_true_;
    _fishbbInfo -> partsMult[18] = (int*) _multParamsBsn1_;
    _fishbbInfo -> partsDerivLog[18] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[18] = NULL;

    _fishbbInfo -> partsLabels[19] = (char*) _idParamsBsn2_;
    _fishbbInfo -> partsExist[19] = (bool*) &_true_;
    _fishbbInfo -> partsMult[19] = (int*) _multParamsBsn2_;
    _fishbbInfo -> partsDerivLog[19] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[19] = NULL;

    /* Linear power spectrum */
    _fishbbInfo -> partsLabels[20] = (char*) _idParamsPk_;
    _fishbbInfo -> partsExist[20] = (bool*) &_true_;
    _fishbbInfo -> partsMult[20] = (int*) _multParamsPk_;
    _fishbbInfo -> partsDerivLog[20] = (bool*) &_false_;
    _fishbbInfo -> partsDerivVals[20] = NULL;


    return 0;
}


static int _ini_info_pb(void)
{
    /*

        Initialise the pb fisher matrix info struct

    */


    _fishpbInfo = malloc(sizeof(_fish_info_t));

    _fishpbInfo -> id = _idFishPB_;

    /**  Spectra Components  **/

    _fishpbInfo -> specSize = 2;

    _fishpbInfo -> specLabels = malloc(sizeof(const char*) * _fishpbInfo -> specSize);
    _fishpbInfo -> specLabels[0] = _idSpecDPnl_;
    _fishpbInfo -> specLabels[1] = _idSpecDBtr_;

    _fishpbInfo -> specOrders = malloc(sizeof(size_t) * _fishpbInfo -> specSize);
    _fishpbInfo -> specOrders[0] = spec_info_get_order(_idSpecDPnl_);
    _fishpbInfo -> specOrders[1] = spec_info_get_order(_idSpecDBtr_);


    /**  Covariance Matrix Components  **/

    _fishpbInfo -> covSize = 3;

    _fishpbInfo -> covLabels = malloc(sizeof(const char*) * _fishpbInfo -> covSize);
    _fishpbInfo -> covLabels[0] = _idCovPP_;
    _fishpbInfo -> covLabels[1] = _idCovPB_;
    _fishpbInfo -> covLabels[2] = _idCovBB_;


    /**  Parts of the Fisher Matrix  **/

    _fishpbInfo -> partsSize = __ID_PARAMS_SIZE__;

    _fishpbInfo -> partsLabels = malloc(sizeof(char*) * _fishpbInfo -> partsSize);
    _fishpbInfo -> partsExist = malloc(sizeof(bool*) * _fishpbInfo -> partsSize);
    _fishpbInfo -> partsMult = malloc(sizeof(int*) * _fishpbInfo -> partsSize);
    _fishpbInfo -> partsDerivLog = malloc(sizeof(bool*) * _fishpbInfo -> partsSize);
    _fishpbInfo -> partsDerivVals = malloc(sizeof(dat_t*) * _fishpbInfo -> partsSize);

    /* Bootstrap */
    _fishpbInfo -> partsLabels[0] = (char*) _idParamsA2Ga_;
    _fishpbInfo -> partsExist[0] = (bool*) &_true_;
    _fishpbInfo -> partsMult[0] = (int*) _multParamsA2Ga_;
    _fishpbInfo -> partsDerivLog[0] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[0] = NULL;

    _fishpbInfo -> partsLabels[1] = (char*) _idParamsD2Ga_;
    _fishpbInfo -> partsExist[1] = (bool*) &_true_;
    _fishpbInfo -> partsMult[1] = (int*) _multParamsD2Ga_;
    _fishpbInfo -> partsDerivLog[1] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[1] = NULL;

    _fishpbInfo -> partsLabels[2] = (char*) _idParamsA3GaA_;
    _fishpbInfo -> partsExist[2] = (bool*) &_true_;
    _fishpbInfo -> partsMult[2] = (int*) _multParamsA3GaA_;
    _fishpbInfo -> partsDerivLog[2] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[2] = NULL;

    _fishpbInfo -> partsLabels[3] = (char*) _idParamsA3GaB_;
    _fishpbInfo -> partsExist[3] = (bool*) &_true_;
    _fishpbInfo -> partsMult[3] = (int*) _multParamsA3GaB_;
    _fishpbInfo -> partsDerivLog[3] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[3] = NULL;

    _fishpbInfo -> partsLabels[4] = (char*) _idParamsD3GaA_;
    _fishpbInfo -> partsExist[4] = (bool*) &_true_;
    _fishpbInfo -> partsMult[4] = (int*) _multParamsD3GaA_;
    _fishpbInfo -> partsDerivLog[4] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[4] = NULL;

    _fishpbInfo -> partsLabels[5] = (char*) _idParamsD3GaB_;
    _fishpbInfo -> partsExist[5] = (bool*) &_true_;
    _fishpbInfo -> partsMult[5] = (int*) _multParamsD3GaB_;
    _fishpbInfo -> partsDerivLog[5] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[5] = NULL;

    /* Bias */
    _fishpbInfo -> partsLabels[6] = (char*) _idParamsB1_;
    _fishpbInfo -> partsExist[6] = (bool*) &_true_;
    _fishpbInfo -> partsMult[6] = (int*) _multParamsB1_;
    _fishpbInfo -> partsDerivLog[6] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[6] = NULL;

    _fishpbInfo -> partsLabels[7] = (char*) _idParamsB2_;
    _fishpbInfo -> partsExist[7] = (bool*) &_true_;
    _fishpbInfo -> partsMult[7] = (int*) _multParamsB2_;
    _fishpbInfo -> partsDerivLog[7] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[7] = NULL;

    _fishpbInfo -> partsLabels[8] = (char*) _idParamsC2Ga_;
    _fishpbInfo -> partsExist[8] = (bool*) &_true_;
    _fishpbInfo -> partsMult[8] = (int*) _multParamsC2Ga_;
    _fishpbInfo -> partsDerivLog[8] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[8] = NULL;

    _fishpbInfo -> partsLabels[9] = (char*) _idParamsBGam3_;
    _fishpbInfo -> partsExist[9] = (bool*) &_true_;
    _fishpbInfo -> partsMult[9] = (int*) _multParamsBGam3_;
    _fishpbInfo -> partsDerivLog[9] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[9] = NULL;

    /* RSD */
    _fishpbInfo -> partsLabels[10] = (char*) _idParamsF_;
    _fishpbInfo -> partsExist[10] = (bool*) &_true_;
    _fishpbInfo -> partsMult[10] = (int*) _multParamsF_;
    _fishpbInfo -> partsDerivLog[10] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[10] = NULL;

    _fishpbInfo -> partsLabels[11] = (char*) _idParamsSigv_;
    _fishpbInfo -> partsExist[11] = (bool*) &_true_;
    _fishpbInfo -> partsMult[11] = (int*) _multParamsSigv_;
    _fishpbInfo -> partsDerivLog[11] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[11] = NULL;

    /* AP */
    _fishpbInfo -> partsLabels[12] = (char*) _idParamsD_;
    _fishpbInfo -> partsExist[12] = (bool*) &_true_;
    _fishpbInfo -> partsMult[12] = (int*) _multParamsD_;
    _fishpbInfo -> partsDerivLog[12] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[12] = NULL;

    _fishpbInfo -> partsLabels[13] = (char*) _idParamsH_;
    _fishpbInfo -> partsExist[13] = (bool*) &_true_;
    _fishpbInfo -> partsMult[13] = (int*) _multParamsH_;
    _fishpbInfo -> partsDerivLog[13] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[13] = NULL;

    /* Ctr */
    _fishpbInfo -> partsLabels[14] = (char*) _idParamsC0_;
    _fishpbInfo -> partsExist[14] = (bool*) &_true_;
    _fishpbInfo -> partsMult[14] = (int*) _multParamsC0_;
    _fishpbInfo -> partsDerivLog[14] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[14] = NULL;

    _fishpbInfo -> partsLabels[15] = (char*) _idParamsC2_;
    _fishpbInfo -> partsExist[15] = (bool*) &_true_;
    _fishpbInfo -> partsMult[15] = (int*) _multParamsC2_;
    _fishpbInfo -> partsDerivLog[15] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[15] = NULL;

    _fishpbInfo -> partsLabels[16] = (char*) _idParamsC4_;
    _fishpbInfo -> partsExist[16] = (bool*) &_true_;
    _fishpbInfo -> partsMult[16] = (int*) _multParamsC4_;
    _fishpbInfo -> partsDerivLog[16] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[16] = NULL;

    /* Shot noise */
    _fishpbInfo -> partsLabels[17] = (char*) _idParamsPsn_;
    _fishpbInfo -> partsExist[17] = (bool*) &_true_;
    _fishpbInfo -> partsMult[17] = (int*) _multParamsPsn_;
    _fishpbInfo -> partsDerivLog[17] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[17] = NULL;

    _fishpbInfo -> partsLabels[18] = (char*) _idParamsBsn1_;
    _fishpbInfo -> partsExist[18] = (bool*) &_true_;
    _fishpbInfo -> partsMult[18] = (int*) _multParamsBsn1_;
    _fishpbInfo -> partsDerivLog[18] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[18] = NULL;

    _fishpbInfo -> partsLabels[19] = (char*) _idParamsBsn2_;
    _fishpbInfo -> partsExist[19] = (bool*) &_true_;
    _fishpbInfo -> partsMult[19] = (int*) _multParamsBsn2_;
    _fishpbInfo -> partsDerivLog[19] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[19] = NULL;

    /* Linear power spectrum */
    _fishpbInfo -> partsLabels[20] = (char*) _idParamsPk_;
    _fishpbInfo -> partsExist[20] = (bool*) &_true_;
    _fishpbInfo -> partsMult[20] = (int*) _multParamsPk_;
    _fishpbInfo -> partsDerivLog[20] = (bool*) &_false_;
    _fishpbInfo -> partsDerivVals[20] = NULL;


    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Free Local Variables   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_info(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int fish_free(void)
{
    /*

        Free local variables

    */

    /* Free info */
    _free_info();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_info(void)
{
    /*

        Free the info structs

    */

    /* PP fisher matrix */
    _fishppInfo = _fish_info_free(_fishppInfo);

    /* BB fisher matrix */
    _fishbbInfo = _fish_info_free(_fishbbInfo);

    /* PB fisher matrix */
    _fishpbInfo = _fish_info_free(_fishpbInfo);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     FISHER MATRIX STRUCTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   Covariance Matrix Output Struct   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fish_out_t *fish_out_new(void)
{
    /*

        Create a new fish_out_t struct

    */

    fish_out_t *out = malloc(sizeof(fish_out_t));

    out -> size = 0;
    out -> labels = NULL;
    out -> log = NULL;

    out -> file = NULL;

    out -> binary = _true_;
    out -> precision = 0;

    return out;
}


fish_out_t *fish_out_free(fish_out_t *out)
{
    /*

        Free out

    */

    /* Check for NULL */
    if (out == NULL)
        return NULL;

    /* Free contents */
    for (size_t i = 0; i < out -> size; i++)
        free(out -> labels[i]);

    free(out -> labels);

    free(out -> log);
    free(out -> file);

    /* Free out itself */
    free(out);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


fish_out_t *fish_out_cp(fish_out_t *out)
{
    /*

        Copy out

    */

    fish_out_t *outCp = malloc(sizeof(fish_out_t));

    outCp -> size = out -> size;
    outCp -> labels = malloc(sizeof(char*) * out -> size);
    outCp -> log = malloc(sizeof(bool) * out -> size);

    for (size_t i = 0; i < out -> size; i++)
      {
        outCp -> labels[i] = misc_scat(1, out -> labels[i]);
        outCp -> log[i] = out -> log[i];
      }

    outCp -> file = misc_scat(1, out -> file);

    outCp -> binary = out -> binary;
    outCp -> precision = out -> precision;

    return outCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int fish_out_add_label(fish_out_t *out, const char *label)
{
    /*

        Add an output fisher matrix element label to out

    */

    /* Label already contained in out -> labels */
    if (misc_ssearch(out -> labels, out -> size, label, NULL))
        return 0;

    /* Increment the output size and reallocate memory */
    out -> size += 1;
    out -> labels = realloc(out -> labels, sizeof(char*) * out -> size);
    out -> log = realloc(out -> log, sizeof(bool) * out -> size);

    /* Add the label */
    out -> labels[out -> size - 1] = misc_scat(1, label);

    /* Set the logarithm flag to _false_ */
    out -> log[out -> size - 1] = _false_;

    return 0;
}


int fish_out_rm_label(fish_out_t *out, const char *label)
{
    /*

        Remove an output fisher matrix element label from out

    */

    bool success;

    /* Get the index */
    size_t index = misc_ssearch_index(out -> labels, out -> size, label, &success);

    /* Not in the array */
    if (!success)
      {
        printf("Couldn't remove the label %s from the fisher matrix element labels list as it did not exist.\n", label);

        return 0;
      }

    /* Decrease the size of the array */
    out -> size -= 1;

    /* Push the string to the end of the array */
    for (size_t i = index; i < out -> size; i++)
      {
        misc_swap(&(out -> labels)[i], &(out -> labels)[i+1], "s");
        misc_swap(&(out -> log)[i], &(out -> log)[i+1], "b");
      }

    /* Free the last label */
    free(out -> labels[out -> size]);

    /* Reallocate memory */
    out -> labels = realloc(out -> labels, sizeof(char*) * out -> size);
    out -> log = realloc(out -> log, sizeof(bool*) * out -> size);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int fish_out_set_log(fish_out_t *out, const char *label, bool log)
{
    /*

        Set the logarithm flag for the output label in out

    */

    bool success;

    /* Get the index */
    size_t index = misc_ssearch_index(out -> labels, out -> size, label, &success);

    /* Not in the array */
    if (!success)
      {
        printf("Cannot set the logarithmic parameter flag for non-existent matrix element label '%s' in a fish_out_t struct!\n", label);

        return 0;
      }

    out -> log[index] = log;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int fish_out_set_file(fish_out_t *out, const char *file)
{
    /*

        Set the file name for out (if file is NULL nothing will be output to file)

    */

    misc_scp(&(out -> file), file);

    return 0;
}


int fish_out_set_binary(fish_out_t *out, bool binary)
{
    /*

        Set the binary flag for out

    */

    out -> binary = binary;

    return 0;
}


int fish_out_set_precision(fish_out_t *out, int precision)
{
    /*

        Set the precision in the output file for out -> also set the binary flag to _false_

    */

    out -> binary = _false_;
    out -> precision = precision;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Fisher Matrix Struct   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fish_mat_t *fish_mat_new(const char *id)
{
    /*

        New fish_mat_t struct

    */

    fish_mat_t *fishMat = malloc(sizeof(fish_mat_t));

    /* Info struct */
    _fish_info_t *info = _fish_info_get_struct(id);

    /* Id */
    fishMat -> id = info -> id;

    /* Output */
    fishMat -> out = NULL;

    /* Print */
    fishMat -> flags = NULL;


    return fishMat;
}


/*  ------------------------------------------------------------------------------------------------------  */


fish_mat_t *fish_mat_free(fish_mat_t *fishMat)
{
    /*

        Free a fish_mat_t struct

    */

    /* Check for NULL */
    if (fishMat == NULL)
        return NULL;


    /* Free fishMat's contents */

    /* Out */
    fishMat -> out = fish_out_free(fishMat -> out);

    /* Flags */
    fishMat -> flags = print_flags_free(fishMat -> flags);


    /* Free fishMat itself */

    free(fishMat);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int fish_mat_set_out(fish_mat_t *fishMat, fish_out_t *out)
{
    /*

        Set the output parameters for spec

    */

    fishMat -> out = fish_out_free(fishMat -> out);
    fishMat -> out = fish_out_cp(out);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int fish_mat_set_flags(fish_mat_t *fishMat, print_flags_t *flags)
{
    /*

        Set the print parameters for spec

    */

    fishMat -> flags = print_flags_free(fishMat -> flags);
    fishMat -> flags = print_flags_cp(flags);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fish_out_t *fish_mat_get_out(fish_mat_t *fishMat)
{
    /*

        Get the out struct from fishMat

    */

    return fishMat -> out;
}


print_flags_t *fish_mat_get_flags(fish_mat_t *fishMat)
{
    /*

        Get the print struct from fishMat

    */

    return fishMat -> flags;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool _fish_mat_test(fish_mat_t *fishMat)
{
    /*

        Test if fishMat can be used in fish_mat_poly

    */

    /* Info struct */
    _fish_info_t *info = _fish_info_get_struct(fishMat -> id);

    /* Spatial variables */
    for (size_t i = 0; i < info -> specSize; i++)
      {
        sample_shape_t *sampleShape1 = flss_get_sample_shape(info -> specLabels[i]);
        sample_arg_t *sampleArgK1 = sampleShape1 -> sampleRawLength -> sampleArg;
        sample_arg_t *sampleArgMu1 = sampleShape1 -> sampleRawOrientation -> sampleArg;

        for (size_t j = i + 1; j < info -> specSize; j++)
          {
            sample_shape_t *sampleShape2 = flss_get_sample_shape(info -> specLabels[j]);
            sample_arg_t *sampleArgK2 = sampleShape2 -> sampleRawLength -> sampleArg;
            sample_arg_t *sampleArgMu2 = sampleShape2 -> sampleRawOrientation -> sampleArg;

            if (!sample_arg_compare_sub(sampleArgK1, sampleArgK2) || !sample_arg_compare_sub(sampleArgMu1, sampleArgMu2))
              {
                printf("The raw sampling of the variables 'k' and 'mu' for the '%s' and '%s' spectra must be atleast subsets of each other.\n", info -> specLabels[i], info -> specLabels[j]);
                exit(1);

                return _false_;
              }
          }
      }

    return _true_;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Default Settings   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fish_mat_t *fish_mat_pp_default(void)
{
    /*

        Default parameters for the pp fisher matrix function

    */

    fish_mat_t *fishMat = fish_mat_new(_idFishPP_);


    /* Output */

    fish_out_t *out = fish_out_new();

    /* Labels and logarithmic flags */
    for (size_t i = 0; i < _fishppInfo -> partsSize; i++)
      {
        fish_out_add_label(out, _fishppInfo -> partsLabels[i]);
        fish_out_set_log(out, _fishppInfo -> partsLabels[i], _false_);
      }

    /* File */
    char *outFile = misc_scat(2, _fishppInfo -> id, _outExt);
    fish_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    fish_out_set_precision(out, 16);
    fish_out_set_binary(out, _true_);

    fish_mat_set_out(fishMat, out);

    out = fish_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    fish_mat_set_flags(fishMat, flags);

    flags = print_flags_free(flags);


    return fishMat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fish_mat_t *fish_mat_bb_default(void)
{
    /*

        Default parameters for the bb fisher matrix function

    */

    fish_mat_t *fishMat = fish_mat_new(_idFishBB_);


    /* Output */

    fish_out_t *out = fish_out_new();

    /* Labels and errors */
    for (size_t i = 0; i < _fishbbInfo -> partsSize; i++)
      {
        fish_out_add_label(out, _fishbbInfo -> partsLabels[i]);
        fish_out_set_log(out, _fishbbInfo -> partsLabels[i], _false_);
      }

    /* File */
    char *outFile = misc_scat(2, _fishbbInfo -> id, _outExt);
    fish_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    fish_out_set_precision(out, 16);
    fish_out_set_binary(out, _true_);

    fish_mat_set_out(fishMat, out);

    out = fish_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    fish_mat_set_flags(fishMat, flags);

    flags = print_flags_free(flags);


    return fishMat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


fish_mat_t *fish_mat_pb_default(void)
{
    /*

        Default parameters for the pb fisher matrix function

    */

    fish_mat_t *fishMat = fish_mat_new(_idFishPB_);


    /* Output */

    fish_out_t *out = fish_out_new();

    /* Labels and logarithmic flags */
    for (size_t i = 0; i < _fishpbInfo -> partsSize; i++)
      {
        fish_out_add_label(out, _fishpbInfo -> partsLabels[i]);
        fish_out_set_log(out, _fishpbInfo -> partsLabels[i], _false_);
      }

    /* File */
    char *outFile = misc_scat(2, _fishpbInfo -> id, _outExt);
    fish_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    fish_out_set_precision(out, 16);
    fish_out_set_binary(out, _true_);

    fish_mat_set_out(fishMat, out);

    out = fish_out_free(out);


    /* Flags */

    print_flags_t *flags = print_flags_new();

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);

    fish_mat_set_flags(fishMat, flags);

    flags = print_flags_free(flags);


    return fishMat;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CALCULATE COVARIANCE MATRICES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Setup Matrix Data Struct   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_fish_setup_mat(fish_mat_t *fishMat)
{
    /*

        Setup the output mat_t struct

    */

    /* Info struct */
    _fish_info_t *info = _fish_info_get_struct(fishMat -> id);


    /* Variables */

    /* Temporal */
    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    /* Spatial */
    sample_raw_t *sampleRawK = sample_raw_new();
    sample_raw_t *sampleRawMu = sample_raw_new();

    sample_arg_t *sampleArgK = sample_arg_new();
    sample_arg_t *sampleArgMu = sample_arg_new();

    for (size_t i = 0; i < info -> specSize; i++)
      {
        sample_shape_t *sampleShape_ = flss_get_sample_shape(info -> specLabels[i]);

        sample_raw_t *sampleRawK_ = sampleShape_ -> sampleRawLength;
        sample_raw_t *sampleRawMu_ = sampleShape_ -> sampleRawOrientation;

        sample_arg_t *sampleArgK_ = sampleRawK_ -> sampleArg;
        sample_arg_t *sampleArgMu_ = sampleRawMu_ -> sampleArg;

        for (size_t n = 0; n < sampleArgK_ -> size; n++)
            misc_insort(&sampleRawK -> array, &sampleArgK -> size, sampleRawK_ -> array[n], __ABSTOL__, NULL);

        for (size_t n = 0; n < sampleArgMu_ -> size; n++)
            misc_insort(&sampleRawMu -> array, &sampleArgMu -> size, sampleRawMu_ -> array[n], __ABSTOL__, NULL);
      }


    /* Sample as data struct */

    /* z sample as dat struct */
    dat_t *sampleZDat = dat_new(1, 0, sampleArgZ -> size);
    dat_set_label(sampleZDat, 0, 'x', (char*) _idVarZ_);

    for (size_t i = 0; i < sampleArgZ -> size; i++)
        dat_set_value(sampleZDat, 0, i, 'x', sampleRawZ -> array[i]);

    /* k sample as dat struct */
    dat_t *sampleKDat = dat_new(1, 0, sampleArgK -> size);
    dat_set_label(sampleKDat, 0, 'x', (char*) _idVarK_);

    for (size_t i = 0; i < sampleArgK -> size; i++)
        dat_set_value(sampleKDat, 0, i, 'x', sampleRawK -> array[i]);

    /* mu sample as dat struct */
    dat_t *sampleMuDat = dat_new(1, 0, sampleArgMu -> size);
    dat_set_label(sampleMuDat, 0, 'x', (char*) _idVarMu_);

    for (size_t i = 0; i < sampleArgMu -> size; i++)
        dat_set_value(sampleMuDat, 0, i, 'x', sampleRawMu -> array[i]);


    /* Parameters */

    /* Dimension of the matrix */
    size_t fishDim = 0;

    for (size_t i = 0; i < info -> partsSize; i++)
      {
        info -> partsExist[i] = (bool*) &_false_;
        info -> partsDerivLog[i] = (bool*) &_false_;
        info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);

        dat_t *derivVals = NULL;

        bool success;
        size_t index = misc_ssearch_index(fishMat -> out -> labels, fishMat -> out -> size, info -> partsLabels[i], &success);

        /* Parameter not requested */
        if (!success)
            continue;

        /* Parameter exists */
        info -> partsExist[i] = (bool*) &_true_;

        /* Logarithmic derivative */
        info -> partsDerivLog[i] = (fishMat -> out -> log[index]) ? (bool*) &_true_ : (bool*) &_false_;

        size_t dim = 1;

        /* z multiplicity */
        if (info -> partsMult[i][0] != 0)
          {
            dim *= sampleArgZ -> size;

            derivVals = dat_cat(2, info -> partsDerivVals[i], sampleZDat);
            info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);
            info -> partsDerivVals[i] = dat_cat(1, derivVals);
            derivVals = dat_free(derivVals);
          }

        /* k multiplicity */
        if (info -> partsMult[i][1] != 0)
          {
            dim *= sampleArgK -> size;

            derivVals = dat_cat(2, info -> partsDerivVals[i], sampleKDat);
            info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);
            info -> partsDerivVals[i] = dat_cat(1, derivVals);
            derivVals = dat_free(derivVals);
          }

        /* mu multiplicity */
        if (info -> partsMult[i][2] != 0)
          {
            dim *= sampleArgMu -> size;

            derivVals = dat_cat(2, info -> partsDerivVals[i], sampleMuDat);
            info -> partsDerivVals[i] = dat_free(info -> partsDerivVals[i]);
            info -> partsDerivVals[i] = dat_cat(1, derivVals);
            derivVals = dat_free(derivVals);
          }

        /* Add the dimensions to the total dimensions */
        fishDim += dim;
      }

    size_t dim[2] = {fishDim, fishDim};
    mat_t *mat = mat_new("f", dim, _false_);
    mat_set_label(mat, (char*) info -> id);

    /* Free memory */
    sampleRawK = sample_raw_free(sampleRawK);
    sampleRawMu = sample_raw_free(sampleRawMu);

    sampleArgK = sample_arg_free(sampleArgK);
    sampleArgMu = sample_arg_free(sampleArgMu);

    sampleZDat = dat_free(sampleZDat);
    sampleKDat = dat_free(sampleKDat);
    sampleMuDat = dat_free(sampleMuDat);

    return mat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Fisher Matrix Function   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *fish_mat_poly(fish_mat_t *fishMat)
{
    /*

        Calculate the fisher matrix

    */


    /* Exit if nothing should be calculated */

    /* fishMat itself is NULL */
    if (fishMat == NULL)
      {
        return NULL;
      }

    /* Outsize is 0 */
    if (fishMat -> out -> size == 0)
      {
        return NULL;
      }

    /* fishMat is faulty */
    if (!_fish_mat_test(fishMat))
      {
        return NULL;
      }


    // TODO: Clean this mess up...

    /* Get parameters from fishMat */

    /* Info on spectrum */
    _fish_info_t *info = _fish_info_get_struct(fishMat -> id);

    /* Kern order */
    size_t kernOrder = fish_info_get_kern_order(info -> id);

    /* Get the output matrix struct */
    mats_t *mats = mats_new(1);
    mat_t *mat = _fish_setup_mat(fishMat);
    mats_set_mat(mats, 0, mat);

    /* Variables */
    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    sample_arg_t *sampleArgK = ((sample_raw_t*) flss_get_sample_shape(info -> specLabels[0]) -> sampleRawLength) -> sampleArg;
    sample_arg_t *sampleArgMu = ((sample_raw_t*) flss_get_sample_shape(info -> specLabels[0]) -> sampleRawOrientation) -> sampleArg;

    /* Shared variables */
    __float128 valFish;

    /* Covariance matrix */
    size_t covDimT[2] = {sampleArgZ -> size, sampleArgZ -> size};
    mat_t *covInv = mat_new("d", covDimT, _true_);

    mat_reduce_t red = {1.e-14}; // TODO: Improve this...

    for (size_t locT[2] = {0, 0}; locT[0] < sampleArgZ -> size; locT[1] = ++locT[0])
      {
        /* Put the single covariance matrices into blocks at each redshift */
        size_t covDimTBlock[2] = {info -> specSize, info -> specSize};
        mat_t *covBlock = mat_new("s", covDimTBlock, _true_);

        for (size_t loc[2] = {0, 0}, index = 0; loc[0] < info -> specSize; loc[0] = (loc[1] < info -> specSize - 1) ? loc[0] : loc[0] + 1, loc[1] = (loc[1] < info -> specSize - 1) ? loc[1] + 1 : loc[0], index++)
          {
            cov_mat_t *covMat = cov_mat_new(info -> covLabels[index]);
            cov_out_add_label(cov_mat_get_out(covMat), cov_in_get_label(info -> covLabels[index]));

            /* Get the block */
            mat_t *covBlockSub = mat_get_block(_covPolyMat_(info -> covLabels[index])(covMat, NULL), locT);

            /* Set the block */
            mat_fset_block(covBlock, loc, covBlockSub);

            /* Free memory */
            covMat = cov_mat_free(covMat);
          }

        /* Reduce the matrix */
        covBlock = mat_reduce(covBlock, &red);

        /* Get the inverse of the covariance matrix */
        mat_t *covBlockInv = mat_inv_cov(mat_inv_lapack, covBlock, NULL);

        /* Set the block */
        mat_fset_block(covInv, locT, covBlockInv);

        /* Free memory */
        covBlock = mat_free(covBlock);
      }

    /* Sizes */
    size_t sizes[5] = {sampleArgZ -> size, sampleArgK -> size, sampleArgMu -> size, mat -> dim[0], mat -> dim[0] * (mat -> dim[1] + 1) / 2};

    /* Print struct */
    print_t *print = print_new();

    print_set_id(print, fishMat -> id);
    print_set_flags(print, fishMat -> flags);

    print_set_sizes(print, sizes);

    /* Store derivatives */
    dat_t **datDeriv = malloc(sizeof(dat_t*) * info -> specSize);

    for (size_t i = 0; i < info -> specSize; i++)
      {
        sample_shape_t *sampleShape = flss_get_sample_shape(info -> specLabels[i]);

        datDeriv[i] = dat_new(0, mat -> dim[0], sampleArgZ -> size * sampleShape -> size);

        for (size_t n = 0; n < datDeriv[i] -> size; n++)
          {
            for (size_t m = 0; m < datDeriv[i] -> yDim; m++)
                dat_set_yvalue(datDeriv[i], m, n, NAN);
          }
      }


    /* Parallel threading */

    #pragma omp parallel shared(valFish)

      { // Start pragma parallel

        /* Declare and define variables */

        /* Spec args struct */
        spec_arg_t *specArg1 = spec_arg_new(NULL);
        spec_arg_t *specArg2 = spec_arg_new(NULL);

        spec_arg_set_dz(specArg1, sampleArgZ -> step);
        spec_arg_set_dz(specArg2, sampleArgZ -> step);

        spec_arg_set_dk(specArg1, sampleArgK -> step);
        spec_arg_set_dk(specArg2, sampleArgK -> step);

        spec_arg_set_dmu(specArg1, sampleArgMu -> step);
        spec_arg_set_dmu(specArg2, sampleArgMu -> step);

        /* Kern struct */
        kern_t *kern1 = kernels_new_order(info -> specOrders[info -> specSize - 1], kernOrder);
        kern_t *kern2 = kernels_new_order(info -> specOrders[info -> specSize - 1], kernOrder);

        spec_arg_set_kern(specArg1, kern1);
        spec_arg_set_kern(specArg2, kern2);

        /* Deriv struct */
        spec_deriv_t *deriv1 = spec_deriv_new();
        spec_deriv_t *deriv2 = spec_deriv_new();

        spec_arg_set_deriv(specArg1, deriv1);
        spec_arg_set_deriv(specArg2, deriv2);

        /* Print stage */
        #pragma omp single
          {
            print_fish(print);
          }

        /* Print progress */
        #pragma omp single
          {
            print_fish(print);
          }


        /* Calculate the fisher matrix */

        for (size_t n = 0, m = 0, f1 = 0, f2 = 0; n < info -> partsSize; n = (m < info -> partsSize - 1) ? n : n + 1, m = (m < info -> partsSize - 1) ? m + 1 : n)

          { // Start info -> partsSize for

            /* Skip if one element does not exist */
            if (!*(info -> partsExist[n]) || !*(info -> partsExist[m]))
                continue;

            /* Set the logarithmic flags */
            spec_deriv_set_log(deriv1, *(info -> partsDerivLog[n]));
            spec_deriv_set_log(deriv2, *(info -> partsDerivLog[m]));

            for (size_t r = 0, s = 0; r < info -> partsDerivVals[n] -> size; r = (s < info -> partsDerivVals[m] -> size - 1) ? r : r + 1, s = (s < info -> partsDerivVals[m] -> size - 1) ? s + 1 : ( (n == m) ? r : 0 ))

              { // Start info -> partsDerivVals -> size for TODO: What about constant parameters?

                /* Set the derivative variables for the first spectrum */
                for (size_t i = 0; i < info -> partsDerivVals[n] -> xDim; i++)
                  {
                    spec_deriv_set_var(deriv1, dat_get_label(info -> partsDerivVals[n], i, 'x'), dat_get_value(info -> partsDerivVals[n], i, r, 'x'));
                  }

                /* Set the derivative variables for the second spectrum */
                for (size_t i = 0; i < info -> partsDerivVals[m] -> xDim; i++)
                  {
                    spec_deriv_set_var(deriv2, dat_get_label(info -> partsDerivVals[m], i, 'x'), dat_get_value(info -> partsDerivVals[m], i, s, 'x'));
                  }

                /* Fisher matrix element */
                #pragma omp single
                valFish = 0.q;

                for (size_t t = 0; t < sampleArgZ -> size; t++)

                  { // Start temporal for

                    /* Skip partial multiplicities with unequal redshifts */
                    if ((info -> partsMult[n][0] == -1 && deriv1 -> z != sampleRawZ -> array[t]) || (info -> partsMult[m][0] == -1 && deriv2 -> z != sampleRawZ -> array[t]))
                        continue;

                    /* Set the redshift + fiducials in the kern structs */
                    kernels_set_z(kern1, sampleRawZ -> array[t]);
                    kernels_set_z(kern2, sampleRawZ -> array[t]);

                    /* Temporal block of the covariance matrix */
                    size_t locT[2] = {t, t};
                    mat_t *covInvT = mat_fget_block(covInv, locT);

                    for (size_t i1 = 0, i2 = 0; i1 < info -> specSize; i1 = (i2 < info -> specSize - 1) ? i1 : i1 + 1, i2 = (i2 < info -> specSize - 1) ? i2 + 1 : i1)

                      { // Start info -> specSize for

                        /* Variables */
                        sample_shape_t *sampleShape1 = flss_get_sample_shape(info -> specLabels[i1]);
                        sample_shape_t *sampleShape2 = flss_get_sample_shape(info -> specLabels[i2]);

                        sample_raw_t *sampleRawK1 = sampleShape1 -> sampleRawLength;
                        sample_raw_t *sampleRawK2 = sampleShape2 -> sampleRawLength;

                        sample_arg_t *sampleArgK1 = sampleRawK1 -> sampleArg;
                        sample_arg_t *sampleArgK2 = sampleRawK2 -> sampleArg;

                        /* Skip if one spectrum's deriv -> k is not included in the sample */
                        if (info -> partsMult[n][1] == 1 && !misc_bsearch(sampleRawK1 -> array, sampleArgK1 -> size, deriv1 -> k, __ABSTOL__, NULL))
                            continue;

                        if (info -> partsMult[m][1] == 1 && !misc_bsearch(sampleRawK2 -> array, sampleArgK2 -> size, deriv2 -> k, __ABSTOL__, NULL))
                            continue;

                        /* Spectral derivative functions */
                        double (*dSpecFunc1)(void*, void*) = _specPoly_(info -> specLabels[i1], info -> partsLabels[n]);
                        double (*dSpecFunc2)(void*, void*) = _specPoly_(info -> specLabels[i2], info -> partsLabels[m]);

                        /* Average flags */
                        bool avrFlag1 = flss_get_avr_shape_flag(info -> specLabels[i1]);
                        bool avrFlag2 = flss_get_avr_shape_flag(info -> specLabels[i2]);

                        /* Average functions */
                        int (*avrFunc1)(double (*)(void*, void*), spec_arg_t *, double*) = (avrFlag1) ? avr_shape_get_func(info -> specOrders[i1]) : avr_shape_inf;
                        int (*avrFunc2)(double (*)(void*, void*), spec_arg_t *, double*) = (avrFlag2) ? avr_shape_get_func(info -> specOrders[i2]) : avr_shape_inf;

                        /* Spectral block of the covariance matrix */
                        size_t locTBlock[2] = {i1, i2};
                        mat_t *covInvTBlock = mat_fget_block(covInvT, locTBlock);

                        /* Block does not exist */
                        if (covInvTBlock == NULL)
                            continue;

                        #pragma omp for schedule(dynamic)

                        for (size_t s1 = 0; s1 < sampleShape1 -> size; s1++)

                          { // Start spatial for 1

                            /* (Spatial) Location in the inverse covariance matrix */
                            size_t locS[2];

                            /* Scalar product between inverse covariance matrix and derivative vector */
                            __float128 valFishTemp = 0.q;

                            /* Shape of first spectrum */
                            shape_t *shape1 = sampleShape1 -> arrayShape[s1];

                            /* Set variables for first shape */
                            for (size_t i = 0; i < shape1 -> dim; i++)
                              {
                                kernels_qset_k(kern1, i, shape1 -> length[i]);
                                kernels_qset_mu(kern1, i, shape1 -> orientation[i]);

                                for (size_t j = i + 1; j < shape1 -> dim; j++)
                                    kernels_qset_nu(kern1, i, j, shape1 -> angle[shape_get_vertex_angle_index(shape1 -> dim, i, j)]);
                              }

                            /* Spectra derivative */
                            double dSpecVal1 = dat_get_yvalue(datDeriv[i1], f1, t * sampleShape1 -> size + s1);

                            if (isnan(dSpecVal1))
                              {
                                dSpecVal1 = avr_shape_direct(avrFunc1, dSpecFunc1, specArg1);
                                dat_set_yvalue(datDeriv[i1], f1, t * sampleShape1 -> size + s1, dSpecVal1);
                              }

                            for (size_t p1 = 0; p1 < 1 + (size_t) (!shape1 -> parity); p1++)

                              { // Start parity for 1

                                /* (Spatial) Row dimension in the inverse covariance matrix */
                                locS[0] = (s1 < sampleShape1 -> sizeParity) ? s1 : sampleShape1 -> sizeParity + 2 * (s1 - sampleShape1 -> sizeParity) + p1;

                                /* Boundaries */
                                size_t spatialBounds2[2];
                                spatialBounds2[0] = (!strcmp(covInvTBlock -> mtype -> id, "d")) ? s1 : 0;
                                spatialBounds2[1] = (!strcmp(covInvTBlock -> mtype -> id, "d")) ? s1 + 1 : sampleShape2 -> size;

                                for (size_t s2 = spatialBounds2[0]; s2 < spatialBounds2[1]; s2++)

                                  { // Start spatial for 2

                                    /* Shape of second spectrum */
                                    shape_t *shape2 = sampleShape2 -> arrayShape[s2];

                                    /* Set variables for second shape */
                                    for (size_t i = 0; i < shape2 -> dim; i++)
                                      {
                                        kernels_qset_k(kern2, i, shape2 -> length[i]);
                                        kernels_qset_mu(kern2, i, shape2 -> orientation[i]);

                                        for (size_t j = i + 1; j < shape2 -> dim; j++)
                                            kernels_qset_nu(kern2, i, j, shape2 -> angle[shape_get_vertex_angle_index(shape2 -> dim, i, j)]);
                                      }

                                    /* Spectra derivative */
                                    double dSpecVal2 = dat_get_yvalue(datDeriv[i2], f2, t * sampleShape2 -> size + s2);

                                    if (isnan(dSpecVal2))
                                      {
                                        dSpecVal2 = avr_shape_direct(avrFunc2, dSpecFunc2, specArg2);
                                        dat_set_yvalue(datDeriv[i2], f2, t * sampleShape2 -> size + s2, dSpecVal2);
                                      }

                                    for (size_t p2 = 0; p2 < 1 + (size_t) (!shape2 -> parity); p2++)

                                      { // Start parity for 2

                                        /* (Spatial) Column dimension in the inverse covariance matrix */
                                        locS[1] = (s2 < sampleShape2 -> sizeParity) ? s2 : sampleShape2 -> sizeParity + 2 * (s2 - sampleShape2 -> sizeParity) + p2;

                                        /* Spatial element of the covariance matrix */
                                        double valCovInv = mat_get_value(covInvTBlock, locS);

                                        /* No need to calculate the spectra */
                                        if (valCovInv == 0.)
                                            continue;

                                        /* Calculate the fisher matrix element */
                                        valFishTemp += ((__float128) valCovInv) * ((__float128) dSpecVal2);

                                      } // End parity for 2

                                  } // End spatial for 2

                              } // End parity for 1

                            /* Update the Fisher matrix element */
                            #pragma omp atomic
                            valFish += ((__float128) dSpecVal1) * valFishTemp;

                          } // End spatial for 1

                        /* Free memory */
                        covInvTBlock = mat_ffree(covInvTBlock);

                      } // End info -> specSize for

                    /* Free memory */
                    covInvT = mat_ffree(covInvT);

                  } // End temporal for

                #pragma omp single
                  {
                    /* Insert the value */
                    size_t locFish1[2] = {f1, f2};
                    size_t locFish2[2] = {f2, f1};

                    mat_set_value(mat, locFish1, (double) valFish);
                    mat_set_value(mat, locFish2, (double) valFish);

                    /* Update the progress */
                    print_update_progress(print, 1. / (double) (print -> sizes[4]));
                    print_fish(print);
                  }

                /* Update the fisher indices */
                if (r < info -> partsDerivVals[n] -> size - 1 || s < info -> partsDerivVals[m] -> size - 1)
                  {
                    f1 = (s < info -> partsDerivVals[m] -> size - 1) ? f1 : f1 + 1;
                    f2 = (s < info -> partsDerivVals[m] -> size - 1) ? f2 + 1 : ( (n == m) ? f1 : f2 + 1 - info -> partsDerivVals[m] -> size );
                  }

              } // End info -> partsDerivVals -> size for

            /* Update the fisher indices */
            f1 = (f2 < mat -> dim[1] - 1) ? f1 + 1 - info -> partsDerivVals[n] -> size : f1 + 1;
            f2 = (f2 < mat -> dim[1] - 1) ? f2 + 1 : f1;

          } // End info -> partsSize for


        /* Stage has finished: print a message and the execution time */

        print -> stageFinished = 1;

        #pragma omp single
          {
            print_set_stagefinished(print, 1);
            print_set_progress(print, 1.);

            print_fish(print);
          }


        /* Free memory */

        specArg1 = spec_arg_free(specArg1);
        specArg2 = spec_arg_free(specArg2);


      } // End pragma parallel


    /* Free memory */

    covInv = mat_free(covInv);

    for (size_t i = 0; i < info -> specSize; i++)
      {
        datDeriv[i] = dat_free(datDeriv[i]);
      }

    free(datDeriv);

    print = print_free(print);


    /* Write the result to file */

    if (fishMat -> out -> file != NULL && fishMat -> out -> precision != 0)
      {
        char *outFile = misc_scat(3, __FISHERLSSDIR__, _outDir, fishMat -> out -> file);

        printf("%s\n", outFile);

        if (fishMat -> out -> binary)
            mats_output(outFile, mats, NULL);

        else
            mats_output(outFile, mats, &fishMat -> out -> precision);

        free(outFile);
      }

    return mats;
}


mats_t *fish_mat_poly_eig(fish_mat_t *fishMat)
{
    /*

        Calculate the fisher matrix

    */


    /* Exit if nothing should be calculated */

    /* fishMat itself is NULL */
    if (fishMat == NULL)
      {
        return NULL;
      }

    /* Outsize is 0 */
    if (fishMat -> out -> size == 0)
      {
        return NULL;
      }

    /* fishMat is faulty */
    if (!_fish_mat_test(fishMat))
      {
        return NULL;
      }


    // TODO: Clean this mess up...

    /* Get parameters from fishMat */

    /* Info on spectrum */
    _fish_info_t *info = _fish_info_get_struct(fishMat -> id);

    /* Kern order */
    size_t kernOrder = fish_info_get_kern_order(info -> id);

    /* Get the output matrix struct */
    mats_t *mats = mats_new(1);
    mat_t *mat = _fish_setup_mat(fishMat);
    mats_set_mat(mats, 0, mat);

    /* Variables */
    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    sample_arg_t *sampleArgK = ((sample_raw_t*) flss_get_sample_shape(info -> specLabels[0]) -> sampleRawLength) -> sampleArg;
    sample_arg_t *sampleArgMu = ((sample_raw_t*) flss_get_sample_shape(info -> specLabels[0]) -> sampleRawOrientation) -> sampleArg;

    /* Eigenvalue and eigenvector matrices */
    size_t covDimT[2] = {sampleArgZ -> size, sampleArgZ -> size};

    mat_t *covDiag = mat_new("d", covDimT, true);
    mat_t *corrEigval = mat_new("d", covDimT, true);
    mat_t *corrEigvec = mat_new("d", covDimT, true);

    /* Derivatives times eigenvectors (dP x U^T) */
    size_t dimDeriv[2] = {sampleArgZ -> size, 1};
    mat_t *matDeriv = mat_new("f", dimDeriv, true);

    mat_reduce_t red = {1.e-14}; // TODO: Improve this...

    for (size_t locT[2] = {0, 0}; locT[0] < sampleArgZ -> size; locT[1] = ++locT[0])
      {
        /* Put the single dP x U^T matrices into blocks at each redshift */
        size_t dimDerivBlock[2] = {info -> specSize, mat -> dim[0]};
        mat_t *matDerivBlock = mat_new("f", dimDerivBlock, true);

        for (size_t loc[2] = {0, 0}; loc[0] < info -> specSize; loc[0] = (loc[1] < mat -> dim[0] - 1) ? loc[0] : loc[0] + 1, loc[1] = (loc[1] < mat -> dim[0] - 1) ? loc[1] + 1 : 0)
          {
            sample_shape_t *sampleShape = flss_get_sample_shape(info -> specLabels[loc[0]]);

            size_t dimDerivBlockSub[2] = {sampleShape -> sizeFull, 1};
            mat_t *matDerivBlockSub = mat_new("f", dimDerivBlockSub, false);

            mat_fset_block(matDerivBlock, loc, matDerivBlockSub);
          }

        size_t locDerivT[2] = {locT[0], 0};
        mat_fset_block(matDeriv, locDerivT, matDerivBlock);

        /* Put the single covariance matrices into blocks at each redshift */
        size_t covDimTBlock[2] = {info -> specSize, info -> specSize};
        mat_t *covBlock = mat_new("s", covDimTBlock, _true_);

        for (size_t loc[2] = {0, 0}, index = 0; loc[0] < info -> specSize; loc[0] = (loc[1] < info -> specSize - 1) ? loc[0] : loc[0] + 1, loc[1] = (loc[1] < info -> specSize - 1) ? loc[1] + 1 : loc[0], index++)
          {
            cov_mat_t *covMat = cov_mat_new(info -> covLabels[index]);
            cov_out_add_label(cov_mat_get_out(covMat), cov_in_get_label(info -> covLabels[index]));

            /* Get the block */
            mat_t *covBlockSub = mat_get_block(_covPolyMat_(info -> covLabels[index])(covMat, NULL), locT);

            /* Set the block */
            mat_fset_block(covBlock, loc, covBlockSub);

            /* Free memory */
            covMat = cov_mat_free(covMat);
          }

        /* Reduce the matrix */
        covBlock = mat_reduce(covBlock, &red);

        /* Diagonal matrix */
        mat_t *covBlockDiag = mat_diag(covBlock);

        mat_fset_block(covDiag, locT, covBlockDiag);

        /* Correlation matrix */
        covBlock = mat_corr_ip(covBlock, NULL);

        /* Get the eigenvalue decomposition */
        mats_t *matsBlockEig = mat_eig(covBlock, NULL);

        /* Set the eigenvalues and eigenvectors */
        mat_fset_block(corrEigval, locT, mats_get_mat(matsBlockEig, 0));
        mat_fset_block(corrEigvec, locT, mats_get_mat(matsBlockEig, 1));

        /* Free memory */
        covBlock = mat_free(covBlock);
        matsBlockEig = mats_free(matsBlockEig);
      }

    /* Sizes */
    size_t sizes[5] = {sampleArgZ -> size, sampleArgK -> size, sampleArgMu -> size, mat -> dim[0], mat -> dim[0] * (mat -> dim[1] + 1) / 2};

    /* Print struct */
    print_t *print = print_new();

    print_set_id(print, fishMat -> id);
    print_set_flags(print, fishMat -> flags);

    print_set_sizes(print, sizes);

    /* Store derivatives */
    dat_t **datDeriv = malloc(sizeof(dat_t*) * info -> specSize);

    for (size_t i = 0; i < info -> specSize; i++)
      {
        sample_shape_t *sampleShape = flss_get_sample_shape(info -> specLabels[i]);

        datDeriv[i] = dat_new(0, mat -> dim[0], sampleArgZ -> size * sampleShape -> size);

        for (size_t n = 0; n < datDeriv[i] -> size; n++)
          {
            for (size_t m = 0; m < datDeriv[i] -> yDim; m++)
                dat_set_yvalue(datDeriv[i], m, n, NAN);
          }
      }


    /**  Prepare the Derivatives  **/

    #pragma omp parallel shared(matDeriv)

      { // Start pragma parallel

        /* Declare and define variables */

        /* Spec args struct */
        spec_arg_t *specArg = spec_arg_new(NULL);

        spec_arg_set_dz(specArg, sampleArgZ -> step);
        spec_arg_set_dk(specArg, sampleArgK -> step);
        spec_arg_set_dmu(specArg, sampleArgMu -> step);

        /* Kern struct */
        kern_t *kern = kernels_new_order(info -> specOrders[info -> specSize - 1], kernOrder);

        spec_arg_set_kern(specArg, kern);

        /* Deriv struct */
        spec_deriv_t *deriv = spec_deriv_new();

        spec_arg_set_deriv(specArg, deriv);

        /* Print stage */
        #pragma omp single
          {
            print_fish(print);
          }

        /* Print progress */
        #pragma omp single
          {
            print_fish(print);
          }


        /* Calculate the fisher matrix */

        for (size_t n = 0, f = 0; n < info -> partsSize; n++)

          { // Start info -> partsSize for

            /* Skip if one element does not exist */
            if (!*(info -> partsExist[n]))
                continue;

            /* Set the logarithmic flags */
            spec_deriv_set_log(deriv, *(info -> partsDerivLog[n]));

            for (size_t r = 0; r < info -> partsDerivVals[n] -> size; r++)

              { // Start info -> partsDerivVals -> size for TODO: What about constant parameters?

                /* Set the derivative variables */
                for (size_t i = 0; i < info -> partsDerivVals[n] -> xDim; i++)
                  {
                    spec_deriv_set_var(deriv, dat_get_label(info -> partsDerivVals[n], i, 'x'), dat_get_value(info -> partsDerivVals[n], i, r, 'x'));
                  }

                for (size_t t = 0; t < sampleArgZ -> size; t++)

                  { // Start temporal for

                    /* Skip partial multiplicities with unequal redshifts */
                    if (info -> partsMult[n][0] == -1 && deriv -> z != sampleRawZ -> array[t])
                        continue;

                    /* Set the redshift + fiducials in the kern structs */
                    kernels_set_z(kern, sampleRawZ -> array[t]);

                    /* Temporal block of the eigenvectors */
                    size_t locCorrEigvecT[2] = {t, t};
                    mat_t *corrEigvecT = mat_fget_block(corrEigvec, locCorrEigvecT);

                    /* Temporal block of the diagonal covariance matrix */
                    size_t locCovDiagT[2] = {t, t};
                    mat_t *covDiagT = mat_fget_block(covDiag, locCovDiagT);

                    /* Temporal block of the derivatives times eigenvectors */
                    size_t locDerivT[2] = {t, 0};
                    mat_t *matDerivT = mat_mget_block(matDeriv, locDerivT);

                    for (size_t i1 = 0, i2 = 0; i1 < info -> specSize; i1 = (i2 < info -> specSize - 1) ? i1 : i1 + 1, i2 = (i2 < info -> specSize - 1) ? i2 + 1 : i1)

                      { // Start info -> specSize for

                        /* Variables */
                        sample_shape_t *sampleShape1 = flss_get_sample_shape(info -> specLabels[i1]);
                        sample_shape_t *sampleShape2 = flss_get_sample_shape(info -> specLabels[i2]);

                        sample_raw_t *sampleRawK1 = sampleShape1 -> sampleRawLength;
                        sample_raw_t *sampleRawK2 = sampleShape2 -> sampleRawLength;

                        sample_arg_t *sampleArgK1 = sampleRawK1 -> sampleArg;
                        sample_arg_t *sampleArgK2 = sampleRawK2 -> sampleArg;

                        /* Skip deriv -> k is not included in the sample */
                        if (info -> partsMult[n][1] == 1 && !(misc_bsearch(sampleRawK1 -> array, sampleArgK1 -> size, deriv -> k, __ABSTOL__, NULL)
                                                              && misc_bsearch(sampleRawK2 -> array, sampleArgK2 -> size, deriv -> k, __ABSTOL__, NULL))) continue;

                        /* Spectral derivative function */
                        double (*dSpecFunc)(void*, void*) = _specPoly_(info -> specLabels[i2], info -> partsLabels[n]);

                        /* Average flag */
                        bool avrFlag = flss_get_avr_shape_flag(info -> specLabels[i2]);

                        /* Average function */
                        int (*avrFunc)(double (*)(void*, void*), spec_arg_t *, double*) = (avrFlag) ? avr_shape_get_func(info -> specOrders[i2]) : avr_shape_inf;

                        /* Spectral block of the eigenvectors */
                        size_t locCorrEigvecTBlock[2] = {i1, i2};
                        mat_t *corrEigvecTBlock = mat_fget_block(corrEigvecT, locCorrEigvecTBlock);

                        /* Spectral block of the diagonal covariance matrix */
                        size_t locCovDiagTBlock[2] = {i2, i2};
                        mat_t *covDiagTBlock = mat_fget_block(covDiagT, locCovDiagTBlock);

                        /* Spectral block of the derivatives times eigenvectors */
                        size_t locDerivTBlock[2] = {i2, f};
                        mat_t *matDerivTBlock = mat_mget_block(matDerivT, locDerivTBlock);

                        /* Block does not exist */
                        if (corrEigvecTBlock == NULL)
                            continue;

                        #pragma omp for schedule(dynamic)

                        for (size_t s1 = 0; s1 < sampleShape1 -> size; s1++)

                          { // Start spatial for 1

                            /* (Spatial) Location in the eigenvectors */
                            size_t locCorrEigvecS[2];

                            /* (Spatial) Location in the diagonal covariance matrix */
                            size_t locCovDiagS[2];

                            /* (Spatial) Location in the derivatives times eigenvectors matrix */
                            size_t locDerivS[2] = {0, 0};

                            /* Shape of first spectrum */
                            shape_t *shape1 = sampleShape1 -> arrayShape[s1];

                            for (size_t p1 = 0; p1 < 1 + (size_t) (!shape1 -> parity); p1++)

                              { // Start parity for 1

                                /* (Spatial) Column location in the eigenvectors */
                                locCorrEigvecS[0] = (s1 < sampleShape1 -> sizeParity) ? s1 : sampleShape1 -> sizeParity + 2 * (s1 - sampleShape1 -> sizeParity) + p1;

                                /* (Spatial) Column location in the derivatives times eigenvectors matrix */
                                locDerivS[0] = locCorrEigvecS[0];

                                /* Boundaries */
                                size_t spatialBounds2[2];
                                spatialBounds2[0] = (!strcmp(corrEigvecTBlock -> mtype -> id, "d")) ? s1 : 0;
                                spatialBounds2[1] = (!strcmp(corrEigvecTBlock -> mtype -> id, "d")) ? s1 + 1 : sampleShape2 -> size;

                                /* Scalar product between eigenvector and derivative */
                                double scalarProductDerivEigvec = 0.;

                                for (size_t s2 = spatialBounds2[0]; s2 < spatialBounds2[1]; s2++)

                                  { // Start spatial for 2

                                    /* Shape of second spectrum */
                                    shape_t *shape2 = sampleShape2 -> arrayShape[s2];

                                    /* Set variables for second shape */
                                    for (size_t i = 0; i < shape2 -> dim; i++)
                                      {
                                        kernels_qset_k(kern, i, shape2 -> length[i]);
                                        kernels_qset_mu(kern, i, shape2 -> orientation[i]);

                                        for (size_t j = i + 1; j < shape2 -> dim; j++)
                                            kernels_qset_nu(kern, i, j, shape2 -> angle[shape_get_vertex_angle_index(shape2 -> dim, i, j)]);
                                      }

                                    /* Spectra derivative */
                                    double dSpecVal = dat_get_yvalue(datDeriv[i2], f, t * sampleShape2 -> size + s2);

                                    if (isnan(dSpecVal))
                                      {
                                        dSpecVal = avr_shape_direct(avrFunc, dSpecFunc, specArg);
                                        dat_set_yvalue(datDeriv[i2], f, t * sampleShape2 -> size + s2, dSpecVal);
                                      }

                                    for (size_t p2 = 0; p2 < 1 + (size_t) (!shape2 -> parity); p2++)

                                      { // Start parity for 2

                                        /* (Spatial) Row location in the eigenvectors */
                                        locCorrEigvecS[1] = (s2 < sampleShape2 -> sizeParity) ? s2 : sampleShape2 -> sizeParity + 2 * (s2 - sampleShape2 -> sizeParity) + p2;

                                        /* Spatial element of the correlation matrix */
                                        double valCorrEigvec = mat_get_value(corrEigvecTBlock, locCorrEigvecS);

                                        /* (Spatial) Row location in the diagonal covariance matrix */
                                        locCovDiagS[0] = locCorrEigvecS[1];
                                        locCovDiagS[1] = locCorrEigvecS[1];

                                        /* (Square root of the) Spatial element of the diagonal covariance matrices */
                                        double valCovDiag = sqrt(mat_get_value(covDiagTBlock, locCovDiagS));

                                        /* No need to calculate the spectra */
                                        if (valCorrEigvec == 0.)
                                            continue;

                                        /* Calculate the scalar product between the eigenvector and the derivative */
                                        scalarProductDerivEigvec += valCorrEigvec * dSpecVal / valCovDiag;

                                      } // End parity for 2

                                  } // End spatial for 2

                                /* Insert the scalar product */
                                mat_set_value(matDerivTBlock, locDerivS, scalarProductDerivEigvec);

                              } // End parity for 1

                          } // End spatial for 1

                        /* Free memory */
                        corrEigvecTBlock = mat_ffree(corrEigvecTBlock);
                        covDiagTBlock = mat_ffree(covDiagTBlock);

                      } // End info -> specSize for

                    /* Free memory */
                    corrEigvecT = mat_ffree(corrEigvecT);
                    covDiagT = mat_ffree(covDiagT);

                  } // End temporal for

                #pragma omp single
                  {
                    /* Update the progress */
                    print_update_progress(print, 1. / (double) (print -> sizes[4]));
                    print_fish(print);
                  }

                /* Update the fisher index */
                f++;

              } // End info -> partsDerivVals -> size for

          } // End info -> partsSize for


        /* Stage has finished: print a message and the execution time */

        print -> stageFinished = 1;

        #pragma omp single
          {
            print_set_stagefinished(print, 1);
            print_set_progress(print, 1.);

            print_fish(print);
          }


        /* Free memory */

        specArg = spec_arg_free(specArg);


      } // End pragma parallel

    /* Free memory */
    for (size_t i = 0; i < info -> specSize; i++)
      {
        datDeriv[i] = dat_free(datDeriv[i]);
      }

    free(datDeriv);


    /**  Calculate the Fisher Matrix  **/

    double valFish;

    #pragma omp parallel shared(valFish)

      { // Start pragma parallel

        /* Declare and define variables */

        /* Deriv struct */
        spec_deriv_t *deriv1 = spec_deriv_new();
        spec_deriv_t *deriv2 = spec_deriv_new();

        /* Print stage */
        #pragma omp single
          {
//            print_fish(print);
          }

        /* Print progress */
        #pragma omp single
          {
//            print_fish(print);
          }


        /*  Calculate the Fisher Matrix  */

        for (size_t n = 0, m = 0, f1 = 0, f2 = 0; n < info -> partsSize; n = (m < info -> partsSize - 1) ? n : n + 1, m = (m < info -> partsSize - 1) ? m + 1 : n)

          { // Start info -> partsSize for

            /* Skip if one element does not exist */
            if (!*(info -> partsExist[n]) || !*(info -> partsExist[m]))
                continue;

            /* Set the logarithmic flags */
            spec_deriv_set_log(deriv1, *(info -> partsDerivLog[n]));
            spec_deriv_set_log(deriv2, *(info -> partsDerivLog[m]));

            for (size_t r = 0, s = 0; r < info -> partsDerivVals[n] -> size; r = (s < info -> partsDerivVals[m] -> size - 1) ? r : r + 1, s = (s < info -> partsDerivVals[m] -> size - 1) ? s + 1 : ( (n == m) ? r : 0 ))

              { // Start info -> partsDerivVals -> size for TODO: What about constant parameters?

                /* Set the derivative variables for the first spectrum */
                for (size_t i = 0; i < info -> partsDerivVals[n] -> xDim; i++)
                  {
                    spec_deriv_set_var(deriv1, dat_get_label(info -> partsDerivVals[n], i, 'x'), dat_get_value(info -> partsDerivVals[n], i, r, 'x'));
                  }

                /* Set the derivative variables for the second spectrum */
                for (size_t i = 0; i < info -> partsDerivVals[m] -> xDim; i++)
                  {
                    spec_deriv_set_var(deriv2, dat_get_label(info -> partsDerivVals[m], i, 'x'), dat_get_value(info -> partsDerivVals[m], i, s, 'x'));
                  }

                /* Fisher matrix element */
                #pragma omp single
                valFish = 0.;

                for (size_t t = 0; t < sampleArgZ -> size; t++)

                  { // Start temporal for

                    /* Skip partial multiplicities with unequal redshifts */
                    if ((info -> partsMult[n][0] == -1 && deriv1 -> z != sampleRawZ -> array[t]) || (info -> partsMult[m][0] == -1 && deriv2 -> z != sampleRawZ -> array[t]))
                        continue;

                    /* Temporal block of the eigenvalues */
                    size_t locCorrEigvalT[2] = {t, t};
                    mat_t *corrEigvalT = mat_fget_block(corrEigval, locCorrEigvalT);

                    /* Temporal block of the derivatives */
                    size_t locDerivT[2] = {t, 0};
                    mat_t *matDerivT = mat_fget_block(matDeriv, locDerivT);

                    for (size_t i = 0; i < info -> specSize; i++)

                      { // Start info -> specSize for

                        /* Variables */
                        sample_shape_t *sampleShape = flss_get_sample_shape(info -> specLabels[i]);

                        sample_raw_t *sampleRawK = sampleShape -> sampleRawLength;
                        sample_arg_t *sampleArgK = sampleRawK -> sampleArg;

                        /* Skip if one spectrum's deriv -> k is not included in the sample */
                        if (info -> partsMult[n][1] == 1 && !misc_bsearch(sampleRawK -> array, sampleArgK -> size, deriv1 -> k, __ABSTOL__, NULL))
                            continue;

                        if (info -> partsMult[m][1] == 1 && !misc_bsearch(sampleRawK -> array, sampleArgK -> size, deriv2 -> k, __ABSTOL__, NULL))
                            continue;

                        /* Spectral block of the eigenvalues */
                        size_t locCorrEigvalTBlock[2] = {i, i};
                        mat_t *corrEigvalTBlock = mat_fget_block(corrEigvalT, locCorrEigvalTBlock);

                        /* Spectral blocks of the derivatives */
                        size_t locDerivTBlock1[2] = {i, f1};
                        mat_t *matDerivTBlock1 = mat_fget_block(matDerivT, locDerivTBlock1);

                        size_t locDerivTBlock2[2] = {i, f2};
                        mat_t *matDerivTBlock2 = mat_fget_block(matDerivT, locDerivTBlock2);

                        /* Block does not exist */
                        if (corrEigvalTBlock == NULL)
                            continue;

                        #pragma omp for
                        for (size_t s = 0; s < sampleShape -> size; s++)

                          { // Start spatial for 1

                            shape_t *shape = sampleShape -> arrayShape[s];

                            for (size_t p = 0; p < 1 + (size_t) (!shape -> parity); p++)

                              { // Start parity for 1

                                /* (Spatial) Location in the eigenvalues */
                                size_t locCorrEigvalS[2];

                                locCorrEigvalS[0] = (s < sampleShape -> sizeParity) ? s : sampleShape -> sizeParity + 2 * (s - sampleShape -> sizeParity) + p;
                                locCorrEigvalS[1] = (s < sampleShape -> sizeParity) ? s : sampleShape -> sizeParity + 2 * (s - sampleShape -> sizeParity) + p;

                                /* Value */
                                double eigVal = mat_get_value(corrEigvalTBlock, locCorrEigvalS);

                                /* Skip negative or zero eigenvalues (covariance matrix should be positive definite) */
                                if (eigVal <= 0.)
                                  {
//                                    printf("Skipped the %ld'th eigenvalue (%e).\n", locCorrEigvalS[0], eigVal);

                                    continue;
                                  }

                                /* (Spatial) Location in the derivatives */
                                size_t locDerivS[2];

                                locDerivS[0] = locCorrEigvalS[0];
                                locDerivS[1] = 0;

                                /* Values */
                                double derivVal1 = mat_get_value(matDerivTBlock1, locDerivS);
                                double derivVal2 = mat_get_value(matDerivTBlock2, locDerivS);

                                /* Update the Fisher matrix element */
                                #pragma omp atomic
                                valFish += derivVal1 * derivVal2 / eigVal;

                              } // End parity for 1

                          } // End spatial for 1

                        /* Free memory */
                        corrEigvalTBlock = mat_ffree(corrEigvalTBlock);

                        matDerivTBlock1 = mat_ffree(matDerivTBlock1);
                        matDerivTBlock2 = mat_ffree(matDerivTBlock2);

                      } // End info -> specSize for

                    /* Free memory */
                    corrEigvalT = mat_ffree(corrEigvalT);

                    matDerivT = mat_ffree(matDerivT);

                  } // End temporal for

                #pragma omp single
                  {
                    /* Insert the value */
                    size_t locFish1[2] = {f1, f2};
                    size_t locFish2[2] = {f2, f1};

                    mat_set_value(mat, locFish1, valFish);
                    mat_set_value(mat, locFish2, valFish);

                    /* Update the progress */
//                    print_update_progress(print, 1. / (double) (print -> sizes[4]));
//                    print_fish(print);
                  }

                /* Update the fisher indices */
                if (r < info -> partsDerivVals[n] -> size - 1 || s < info -> partsDerivVals[m] -> size - 1)
                  {
                    f1 = (s < info -> partsDerivVals[m] -> size - 1) ? f1 : f1 + 1;
                    f2 = (s < info -> partsDerivVals[m] -> size - 1) ? f2 + 1 : ( (n == m) ? f1 : f2 + 1 - info -> partsDerivVals[m] -> size );
                  }

              } // End info -> partsDerivVals -> size for

            /* Update the fisher indices */
            f1 = (f2 < mat -> dim[1] - 1) ? f1 + 1 - info -> partsDerivVals[n] -> size : f1 + 1;
            f2 = (f2 < mat -> dim[1] - 1) ? f2 + 1 : f1;

          } // End info -> partsSize for


        /* Stage has finished: print a message and the execution time */

        print -> stageFinished = 1;

        #pragma omp single
          {
            print_set_stagefinished(print, 1);
            print_set_progress(print, 1.);

//            print_fish(print);
          }


        /* Free memory */

        deriv1 = spec_deriv_free(deriv1);
        deriv2 = spec_deriv_free(deriv2);


      } // End pragma parallel


    /* Free memory */

    corrEigval = mat_free(corrEigval);
    corrEigvec = mat_free(corrEigvec);

    matDeriv = mat_free(matDeriv);

    print = print_free(print);


    /* Write the result to file */

    if (fishMat -> out -> file != NULL && fishMat -> out -> precision != 0)
      {
        char *outFile = misc_scat(3, __FISHERLSSDIR__, _outDir, fishMat -> out -> file);

        printf("%s\n", outFile);

        if (fishMat -> out -> binary)
            mats_output(outFile, mats, NULL);

        else
            mats_output(outFile, mats, &fishMat -> out -> precision);

        free(outFile);
      }

    return mats;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


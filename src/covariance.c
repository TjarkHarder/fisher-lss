/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COVARIANCE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "covariance.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE / SETUP / FREE VARIABLES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Local Variables   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Zero function for covariance matrix  **/

static int _cov_zero(cov_arg_t *covArg, double *result)
{
    /*

        Fills result with zeros (always assuming that result has size 2 !!!)

    */

    (void) covArg;

    result[0] = 0.;
    result[1] = 0.;

    return 0;
}


/**  Output directory and file extension  **/

static const char *_outDir = "/output/cov/";
static const char *_outExt = ".mat";


/**  Read in covariance matrix interpolation functions  **/

/* PP Covariance Matrix */

static char *_covppInFile;
static char *_covppInLabel;

static mat_t *_covppInMat;


/* BB Covariance Matrix */

static char *_covbbInFile;
static char *_covbbInLabel;

static mat_t *_covbbInMat;


/* PB Covariance Matrix */

static char *_covpbInFile;
static char *_covpbInLabel;

static mat_t *_covpbInMat;


/**  Information about covariance matrix  **/

typedef struct
{
    /*

        Information about the covariance matrix in condensed form

    */

    const char *id;


    size_t specSize;

    const char **specLabels;
    size_t *specOrders;


    size_t partsSize;

    char **partsLabels;
    bool **partsExist;
    bool **partsHasErr;
    char **partsType;

} _cov_info_t;


/*  ------------------------------------------------------------------------------------------------------  */


static _cov_info_t *_covppInfo = NULL;
static _cov_info_t *_covbbInfo = NULL;
static _cov_info_t *_covpbInfo = NULL;


static _cov_info_t *_cov_info_get_struct(const char *id)
{
    /*

        Get the info struct with given id

    */

    /* PP Covariance Matrix */
    if (!strcmp(id, _idCovPP_))
        return _covppInfo;

    /* BB Covariance Matrix */
    if (!strcmp(id, _idCovBB_))
        return _covbbInfo;

    /* PB Covariance Matrix */
    if (!strcmp(id, _idCovPB_))
        return _covpbInfo;

    printf("Covariance matrix with ID ’%s’ does not exist.\n", id);
    exit(1);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


static _cov_info_t *_cov_info_free(_cov_info_t *info)
{
    /*

        Free a _cov_info_t struct

    */

    if (info == NULL)
        return 0;

    free(info -> specLabels);
    free(info -> specOrders);

    free(info -> partsLabels);
    free(info -> partsExist);
    free(info -> partsHasErr);
    free(info -> partsType);

    free(info);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


size_t cov_info_get_kern_order(const char *id)
{
    /*

        Get the kern order of a covariance matrix

    */

    _cov_info_t *info = _cov_info_get_struct(id);

    size_t kernOrder1 = spec_info_get_kern_order(info -> specLabels[0]) + info -> specOrders[0] % 2;
    size_t kernOrder2 = spec_info_get_kern_order(info -> specLabels[info -> specSize - 1]) + info -> specOrders[info -> specSize - 1] % 2;

    return (kernOrder1 > kernOrder2) ? kernOrder2 : kernOrder1;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  PP Covariance Matrix  **/

/* Part Labels and Types */
static const char *_covppLabelGauss = NULL;
static const char *_covppTypeGauss = "d";

static const char *_covppLabelNGauss = NULL;
static const char *_covppTypeNGauss = "s";

static const char *_covppLabelFull = NULL;
static const char *_covppTypeFull = "s";


/**  BB Covariance Matrix  **/

/* Part Labels and Types */
static const char *_covbbLabelGauss = NULL;
static const char *_covbbTypeGauss = "d";

static const char *_covbbLabelNGauss = NULL;
static const char *_covbbTypeNGauss = "s";

static const char *_covbbLabelFull = NULL;
static const char *_covbbTypeFull = "s";


/* PB Covariance Matrix */

/* Part Labels and Types */
static const char *_covpbLabelGauss = NULL;
static const char *_covpbTypeGauss = "null";

static const char *_covpbLabelNGauss = NULL;
static const char *_covpbTypeNGauss = "f";

static const char *_covpbLabelFull = NULL;
static const char *_covpbTypeFull = "f";



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Initialise Local Variables   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_labels(void);
static int _ini_info(void);
static int _ini_input(void);
static int _ini_func(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_ini(void)
{
    /*


        Initialise covariance matrix modules

    */

    /* Initialise labels */
    _ini_labels();

    /* Initialise info */
    _ini_info();

    /* Initialise input (e.g. covariance matrix from file) */
    _ini_input();

    /* Initialise functions */
    _ini_func();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_labels_covpp();
static int _ini_labels_covbb();
static int _ini_labels_covpb();

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_labels(void)
{
    /*

        Initialise the labels

    */

    /* PP Covariance Matrix */
    _ini_labels_covpp();

    /* BB Covariance Matrix */
    _ini_labels_covbb();

    /* PB Covariance Matrix */
    _ini_labels_covpb();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_labels_covpp(void)
{
    /*

        Initialise the pp covariance matrix labels

    */

    _covppLabelGauss = misc_scat(5, "'", _idCovPP_, "-", _idGauss_, "'");
    _covppLabelNGauss = misc_scat(5, "'", _idCovPP_, "-", _idNGauss_, "'");
    _covppLabelFull = misc_scat(5, "'", _idCovPP_, "-", _idFull_, "'");

    return 0;
}


static int _ini_labels_covbb(void)
{
    /*

        Initialise the bb covariance matrix labels

    */

    _covbbLabelGauss = misc_scat(5, "'", _idCovBB_, "-", _idGauss_, "'");
    _covbbLabelNGauss = misc_scat(5, "'", _idCovBB_, "-", _idNGauss_, "'");
    _covbbLabelFull = misc_scat(5, "'", _idCovBB_, "-", _idFull_, "'");

    return 0;
}


static int _ini_labels_covpb(void)
{
    /*

        Initialise the pb covariance matrix labels

    */

    _covpbLabelGauss = misc_scat(5, "'", _idCovPB_, "-", _idGauss_, "'");
    _covpbLabelNGauss = misc_scat(5, "'", _idCovPB_, "-", _idNGauss_, "'");
    _covpbLabelFull = misc_scat(5, "'", _idCovPB_, "-", _idFull_, "'");

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

    /* PP Covariance Matrix */
    _ini_info_pp();

    /* BB Covariance Matrix */
    _ini_info_bb();

    /* PB Covariance Matrix */
    _ini_info_pb();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_info_pp(void)
{
    /*

        Initialise the pp info struct

    */


    _covppInfo = malloc(sizeof(_cov_info_t));

    _covppInfo -> id = _idCovPP_;


    /**  PP Covariance Matrix Components  **/

    _covppInfo -> specSize = 1;

    _covppInfo -> specLabels = malloc(sizeof(const char*) * _covppInfo -> specSize);
    _covppInfo -> specLabels[0] = _idSpecPnl_;

    _covppInfo -> specOrders = malloc(sizeof(size_t) * _covppInfo -> specSize);
    _covppInfo -> specOrders[0] = spec_info_get_order(_idSpecPnl_);


    /**  Parts of the PP Covariance Matrix  **/

    _covppInfo -> partsSize = 3;

    _covppInfo -> partsLabels = malloc(sizeof(char*) * _covppInfo -> partsSize);
    _covppInfo -> partsHasErr = malloc(sizeof(bool*) * _covppInfo -> partsSize);
    _covppInfo -> partsExist = malloc(sizeof(bool*) * _covppInfo -> partsSize);
    _covppInfo -> partsType = malloc(sizeof(char*) * _covppInfo -> partsSize);

    /* Gaussian */
    _covppInfo -> partsLabels[0] = (char*) _covppLabelGauss;
    _covppInfo -> partsHasErr[0] = (bool*) &_avrShapeLine_;
    _covppInfo -> partsExist[0] = (bool*) &_true_;
    _covppInfo -> partsType[0] = (char*) _covppTypeGauss;

    /* Non-Gaussian */
    _covppInfo -> partsLabels[1] = (char*) _covppLabelNGauss;
    _covppInfo -> partsHasErr[1] = (bool*) &_true_;
    _covppInfo -> partsExist[1] = (bool*) &_true_;
    _covppInfo -> partsType[1] = (char*) _covppTypeNGauss;

    /* Non-linear power spectrum */
    _covppInfo -> partsLabels[2] = (char*) _covppLabelFull;
    _covppInfo -> partsHasErr[2] = (bool*) &_true_;
    _covppInfo -> partsExist[2] = (bool*) &_true_;
    _covppInfo -> partsType[2] = (char*) _covppTypeFull;


    return 0;
}


static int _ini_info_bb(void)
{
    /*

        Initialise the bb info struct

    */


    _covbbInfo = malloc(sizeof(_cov_info_t));

    _covbbInfo -> id = _idCovBB_;


    /**  PP Covariance Matrix Components  **/

    _covbbInfo -> specSize = 1;

    _covbbInfo -> specLabels = malloc(sizeof(const char*) * _covbbInfo -> specSize);
    _covbbInfo -> specLabels[0] = _idSpecBtr_;

    _covbbInfo -> specOrders = malloc(sizeof(size_t) * _covbbInfo -> specSize);
    _covbbInfo -> specOrders[0] = spec_info_get_order(_idSpecBtr_);


    /**  Parts of the PP Covariance Matrix  **/

    _covbbInfo -> partsSize = 3;

    _covbbInfo -> partsLabels = malloc(sizeof(char*) * _covbbInfo -> partsSize);
    _covbbInfo -> partsHasErr = malloc(sizeof(bool*) * _covbbInfo -> partsSize);
    _covbbInfo -> partsExist = malloc(sizeof(bool*) * _covbbInfo -> partsSize);
    _covbbInfo -> partsType = malloc(sizeof(char*) * _covbbInfo -> partsSize);

    /* Gaussian */
    _covbbInfo -> partsLabels[0] = (char*) _covbbLabelGauss;
    _covbbInfo -> partsHasErr[0] = (bool*) &_avrShapeTri_;
    _covbbInfo -> partsExist[0] = (bool*) &_true_;
    _covbbInfo -> partsType[0] = (char*) _covbbTypeGauss;

    /* Non-Gaussian */
    _covbbInfo -> partsLabels[1] = (char*) _covbbLabelNGauss;
    _covbbInfo -> partsHasErr[1] = (bool*) &_avrShapeTri_;
    _covbbInfo -> partsExist[1] = (bool*) &_true_;
    _covbbInfo -> partsType[1] = (char*) _covbbTypeNGauss;

    /* Non-linear power spectrum */
    _covbbInfo -> partsLabels[2] = (char*) _covbbLabelFull;
    _covbbInfo -> partsHasErr[2] = (bool*) &_avrShapeTri_;
    _covbbInfo -> partsExist[2] = (bool*) &_true_;
    _covbbInfo -> partsType[2] = (char*) _covbbTypeFull;


    return 0;
}


static int _ini_info_pb(void)
{
    /*

        Initialise the pb info struct

    */


    _covpbInfo = malloc(sizeof(_cov_info_t));

    _covpbInfo -> id = _idCovPB_;


    /**  PP Covariance Matrix Components  **/

    _covpbInfo -> specSize = 2;

    _covpbInfo -> specLabels = malloc(sizeof(const char*) * _covpbInfo -> specSize);
    _covpbInfo -> specLabels[0] = _idSpecPnl_;
    _covpbInfo -> specLabels[1] = _idSpecBtr_;

    _covpbInfo -> specOrders = malloc(sizeof(size_t) * _covpbInfo -> specSize);
    _covpbInfo -> specOrders[0] = spec_info_get_order(_idSpecPnl_);
    _covpbInfo -> specOrders[1] = spec_info_get_order(_idSpecBtr_);


    /**  Parts of the PP Covariance Matrix  **/

    _covpbInfo -> partsSize = 3;

    _covpbInfo -> partsLabels = malloc(sizeof(char*) * _covpbInfo -> partsSize);
    _covpbInfo -> partsHasErr = malloc(sizeof(bool*) * _covpbInfo -> partsSize);
    _covpbInfo -> partsExist = malloc(sizeof(bool*) * _covpbInfo -> partsSize);
    _covpbInfo -> partsType = malloc(sizeof(char*) * _covpbInfo -> partsSize);

    /* Gaussian */
    _covpbInfo -> partsLabels[0] = (char*) _covpbLabelGauss;
    _covpbInfo -> partsHasErr[0] = (bool*) &_false_;
    _covpbInfo -> partsExist[0] = (bool*) &_true_;
    _covpbInfo -> partsType[0] = (char*) _covpbTypeGauss;

    /* Non-Gaussian */
    _covpbInfo -> partsLabels[1] = (char*) _covpbLabelNGauss;
    _covpbInfo -> partsHasErr[1] = (bool*) &_false_;
    _covpbInfo -> partsExist[1] = (bool*) &_true_;
    _covpbInfo -> partsType[1] = (char*) _covpbTypeNGauss;

    /* Full pb covariance matrix */
    _covpbInfo -> partsLabels[2] = (char*) _covpbLabelFull;
    _covpbInfo -> partsHasErr[2] = (bool*) &_false_;
    _covpbInfo -> partsExist[2] = (bool*) &_true_;
    _covpbInfo -> partsType[2] = (char*) _covpbTypeFull;


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _ini_input_pp(void);
static int _ini_input_bb(void);
static int _ini_input_pb(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_input(void)
{
    /*

        Initialiase variables to input covariance matrices

    */

    /* PP Covariance Matrix */
    _ini_input_pp();

    /* BB Covariance Matrix */
    _ini_input_bb();

    /* PB Covariance Matrix */
    _ini_input_pb();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_input_pp(void)
{
    /*

        Initialise the input variable for the pp covariance matrix

    */

    /* File */
    _covppInFile = misc_scat(2, _covppInfo -> id, _outExt);

    /* Label */
    _covppInLabel = misc_scat(1, _covppLabelGauss);

    /* Matrix */
    _covppInMat = NULL;

    return 0;
}


static int _ini_input_bb(void)
{
    /*

        Initialise the input variable for the bb covariance matrix

    */

    /* File */
    _covbbInFile = misc_scat(2, _covbbInfo -> id, _outExt);

    /* Label */
    _covbbInLabel = misc_scat(1, _covbbLabelGauss);

    /* Matrix */
    _covbbInMat = NULL;

    return 0;
}


static int _ini_input_pb(void)
{
    /*

        Initialise the input variable for the pb covariance matrix

    */

    /* File */
    _covpbInFile = misc_scat(2, _covppInfo -> id, _outExt);

    /* Label */
    _covpbInLabel = misc_scat(1, _covpbLabelGauss);

    /* Matrix */
    _covpbInMat = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _cov_pp(void *var, void *params);
static mat_t *_cov_pp_mat(void *var, void *params);

static double _cov_bb(void *var, void *params);
static mat_t *_cov_bb_mat(void *var, void *params);

static double _cov_pb(void *var, void *params);
static mat_t *_cov_pb_mat(void *var, void *params);

static double (*_cov_poly(const char *covLabel, const char *label))(void*, void*);
static mat_t *(*_cov_poly_mat(const char *covLabel))(void*, void*);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _ini_func(void)
{
    /*

        Initialise covariance matrix functions

    */

    /* PP covariance matrix */
    _covPP_ = _cov_pp;
    _covPPMat_ = _cov_pp_mat;

    /* BB covariance matrix */
    _covBB_ = _cov_bb;
    _covBBMat_ = _cov_bb_mat;

    /* PB covariance matrix */
    _covPB_ = _cov_pb;
    _covPBMat_ = _cov_pb_mat;

    /* Poly covariance matrix */
    _covPoly_ = _cov_poly;
    _covPolyMat_ = _cov_poly_mat;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Free Local Variables   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_labels(void);
static int _free_info(void);
static int _free_input(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_free(void)
{
    /*

        Free local variables

    */

    /* Free labels */
    _free_labels();

    /* Free info */
    _free_info();

    /* Free input */
    _free_input();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_labels_covpp(void);
static int _free_labels_covbb(void);
static int _free_labels_covpb(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_labels(void)
{
    /*

        Free the labels

    */

    /* PP Covariance Matrix */
    _free_labels_covpp();

    /* BB Covariance Matrix */
    _free_labels_covbb();

    /* PB Covariance Matrix */
    _free_labels_covpb();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _free_labels_covpp(void)
{
    /*

        Initialise the pp covariance matrix labels

    */

    /* Free labels */
    free((char*) _covppLabelGauss);
    free((char*) _covppLabelNGauss);
    free((char*) _covppLabelFull);

    /* Reset labels */
    _covppLabelGauss = NULL;
    _covppLabelNGauss = NULL;
    _covppLabelFull = NULL;

    return 0;
}


static int _free_labels_covbb(void)
{
    /*

        Initialise the bb covariance matrix labels

    */

    /* Free labels */
    free((char*) _covbbLabelGauss);
    free((char*) _covbbLabelNGauss);
    free((char*) _covbbLabelFull);

    /* Reset labels */
    _covbbLabelGauss = NULL;
    _covbbLabelNGauss = NULL;
    _covbbLabelFull = NULL;

    return 0;
}


static int _free_labels_covpb(void)
{
    /*

        Initialise the pb covariance matrix labels

    */

    /* Free labels */
    free((char*) _covpbLabelGauss);
    free((char*) _covpbLabelNGauss);
    free((char*) _covpbLabelFull);

    /* Reset labels */
    _covpbLabelGauss = NULL;
    _covpbLabelNGauss = NULL;
    _covpbLabelFull = NULL;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_info(void)
{
    /*

        Free the info structs

    */

    /* PP Covariance Matrix */
    _covppInfo = _cov_info_free(_covppInfo);

    /* BB Covariance Matrix */
    _covbbInfo = _cov_info_free(_covbbInfo);

    /* PB Covariance Matrix */
    _covpbInfo = _cov_info_free(_covpbInfo);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _free_input_covpp(void);
static int _free_input_covbb(void);
static int _free_input_covpb(void);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _free_input(void)
{
    /*

        Free input variables

    */

    /* PP Covariance Matrix */
    _free_input_covpp();

    /* BB Covariance Matrix */
    _free_input_covbb();

    /* PB Covariance Matrix */
    _free_input_covpb();

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _free_input_covpp(void)
{
    /*

        Free input variables for the pp covariance matrix

    */

    free(_covppInFile);
    free(_covppInLabel);

    _covppInMat = mat_free(_covppInMat);

    return 0;
}


static int _free_input_covbb(void)
{
    /*

        Free input variables for the bb covariance matrix

    */

    free(_covbbInFile);
    free(_covbbInLabel);

    _covbbInMat = mat_free(_covbbInMat);

    return 0;
}


static int _free_input_covpb(void)
{
    /*

        Free input variables for the pb covariance matrix

    */

    free(_covpbInFile);
    free(_covpbInLabel);

    _covpbInMat = mat_free(_covpbInMat);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COVARIANCE MATRIX FROM FILE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------------   Setters   ----------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_in_pp_set_file(const char *file)
{
    /*

        Set the input file of the pp covariance matrix

    */

    _covppInFile = realloc(_covppInFile, strlen(file) + 1);
    snprintf(_covppInFile, strlen(file) + 1, "%s", file);

    return 0;
}


int cov_in_pp_set_label(const char *label)
{
    /*

        Set the label of the pp covariance matrix in the input file

    */

    _covppInLabel = realloc(_covppInLabel, strlen(label) + 1);
    snprintf(_covppInLabel, strlen(label) + 1, "%s", label);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_in_bb_set_file(const char *file)
{
    /*

        Set the input file of the bb covariance matrix

    */

    _covbbInFile = realloc(_covbbInFile, strlen(file) + 1);
    snprintf(_covbbInFile, strlen(file) + 1, "%s", file);

    return 0;
}


int cov_in_bb_set_label(const char *label)
{
    /*

        Set the label of the bb covariance matrix in the input file

    */

    _covbbInLabel = realloc(_covbbInLabel, strlen(label) + 1);
    snprintf(_covbbInLabel, strlen(label) + 1, "%s", label);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_in_pb_set_file(const char *file)
{
    /*

        Set the input file of the pb covariance matrix

    */

    _covpbInFile = realloc(_covpbInFile, strlen(file) + 1);
    snprintf(_covpbInFile, strlen(file) + 1, "%s", file);

    return 0;
}


int cov_in_pb_set_label(const char *label)
{
    /*

        Set the label of the pb covariance matrix in the input file

    */

    _covpbInLabel = realloc(_covpbInLabel, strlen(label) + 1);
    snprintf(_covpbInLabel, strlen(label) + 1, "%s", label);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------------   Getters   ----------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *cov_in_get_label(const char *id)
{
    /*

        Get the label of the covariance matrix with given id in the input file

    */

    /* Info struct */
    _cov_info_t *info = _cov_info_get_struct(id);

    if (!strcmp(info -> id, _idCovPP_))
        return _covppInLabel;

    if (!strcmp(info -> id, _idCovBB_))
        return _covbbInLabel;

    if (!strcmp(info -> id, _idCovPB_))
        return _covpbInLabel;

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Setup Covariance Matrix   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_in_pp_setup(size_t *size)
{
    /*

        Setup the pp covariance matrix from file

    */

    /* Path to the file */
    char *path = misc_scat(3, __FISHERLSSDIR__, _outDir, _covppInFile);

    /* Read the matrix */
    _covppInMat = mat_free(_covppInMat);

    mats_t *matsIn = mats_input(path, _covppInLabel, _true_);
    _covppInMat = mats_get_mat(matsIn, 0);

    if (size != NULL)
      {
        /* Dimensions of blocks */
        size_t bdim[2];

        size_t depth = 0;
        size_t *loc = calloc(2, sizeof(size_t));

        mat_get_bdim(_covppInMat, depth, &loc, bdim);

        free(loc);

        /* Generate random numbers for submatrix indices */
        size_t *indices = misc_random_sample(0, bdim[0] - 1, size, false);

        /* Submatrix */
        mat_t *matIn = _covppInMat;
        mat_reduce_t reduce = {1.e-10};

        _covppInMat = mat_new("d", matIn -> dim, true);

        for (size_t i = 0; i < matIn -> dim[2]; i++)
          {
            size_t loc[2] = {i, i};
            mat_t *matInBlock = mat_mget_block(matIn, loc);

            mat_fset_block(_covppInMat, loc, mat_sub(matInBlock, indices, *size, indices, *size, &reduce));
          }

        /* Free memory */
        free(indices);
        mat_free(matIn); // Will not be freed in mats_free call
      }

    /* Free memory */
    free(path);
    matsIn = mats_free(matsIn);

    return 0;
}


int cov_in_pp_reset(void)
{
    /*

        Reset the pp covariance matrix from file

    */

    _covppInMat = mat_free(_covppInMat);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_in_bb_setup(size_t *size)
{
    /*

        Setup the bb covariance matrix from file

    */

    /* Path to the file */
    char *path = misc_scat(3, __FISHERLSSDIR__, _outDir, _covbbInFile);

    /* Read the matrix */
    _covbbInMat = mat_free(_covbbInMat);

    mats_t *matsIn = mats_input(path, _covbbInLabel, _true_);
    _covbbInMat = mats_get_mat(matsIn, 0);

    if (size != NULL)
      {
        /* Dimensions of blocks */
        size_t bdim[2];

        size_t depth = 0;
        size_t *loc = calloc(2, sizeof(size_t));

        mat_get_bdim(_covbbInMat, depth, &loc, bdim);

        free(loc);

        /* Generate random numbers for submatrix indices */
        size_t *indices = misc_random_sample(0, bdim[0] - 1, size, false);

        /* Submatrix */
        mat_t *matIn = _covbbInMat;
        mat_reduce_t reduce = {1.e-10};

        _covbbInMat = mat_new("d", matIn -> dim, true);

        for (size_t i = 0; i < matIn -> dim[2]; i++)
          {
            size_t loc[2] = {i, i};
            mat_t *matInBlock = mat_mget_block(matIn, loc);

            mat_fset_block(_covbbInMat, loc, mat_sub(matInBlock, indices, *size, indices, *size, &reduce));
          }

        /* Free memory */
        free(indices);
        mat_free(matIn); // Will not be freed in mats_free call
      }

    /* Free memory */
    free(path);
    matsIn = mats_free(matsIn);

    return 0;
}


int cov_in_bb_reset(void)
{
    /*

        Reset the bb covariance matrix from file

    */

    _covbbInMat = mat_free(_covbbInMat);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_in_pb_setup(size_t *size1, size_t *size2)
{
    /*

        Setup the pb covariance matrix from file

    */

    /* Path to the file */
    char *path = misc_scat(3, __FISHERLSSDIR__, _outDir, _covpbInFile);

    /* Read the matrix */
    _covpbInMat = mat_free(_covpbInMat);

    mats_t *matsIn = mats_input(path, _covpbInLabel, _true_);
    _covpbInMat = mats_get_mat(matsIn, 0);

    if (size1 != NULL || size2 != NULL)
      {
        /* Dimensions of blocks */
        size_t bdim[2];

        size_t depth = 0;
        size_t *loc = calloc(2, sizeof(size_t));

        mat_get_bdim(_covpbInMat, depth, &loc, bdim);

        free(loc);

        /* Submatrix dimensions */
        size_t size1_ = (size1 != NULL) ? *size1 : bdim[0];
        size_t size2_ = (size2 != NULL) ? *size2 : bdim[1];

        /* generate random numbers for submatrix indices */
        size_t *indices1 = misc_random_sample(0, bdim[0] - 1, &size1_, false);
        size_t *indices2 = misc_random_sample(0, bdim[1] - 1, &size2_, false);

        /* Submatrix */
        mat_t *matIn = _covpbInMat;
        mat_reduce_t reduce = {1.e-10};

        _covpbInMat = mat_new("d", matIn -> dim, true);

        for (size_t i = 0; i < matIn -> dim[2]; i++)
          {
            size_t loc[2] = {i, i};
            mat_t *matInBlock = mat_mget_block(matIn, loc);

            mat_fset_block(_covpbInMat, loc, mat_sub(matInBlock, indices1, *size1, indices2, *size2, &reduce));
          }

        /* Free memory */
        free(indices1);
        free(indices2);
        mat_free(matIn); // Will not be freed in mats_free call
      }

    /* Free memory */
    free(path);
    matsIn = mats_free(matsIn);

    return 0;
}


int cov_in_pb_reset(void)
{
    /*

        Reset the pb covariance matrix from file

    */

    _covpbInMat = mat_free(_covpbInMat);

    return 0;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COVARIANCE MATRIX STRUCTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------   Covariance Matrix Arguments Struct   --------------------------------  */


cov_arg_t *cov_arg_new(const char *id)
{
    /*

        Create a new cov_arg_t struct

    */

    cov_arg_t *covArg = malloc(sizeof(cov_arg_t));

    /* Info */
    _cov_info_t *info = _cov_info_get_struct(id);

    /* Id */
    covArg -> id = info -> id;

    /* specArg */
    covArg -> specArg = spec_arg_new(NULL); // Do not allocate any memory for kern in specArg (ie use NULL as argument)

    /* Kern */
    kern_t *kern = kernels_new_order(info -> specOrders[0] + info -> specOrders[info -> specSize - 1], cov_info_get_kern_order(info -> id));
    spec_arg_set_kern(covArg -> specArg, kern);

    /* Shapes */
    covArg -> shape1 = NULL;
    covArg -> shape2 = NULL;

    return covArg;
}


cov_arg_t *cov_arg_free(cov_arg_t *covArg)
{
    /*

        Free covArg

    */

    /* Check for NULL */
    if (covArg == NULL)
        return NULL;

    /* Free contents */
    covArg -> specArg = spec_arg_free(covArg -> specArg);

    /* Free covArg itself */
    free(covArg);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_arg_set_shape(cov_arg_t *covArg, shape_t *shape, size_t index)
{
    /*

        Set the index'th shape in covArg (can only have 0 or 1)

    */

    if (covArg == NULL)
        return 1;


    if (index == 0)
        covArg -> shape1 = shape;

    else if (index == 1)
        covArg -> shape2 = shape;

    else
      {
        printf("Can only set the shape in a 'cov_arg_t' struct if 'index' is either '0' or '1'. Instead got %ld.\n", index);
        exit(1);

        return 1;
      }

    return 0;
}


int cov_arg_set_dk(cov_arg_t *covArg, double dk, size_t index)
{
    /*

        Set the index'th dk in covArg (can only have 0 or 1)

    */

    if (covArg == NULL)
        return 1;


    if (index == 0)
        covArg -> dk1 = dk;

    else if (index == 1)
        covArg -> dk2 = dk;

    else
      {
        printf("Can only set the k-bin-width in a 'cov_arg_t' struct if 'index' is either '0' or '1'. Instead got %ld.\n", index);
        exit(1);

        return 1;
      }

    return 0;
}


int cov_arg_set_dmu(cov_arg_t *covArg, double dmu, size_t index)
{
    /*

        Set the index'th dmu in covArg (can only have 0 or 1)

    */

    if (covArg == NULL)
        return 1;


    if (index == 0)
        covArg -> dmu1 = dmu;

    else if (index == 1)
        covArg -> dmu2 = dmu;

    else
      {
        printf("Can only set the mu-bin-width in a 'cov_arg_t' struct if 'index' is either '0' or '1'. Instead got %ld.\n", index);
        exit(1);

        return 1;
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


spec_arg_t *cov_arg_get_spec_arg(cov_arg_t *covArg)
{
    /*

        Get the specArg struct from covArg

    */

    return covArg -> specArg;
}


kern_t *cov_arg_get_kern(cov_arg_t *covArg)
{
    /*

        Get the kern struct from covArg

    */

    return covArg -> specArg -> kern;
}


/*  ------------------------------------------------------------------------------------------------------  */





/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   Covariance Matrix Output Struct   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


cov_out_t *cov_out_new(void)
{
    /*

        Create a new cov_out_t struct

    */

    cov_out_t *out = malloc(sizeof(cov_out_t));

    out -> size = 0;
    out -> labels = NULL;

    out -> err = NULL;

    out -> file = NULL;

    out -> binary = _true_;
    out -> precision = 0;

    return out;
}


cov_out_t *cov_out_free(cov_out_t *out)
{
    /*

        Free out

    */

    /* Check for NULL */
    if (out == NULL)
        return NULL;

    /* Free contents */
    for (size_t i = 0; i < out -> size; i++)
      {
        free(out -> labels[i]);
      }

    free(out -> labels);
    free(out -> err);
    free(out -> file);

    /* Free out itself */
    free(out);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


cov_out_t *cov_out_cp(cov_out_t *out)
{
    /*

        Copy out

    */

    cov_out_t *outCp = malloc(sizeof(cov_out_t));

    outCp -> size = out -> size;

    outCp -> labels = malloc(sizeof(char*) * out -> size);
    outCp -> err = malloc(sizeof(bool) * out -> size);

    for (size_t i = 0; i < out -> size; i++)
      {
        outCp -> labels[i] = misc_scat(1, out -> labels[i]);
        outCp -> err[i] = out -> err[i];
      }

    outCp -> file = misc_scat(1, out -> file);

    outCp -> binary = out -> binary;
    outCp -> precision = out -> precision;

    return outCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_out_add_label(cov_out_t *out, const char *label)
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

    /* Add the label */
    out -> labels[out -> size - 1] = misc_scat(1, label);

    /* Set the error flag to _false_ */
    out -> err[out -> size - 1] = _false_;

    return 0;
}


int cov_out_rm_label(cov_out_t *out, const char *label)
{
    /*

        Remove an output label from out

    */

    bool success;

    /* Get the index */
    size_t index = misc_ssearch_index(out -> labels, out -> size, label, &success);

    /* Not in the array */
    if (!success)
      {
        printf("Cannot remove the output label '%s' if it does not exist in the list of output labels.\n", label);

        return 0;
      }

    /* Decrease the size of the array */
    out -> size -= 1;

    /* Push the string to the end of the array */
    for (size_t i = index; i < out -> size; i++)
      {
        misc_swap(&(out -> labels)[i], &(out -> labels)[i+1], "s");
        misc_swap(&(out -> err)[i], &(out -> err)[i+1], "b");
      }

    /* Free the last label */
    free(out -> labels[out -> size]);

    /* Reallocate memory */
    out -> labels = realloc(out -> labels, sizeof(char*) * out -> size);
    out -> err = realloc(out -> err, sizeof(bool*) * out -> size);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int cov_out_set_err(cov_out_t *out, const char *label, bool err)
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
        printf("Cannot set the integration error flag for non-existent label '%s' in a cov_out_t struct!\n", label);

        return 0;
      }

    out -> err[index] = err;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int cov_out_set_file(cov_out_t *out, const char *file)
{
    /*

        Set the file name for out (if file is NULL nothing will be output to file)

    */

    misc_scp(&(out -> file), file);

    return 0;
}


int cov_out_set_binary(cov_out_t *out, bool binary)
{
    /*

        Set the binary flag for out

    */

    out -> binary = binary;

    return 0;
}


int cov_out_set_precision(cov_out_t *out, int precision)
{
    /*

        Set the precision in the output file for out -> also set the binary flag to _false_

    */

    out -> binary = _false_;
    out -> precision = precision;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Covariance Matrix Struct   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


cov_mat_t *cov_mat_new(const char *id)
{
    /*

        New cov_mat_t struct

    */

    cov_mat_t *covMat = malloc(sizeof(cov_mat_t));

    /* Info struct */
    _cov_info_t *info = _cov_info_get_struct(id);

    /* Id */
    covMat -> id = info -> id;

    /* Output */
    covMat -> out = cov_out_new();

    /* Print */
    covMat -> flags = print_flags_new();

    return covMat;
}


/*  ------------------------------------------------------------------------------------------------------  */


cov_mat_t *cov_mat_free(cov_mat_t *covMat)
{
    /*

        Free a cov_mat_t struct

    */

    /* Check for NULL */
    if (covMat == NULL)
        return NULL;


    /* Free contents */

    /* Out */
    covMat -> out = cov_out_free(covMat -> out);

    /* Flags */
    covMat -> flags = print_flags_free(covMat -> flags);


    /* Free covMat itself */
    free(covMat);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


cov_out_t *cov_mat_get_out(cov_mat_t *covMat)
{
    /*

        Get the out struct from covMat

    */

    return covMat -> out;
}


print_flags_t *cov_mat_get_flags(cov_mat_t *covMat)
{
    /*

        Get the print struct from covMat

    */

    return covMat -> flags;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool _cov_mat_test(cov_mat_t *covMat)
{
    /*

        Test if covMat can be used in cov_mat_poly

    */

    /* Info struct */
    _cov_info_t *info = _cov_info_get_struct(covMat -> id);

    /* Spatial variables */
    for (size_t i = 0; i < info -> specSize; i++)
      {
        sample_shape_t *sampleSpatial1 = flss_get_sample_shape(info -> specLabels[i]);
        sample_arg_t *sampleArgK1 = sampleSpatial1 -> sampleRawLength -> sampleArg;
        sample_arg_t *sampleArgMu1 = sampleSpatial1 -> sampleRawOrientation -> sampleArg;

        for (size_t j = i + 1; j < info -> specSize; j++)
          {
            sample_shape_t *sampleSpatial2 = flss_get_sample_shape(info -> specLabels[j]);
            sample_arg_t *sampleArgK2 = sampleSpatial2 -> sampleRawLength -> sampleArg;
            sample_arg_t *sampleArgMu2 = sampleSpatial2 -> sampleRawOrientation -> sampleArg;

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


cov_mat_t *cov_mat_pp_default(void)
{
    /*

        Default parameters for the PP covariance matrix function

    */

    cov_mat_t *covMat = cov_mat_new(_idCovPP_);


    /* Output */

    cov_out_t *out = cov_mat_get_out(covMat);

    /* Labels and errors */
    for (size_t i = 0; i < _covppInfo -> partsSize; i++)
      {
        cov_out_add_label(out, _covppInfo -> partsLabels[i]);
        cov_out_set_err(out, _covppInfo -> partsLabels[i], _covppInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _covppInfo -> id, _outExt);
    cov_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    cov_out_set_precision(out, 16);
    cov_out_set_binary(out, _true_);


    /* Flags */

    print_flags_t *flags = cov_mat_get_flags(covMat);

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);


    return covMat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


cov_mat_t *cov_mat_bb_default(void)
{
    /*

        Default parameters for the BB covariance matrix function

    */

    cov_mat_t *covMat = cov_mat_new(_idCovBB_);


    /* Output */

    cov_out_t *out = cov_mat_get_out(covMat);

    /* Labels and errors */
    for (size_t i = 0; i < _covbbInfo -> partsSize; i++)
      {
        cov_out_add_label(out, _covbbInfo -> partsLabels[i]);
        cov_out_set_err(out, _covbbInfo -> partsLabels[i], _covbbInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _covbbInfo -> id, _outExt);
    cov_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    cov_out_set_precision(out, 16);
    cov_out_set_binary(out, _true_);


    /* Flags */

    print_flags_t *flags = cov_mat_get_flags(covMat);

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);


    return covMat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


cov_mat_t *cov_mat_pb_default(void)
{
    /*

        Default parameters for the PB covariance matrix function

    */

    cov_mat_t *covMat = cov_mat_new(_idCovPB_);


    /* Output */

    cov_out_t *out = cov_mat_get_out(covMat);

    /* Labels and errors */
    for (size_t i = 0; i < _covpbInfo -> partsSize; i++)
      {
        cov_out_add_label(out, _covpbInfo -> partsLabels[i]);
        cov_out_set_err(out, _covpbInfo -> partsLabels[i], _covpbInfo -> partsHasErr[i]);
      }

    /* File */
    char *outFile = misc_scat(2, _covpbInfo -> id, _outExt);
    cov_out_set_file(out, outFile);
    free(outFile);

    /* Precision / Binary */
    cov_out_set_precision(out, 16);
    cov_out_set_binary(out, _true_);


    /* Flags */

    print_flags_t *flags = cov_mat_get_flags(covMat);

    print_flags_set_main(flags, _true_);
    print_flags_set_sub(flags, _false_);
    print_flags_set_fid(flags, _false_);


    return covMat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CALCULATE COVARIANCE MATRICES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Setup Matrix Data Struct   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _setup_mat_append_block(_cov_info_t *info, mats_t *mats, char *label, bool err, char *type);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static mats_t *_cov_setup_mat(cov_mat_t *covMat)
{
    /*

        Setup the output mats_t struct

        TODO: Do not store everything in info... Make it parallelisable

    */

    /* Info struct */
    _cov_info_t *info = _cov_info_get_struct(covMat -> id);

    /* Matrix collection */
    mats_t *mats = mats_new(0);

    /* Variables */
    sample_raw_t *sampleRawZ = flss_get_sample_redshift();

    /* Block matrix dimensions */
    size_t dimT[2] = {sampleRawZ -> sampleArg -> size, sampleRawZ -> sampleArg -> size};

    /* If full matrix is not output, can drop parts if also not output */
    bool fullExist = misc_sins(covMat -> out -> labels, covMat -> out -> size, _idFull_, NULL);

    for (size_t i = 0; i < info -> partsSize; i++)
      {
        bool success;
        misc_ssearch_index(covMat -> out -> labels, covMat -> out -> size, info -> partsLabels[i], &success);

        if (!fullExist && !success)
            info -> partsExist[i] = (bool*) &_false_;

        /* (Non-Null) Cross-Covariance Matrices have errors, only if both spectra are averaged over */
        if (strcmp(info -> specLabels[0], info -> specLabels[info -> specSize - 1]) && strcmp(info -> partsType[i], "null"))
          {
            if (flss_get_avr_shape_flag(info -> specLabels[0]) && flss_get_avr_shape_flag(info -> specLabels[info -> specSize - 1]))
                info -> partsHasErr[i] = (bool*) &_true_;

            else
                info -> partsHasErr[i] = (bool*) &_false_;
          }
      }


    /* Matrices */

    /* Search for the labels and append new matrices */
    for (size_t i = 0; i < info -> partsSize; i++)
      {
        if (!(*(info -> partsExist[i])))
            continue;

        bool success;
        size_t index = misc_ssearch_index(covMat -> out -> labels, covMat -> out -> size, info -> partsLabels[i], &success);

        if (success)
          {
            /* Append another (block) matrix */
            mats_append_empty(mats, "d", dimT, _true_, info -> partsLabels[i]);
            _setup_mat_append_block(info, mats, info -> partsLabels[i], _false_, info -> partsType[i]);

            /* Append another (block) matrix for errors */
            if (*info -> partsHasErr[i] && covMat -> out -> err[index])
              {
                char *errLabel = misc_scat(3, info -> partsLabels[i], " ", _idIntError_);

                mats_append_empty(mats, "d", dimT, _true_, errLabel);
                _setup_mat_append_block(info, mats, info -> partsLabels[i], _true_, info -> partsType[i]);

                free(errLabel);
              }
          }
      }

    /* Did not find any matching label */
    if (mats -> size == 0)
      {
        printf("Labels for the '%s' function were either not provided or did not fit any of the expected labels.\n", info -> id);
      }

    return mats;
}


static int _setup_mat_append_block(_cov_info_t *info, mats_t *mats, char *label, bool err, char *type)
{
    /*

        Setup the blocks of mats at redshift index n

    */

    /* Variables */
    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_shape_t *sampleShape1 = flss_get_sample_shape(info -> specLabels[0]);
    sample_shape_t *sampleShape2 = flss_get_sample_shape(info -> specLabels[info -> specSize - 1]);

    /* Block matrix dimensions */
    size_t dimS[2] = {sampleShape1 -> sizeFull, sampleShape2 -> sizeFull};

    /* Last matrix */
    mat_t *mat = mats_get_mat(mats, mats -> size - 1);

    /* Add the blocks, their labels and their type */
    for (size_t i = 0; i < sampleRawZ -> sampleArg -> size; i++)
      {
        /* Location */
        size_t loc[2] = {i, i};

        /* Get the z value as a string */
        char *zToString = misc_dtos(sampleRawZ -> array[i], __XSBUFFER__);

        /* Get the label */
        char *totLabel = (err) ? misc_scat(7, label, " [", _idVarZ_, " = ", zToString, "] ", _idIntError_) : misc_scat(6, label, " [", _idVarZ_, " = ", zToString, "]");

        /* Allocate memory */
        mat_t *matBlock = mat_new(type, dimS, false);
        mat_set_label(matBlock, totLabel);

        /* Insert the block */
        mat_fset_block(mat, loc, matBlock);

        /* Free memory */
        free(zToString);
        free(totLabel);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_cov_get_mat(mats_t *mats, const char *label)
{
    /*

        Get the covariance matrix with some label

    */

    return mats_get_mat_by_label(mats, label);
}


static mat_t *_cov_get_mat_err(mats_t *mats, const char *label)
{
    /*

        Get the integration error covariance matrix with some label

    */

    /* Label of the matrix */
    char *labelErr = misc_scat(3, label, " ", _idIntError_);

    /* Get the matrix */
    mat_t *mat = mats_get_mat_by_label(mats, labelErr);

    /* Free memory */
    free(labelErr);

    return mat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------   Covariance Matrix Function   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *cov_mat_poly(cov_mat_t *covMat)
{
    /*

        Calculate the covariance matrix

    */

    /* Exit if nothing should be calculated */

    /* covMat itself is NULL */
    if (covMat == NULL)
      {
        return NULL;
      }

    /* Outsize is 0 */
    if (covMat -> out -> size == 0)
      {
        return NULL;
      }

    /* covMat is faulty */
    if (!_cov_mat_test(covMat))
      {
        return NULL;
      }


    // TODO: Clean this mess up...

    /* Get parameters from covMat */

    /* Info struct */
    _cov_info_t *info = _cov_info_get_struct(covMat -> id);

    /* Get the output matrices struct */
    mats_t *mats = _cov_setup_mat(covMat);

    /* Variables */
    sample_shape_t *sampleShape1 = flss_get_sample_shape(info -> specLabels[0]);
    sample_shape_t *sampleShape2 = flss_get_sample_shape(info -> specLabels[info -> specSize - 1]);

    sample_raw_t *sampleRawZ = flss_get_sample_redshift();
    sample_arg_t *sampleArgZ = sampleRawZ -> sampleArg;

    sample_raw_t *sampleRawK1 = (sample_raw_t*) sampleShape1 -> sampleRawLength;
    sample_arg_t *sampleArgK1 = sampleRawK1 -> sampleArg;

    sample_raw_t *sampleRawK2 = (sample_raw_t*) sampleShape2 -> sampleRawLength;
    sample_arg_t *sampleArgK2 = sampleRawK2 -> sampleArg;

    sample_raw_t *sampleRawMu1 = (sample_raw_t*) sampleShape1 -> sampleRawOrientation;
    sample_arg_t *sampleArgMu1 = sampleRawMu1 -> sampleArg;

    sample_raw_t *sampleRawMu2 = (sample_raw_t*) sampleShape2 -> sampleRawOrientation;
    sample_arg_t *sampleArgMu2 = sampleRawMu2 -> sampleArg;

    size_t sizes[6] = {sampleArgZ -> size, sampleArgK1 -> size, sampleArgMu1 -> size, sampleArgK2 -> size, sampleArgMu2 -> size, 0};

    for (size_t i = 0; i < info -> partsSize; i++)
      {
        if (!(*info -> partsExist[i]) || (misc_sin(info -> partsLabels[i], _idFull_) && info -> partsSize != 1))
            continue;

        if (!strcmp(info -> partsType[i], "d"))
            sizes[5] += sampleArgZ -> size * sampleShape1 -> size;

        else if (!strcmp(info -> partsType[i], "s"))
            sizes[5] += sampleArgZ -> size * sampleShape1 -> size * (sampleShape1 -> size + 1) / 2;

        else
            sizes[5] += sampleArgZ -> size * sampleShape1 -> size * sampleShape2 -> size;
      }

    /* Get the full covariance matrix */
    mat_t *matFull = NULL;
    mat_t *matFullErr = NULL;

    for (size_t i = 0; i < info -> partsSize; i++)
      {
        if (misc_sin(info -> partsLabels[i], _idFull_))
          {
            matFull = _cov_get_mat(mats, info -> partsLabels[i]);
            matFullErr = _cov_get_mat_err(mats, info -> partsLabels[i]);
          }
      }

    /* Print struct */
    print_t *print = print_new();

    print_set_id(print, covMat -> id);
    print_set_flags(print, covMat -> flags);

    print_set_sizes(print, sizes);


        /* Integrate */

    /* Parallel threading */

    #pragma omp parallel

      { // Start pragma parallel

        /* Declare and define variables */

        /* Covariance matrix arguments struct */
        cov_arg_t *covArg = cov_arg_new(covMat -> id);

        covArg -> dk1 = sampleArgK1 -> step;
        covArg -> dmu1 = sampleArgMu1 -> step;

        covArg -> dk2 = sampleArgK2 -> step;
        covArg -> dmu2 = sampleArgMu2 -> step;

        /* specArg struct */
        spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

        spec_arg_set_dz(specArg, sampleArgZ -> step);

        /* Kern struct */
        kern_t *kern = cov_arg_get_kern(covArg);

        /* Print stage */
        #pragma omp single
          {
            print_cov(print);
          }


        /* Calculate the Covariance Matrix */

        for (size_t i = 0; i < info -> partsSize; i++)

          { // Start info -> partsSize for

            /* Skip matrices that don't exist */
            if (!(*(info -> partsExist[i])))
                continue;

            /* Full Matrix is calculates in parts */
            if (misc_sin(info -> partsLabels[i], _idFull_))
                continue;

            /* Current Part's Matrix and Error Matrix */
            mat_t *mat = _cov_get_mat(mats, info -> partsLabels[i]);
            mat_t *matErr = _cov_get_mat_err(mats, info -> partsLabels[i]);

            for (size_t t = 0; t < sampleArgZ -> size; t++)

              { // Start temporal for

                /* (Temporal) Location in the matrix */
                size_t locT[2] = {t, t};

                /* Redshift + Fiducials */
                kernels_set_z(kern, sampleRawZ -> array[t]);

                /* Print fiducials */
                #pragma omp single
                  {
                    print_set_kern(print, kern);

                    print_cov(print);
                  }

                /* Temporal Blocks of the Current Part's Matrices */
                mat_t *matT = (mat == NULL) ? NULL : mat_mget_block(mat, locT);
                mat_t *matTErr = (matErr == NULL) ? NULL : mat_mget_block(matErr, locT);

                /* Temporal Blocks of the Full Matrices */
                mat_t *matTFull = (matFull == NULL) ? NULL : mat_mget_block(matFull, locT);
                mat_t *matTFullErr = (matFullErr == NULL) ? NULL : mat_mget_block(matFullErr, locT);


                /* Calculate the covariance matrix for the current redshift */

                /* Bounds */
                size_t spatialBounds1[2];
                spatialBounds1[0] = 0;
                spatialBounds1[1] = (!strcmp(info -> partsType[i], "null")) ? 0 : sampleShape1 -> size;

                #pragma omp for schedule(dynamic)

                for (size_t s1 = spatialBounds1[0]; s1 < spatialBounds1[1]; s1++)

                  { // Start spatial for 1

                    /* (Spatial) Location in the matrix */
                    size_t locS[2], locSP[2];

                    /* Shape of first spectrum */
                    shape_t *shape1 = shape_cp(sampleShape1 -> arrayShape[s1]);

                    /* Set the first shape */
                    covArg -> i1 = s1;
                    covArg -> shape1 = shape1;

                    /* (Spatial) Row dimension in the matrix */
                    locS[0] = (s1 < sampleShape1 -> sizeParity) ? s1 : sampleShape1 -> sizeParity + 2 * (s1 - sampleShape1 -> sizeParity) + 0; // location of shape1
                    locSP[0] = (s1 < sampleShape1 -> sizeParity) ? s1 : sampleShape1 -> sizeParity + 2 * (s1 - sampleShape1 -> sizeParity) + 1; // location of parity transformed shape1

                    /* Set variables of first shape */
                    for (size_t i = 0; i < shape1 -> dim; i++)
                      {
                        kernels_qset_k(kern, i, shape_get_vertex_length(shape1, i));
                        kernels_qset_mu(kern, i, shape_get_vertex_orientation(shape1, i));

                        for (size_t j = i + 1; j < shape1 -> dim; j++)
                            kernels_qset_nu(kern, i, j, shape_get_vertex_angle(shape1, i, j));
                      }

                    /* Bounds */
                    size_t spatialBounds2[2];
                    spatialBounds2[0] = (!strcmp(info -> partsType[i], "d") || !strcmp(info -> partsType[i], "s")) ? s1 : 0;
                    spatialBounds2[1] = (!strcmp(info -> partsType[i], "d")) ? s1 + 1 : sampleShape2 -> size;

//                    #pragma omp for schedule(dynamic)

                    for (size_t s2 = spatialBounds2[0]; s2 < spatialBounds2[1]; s2++)

                      { // Start spatial for 2

                        /* Shape of second spectrum */
                        shape_t *shape2 = shape_cp(sampleShape2 -> arrayShape[s2]);

                        /* Set the second shape */
                        covArg -> i2 = s2;
                        covArg -> shape2 = shape2;

                        for (int p2 = 0; p2 < 1 + (!shape2 -> parity); p2++)

                          { // Start parity for 2

                            /* Set variables for second shape */
                            for (size_t i = 0; i < shape2 -> dim; i++)
                              {
                                kernels_qset_k(kern, i + shape1 -> dim, shape_get_vertex_length(shape2, i));
                                kernels_qset_mu(kern, i + shape1 -> dim, shape_get_vertex_orientation(shape2, i));

                                for (size_t j = i + 1; j < shape2 -> dim; j++)
                                    kernels_qset_nu(kern, i + shape1 -> dim, j + shape1 -> dim, shape_get_vertex_angle(shape2, i, j));
                              }

                            /* (Spatial) Column dimension in the matrix */
                            locS[1] = (s2 < sampleShape2 -> sizeParity) ? s2 : sampleShape2 -> sizeParity + 2 * (s2 - sampleShape2 -> sizeParity) + (size_t) p2; // location of shape2
                            locSP[1] = (s2 < sampleShape2 -> sizeParity) ? s2 : sampleShape2 -> sizeParity + 2 * (s2 - sampleShape2 -> sizeParity) + (size_t) (1 + p2) % 2; // location of parity transformed shape2

                            /* Print what every thread is calculating. */
                            #pragma omp critical
                              {
                                print_set_thread(print, omp_get_thread_num());
                                print_set_kern(print, kern);

                                print_set_subfinished(print, 0);

                                print_cov(print);
                              }


                            /* Calculate the covariance matrix contributions */

                            /* Covariance matrix results */
                            double result[3];
                            (cov_poly(covMat -> id, info -> partsLabels[i]))(covArg, result);

                            /* Insert the result */
                            if (matT != NULL)
                              {
                                mat_set_value(matT, locS, result[0]);
                                mat_set_value(matT, locSP, result[0]);

                                /* Insert the integral error */
                                if (matTErr != NULL)
                                  {
                                    mat_set_value(matTErr, locS, (result[0] == 0.) ? result[1] : result[1] / result[0]);
                                    mat_set_value(matTErr, locSP, (result[0] == 0.) ? result[1] : result[1] / result[0]);
                                  }
                              }

                            /* Add the result to the full matrix */
                            if (matTFull != NULL)
                              {
                                /* Get the previous value */
                                double val = mat_get_value(matTFull, locS);

                                /* Add the current result */
                                val += result[0];

                                /* Insert the value */
                                mat_set_value(matTFull, locS, val);
                                mat_set_value(matTFull, locSP, val);

                                /* Add the integral error to the full matrix */
                                if (matTFullErr != NULL)
                                  {
                                    /* Get the previous value */
                                    double valErr = mat_get_value(matTFullErr, locS);

                                    /* Add the current error */
                                    valErr = sqrt(valErr*valErr * (val-result[0])*(val-result[0]) + result[1]*result[1]);
                                    valErr = (val == 0.) ? valErr : valErr / val;

                                    /* Insert the value */
                                    mat_set_value(matTFullErr, locS, valErr);
                                    mat_set_value(matTFullErr, locSP, valErr);
                                  }
                              }

                            /* Parity transform second shape */
                            shape_parity(shape2);

                          } // End parity for 2

                        /* Thread has finished */
                        #pragma omp critical
                          {
                            print_set_subfinished(print, 1);
                            print_update_progress(print, 1. / ((double) print -> sizes[5]));

                            print_cov(print);
                          }

                        /* Free memory */
                        shape2 = shape_free(shape2);

                      } // End spatial for 2

                    /* Free memory */
                    shape1 = shape_free(shape1);

                  } // End spatial for 1

              } // End temporal for

          } // End info -> partsSize for


        /* Stage has finished */

        #pragma omp single
          {
            print_set_stagefinished(print, 1);
            print_set_progress(print, 1.);

            print_cov(print);
          }


        /* Free memory */

        covArg = cov_arg_free(covArg);

      } // End pragma parallel


    /* Free memory */

    print = print_free(print);


    /* Write the result to file */

    if (covMat -> out -> file != NULL)
      {
        char *outFile = misc_scat(3, __FISHERLSSDIR__, _outDir, covMat -> out -> file);

        printf("%s\n", outFile);

        if (covMat -> out -> binary)
            mats_output(outFile, mats, NULL);

        else if (covMat -> out -> precision != 0)
            mats_output(outFile, mats, &covMat -> out -> precision);

        free(outFile);
      }


    return mats;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PP COVARIANCE MATRIX     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   (Direct) PP Covariance Matrix   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_pp_gauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the gaussian part of the PP covariance matrix element

    */

    /* specArg struct */
    spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

    /* Kern struct */
    kern_t *kern = spec_arg_get_kern(specArg);

    /* Only calculate diagonal terms */
    if (shape_comp(covArg -> shape1, covArg -> shape2) == 0)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Additional factor */
    double factor = 1. / (1.e9 * shape_get_volume(covArg -> shape1) * kern -> surv -> V);

    /* Result */
    if (_avrShapeLine_)
        avr_cov_pp_gauss(covArg, result);

    else
        avr_cov_pp_gauss_inf(covArg, result);

    result[0] *= factor;
    result[1] *= factor;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int cov_pp_ngauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the non-gaussian part of the PP covariance matrix element

    */

    /* specArg struct */
    spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

    /* Kern struct */
    kern_t *kern = spec_arg_get_kern(specArg);

    /* Additional factor */
    double factor = 1. / (kern -> surv -> V * 1.e9);

    /* Result */
    if (_avrShapeLine_)
        avr_cov_pp_ngauss(covArg, result);

    else
        avr_cov_pp_ngauss_inf(covArg, result);

    result[0] *= factor;
    result[1] *= factor;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int cov_pp_full(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the full PP covariance matrix element

    */

    /* Gaussian pp covariance matrix */
    double covppGauss[3];
    cov_pp_gauss(covArg, covppGauss);

    /* Non-Gaussian pp covariance matrix */
    double covppNGauss[3];
    cov_pp_ngauss(covArg, covppNGauss);

    /* Final result */
    result[0] = covppGauss[0] + covppNGauss[0];
    result[1] = sqrt(covppGauss[1]*covppGauss[1] + covppNGauss[1]*covppNGauss[1]);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*cov_pp(const char *label))(cov_arg_t*, double*)
{
    /*

        Return the requested function that calculates part (or full) of the pp covariance matrix

    */

    /* Gaussian PP Covariance Matrix */
    if (!strcmp(label, _covppLabelGauss))
      {
        return cov_pp_gauss;
      }

    /* Non-Gaussian PP Covariance Matrix */
    if (!strcmp(label, _covppLabelNGauss))
      {
        return cov_pp_ngauss;
      }

    /* Full PP Covariance Matrix */
    if (!strcmp(label, _covppLabelFull))
      {
        return cov_pp_full;
      }

    /* Label not found */
    return _cov_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   (Indirect) PP Covariance Matrix   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _cov_pp(void *var, void *params)
{
    /*

        Calculate the pp covariance matrix

    */

    (void) params;

//    if (_covppInMat == NULL)
      {
        double result[3];
        cov_pp_full(var, result);

        return result[0];
      }

    /* TODO : Must get sample args from file as well... (to know which index the ki, mui belong to in the matrix) */

    printf("Not yet supported!\n");
    exit(1);

    return NAN;
}


static mat_t *_cov_pp_mat(void *var, void *params)
{
    /*

        TODO: Get the full PP covariance matrix

    */

    (void) params;

    cov_mat_t *covMat = (cov_mat_t*) var;

    /* Must first get the covariance matrix */
    if (_covppInMat == NULL || (_covppInMat -> dim[0] != _sampleRedshift_ -> sampleArg -> size || _covppInMat -> mblock -> bdim[0][0] != _sampleShapeLine_ -> sizeFull))
      {
        mats_t *cov = cov_mat_poly(covMat);

        _covppInMat = mat_free(_covppInMat);
        _covppInMat = mat_cp(mats_get_mat(cov, 0));

        cov = mats_free_full(cov);
      }

    return _covppInMat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     BB COVARIANCE MATRIX     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   (Direct) Tree-Level Trispectrum   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_bb_gauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the gaussian part of the BB covariance matrix element

    */

    /* specArg struct */
    spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

    /* Kern struct */
    kern_t *kern = spec_arg_get_kern(specArg);

    /* Variables */
    shape_t *shape1 = covArg -> shape1;
    shape_t *shape2 = covArg -> shape2;

    /* Triangles must be the same */
    int comp = shape_comp(shape1, shape2);

    if (comp == -1 || comp == 0)
      {
        result[0] = 0.;
        result[1] = 0.;

        return 0;
      }

    /* Get the degeneracies */
    size_t *degen = shape_degen_vertexes(shape1);

    /* Additional factor */
    double factor = ((degen[0] == 1) ? 1.
                                     : ((degen[0] == 2) ? 2.
                                                        : 6.)) / (1.e9 * shape_get_volume(shape1) * kern -> surv -> V);

    /* Result */
    if (_avrShapeTri_)
        avr_cov_bb_gauss(covArg, result);

    else
        avr_cov_bb_gauss_inf(covArg, result);

    result[0] *= factor;
    result[1] *= factor;

    /* Reset variables */
    kernels_qset_k(kern, 0, shape_get_vertex_length(shape1, 0));
    kernels_qset_k(kern, 1, shape_get_vertex_length(shape1, 1));

    kernels_qset_mu(kern, 0, shape_get_vertex_orientation(shape1, 0));
    kernels_qset_mu(kern, 1, shape_get_vertex_orientation(shape1, 1));

    kernels_qset_nu(kern, 0, 1, shape_get_vertex_angle(shape1, 0, 1));

    /* Free memory */
    free(degen);

    return 0;
}


int cov_bb_ngauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the non-gaussian part of the BB covariance matrix element

    */

    /* specArg struct */
    spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

    /* Kern struct */
    kern_t *kern = spec_arg_get_kern(specArg);

    /* Variables */
    size_t dim = 3;
    shape_t *shape1 = covArg -> shape1;
    shape_t *shape2 = covArg -> shape2;

    /* Get the degeneracies */
    size_t *degen1 = shape_degen_vertexes(shape1);
    size_t *degen2 = shape_degen_vertexes(shape2);

    /* Initialise the result */
    result[0] = 0.;
    result[1] = 0.;

    /* Additional factor */
    double factor = 1. / (1.e9 * kern -> surv -> V);


    /**  Trispectrum x Power Spectrum Contributions  **/

    double resultTP[3];

    for (size_t i = 0; i < dim; i++)

      { // Start loop over shape1

        /* Skip vertexes that have already been included */
        if (degen1[i] == 0)
            continue;

        for (size_t j = 0; j < dim; j++)

          { // Start loop over shape2

            /* Skip vertexes that have already been included */
            if (degen2[j] == 0)
                continue;

            /* Edges don't match */
            if (_avrShapeTri_)
              {
                if (!shape_comp_vertex_bin(shape1, i, shape2, j, (covArg -> dk1 + covArg -> dk2) / 2., (covArg -> dmu1 + covArg -> dmu2) / 2.))
                    continue;
              }

            else
              {
                if (!shape_comp_vertex(shape1, i, shape2, j))
                    continue;
              }

            /* Set the variables in the correct order */
            for (size_t n1 = 0; n1 < dim; n1++)
              {
                /* Permutate the first triangle */
                kernels_qset_k(kern, n1, shape_get_vertex_length(shape1, (i + n1) % dim));
                kernels_qset_mu(kern, n1, shape_get_vertex_orientation(shape1, (i + n1) % dim));

                for (size_t n2 = n1 + 1; n2 < dim; n2++)
                    kernels_qset_nu(kern, n1, n2, shape_get_vertex_angle(shape1, (i + n1) % dim, (i + n2) % dim));

                /* Permutate the second triangle */
                kernels_qset_k(kern, n1 + dim, shape_get_vertex_length(shape2, (j + n1) % dim));
                kernels_qset_mu(kern, n1 + dim, shape_get_vertex_orientation(shape2, (j + n1) % dim));

                for (size_t n2 = n1 + 1; n2 < dim; n2++)
                    kernels_qset_nu(kern, n1 + dim, n2 + dim, shape_get_vertex_angle(shape2, (j + n1) % dim, (j + n2) % dim));
              }

            /* Calculate the result */
            if (_avrShapeTri_)
                avr_cov_bb_ngauss_tp(covArg, resultTP);

            else
                avr_cov_bb_ngauss_tp_inf(covArg, resultTP);

            resultTP[0] *= (double) (degen1[i] * degen2[j]) * factor;
            resultTP[1] *= (double) (degen1[i] * degen2[j]) * factor;

            /* Increment the total result */
            result[0] += resultTP[0];
            result[1] = sqrt(result[1]*result[1] + resultTP[1]*resultTP[1]);

          } // End loop over shape2

      } // End loop over shape1


    /**  Bispectrum x Bispectrum Contributions  **/

    double resultBB[3];

    for (size_t i = 0; i < dim; i++)

      { // Start loop over shape1

        /* Skip vertexes that have already been included */
        if (degen1[i] == 0)
            continue;

        for (size_t j = 0; j < dim; j++)

          { // Start loop over shape2

            /* Skip vertexes that have already been included */
            if (degen2[j] == 0)
                continue;

            /* Edges don't match */
            if (_avrShapeTri_)
              {
                if (!shape_comp_vertex_parity_bin(shape1, i, shape2, j, (covArg -> dk1 + covArg -> dk2) / 2., (covArg -> dmu1 + covArg -> dmu2) / 2.))
                    continue;
              }

            else
              {
                if (!shape_comp_vertex_parity(shape1, i, shape2, j))
                    continue;
              }

            /* Set the variables in the correct order */
            for (size_t n1 = 0; n1 < dim; n1++)
              {
                /* Permutate the first triangle */
                kernels_qset_k(kern, n1, shape_get_vertex_length(shape1, (i + n1) % dim));
                kernels_qset_mu(kern, n1, shape_get_vertex_orientation(shape1, (i + n1) % dim));

                for (size_t n2 = n1 + 1; n2 < dim; n2++)
                    kernels_qset_nu(kern, n1, n2, shape_get_vertex_angle(shape1, (i + n1) % dim, (i + n2) % dim));

                /* Permutate the second triangle */
                kernels_qset_k(kern, n1 + dim, shape_get_vertex_length(shape2, (j + n1) % dim));
                kernels_qset_mu(kern, n1 + dim, shape_get_vertex_orientation(shape2, (j + n1) % dim));

                for (size_t n2 = n1 + 1; n2 < dim; n2++)
                    kernels_qset_nu(kern, n1 + dim, n2 + dim, shape_get_vertex_angle(shape2, (j + n1) % dim, (j + n2) % dim));
              }

            /* Calculate the result */
            if (_avrShapeTri_)
                avr_cov_bb_ngauss_bb(covArg, resultBB);

            else
                avr_cov_bb_ngauss_bb_inf(covArg, resultBB);

            resultBB[0] *= (double) (degen1[i] * degen2[j]) * factor;
            resultBB[1] *= (double) (degen1[i] * degen2[j]) * factor;

            /* Increment the total result */
            result[0] += resultBB[0];
            result[1] = sqrt(result[1]*result[1] + resultBB[1]*resultBB[1]);

          } // End loop over shape2

      } // End loop over shape1

    /* Reset variables */
    for (size_t i = 0; i < dim; i++)
      {
        kernels_qset_k(kern, i, shape_get_vertex_length(shape1, i));
        kernels_qset_k(kern, i + dim, shape_get_vertex_length(shape2, i));

        kernels_qset_mu(kern, i, shape_get_vertex_orientation(shape1, i));
        kernels_qset_mu(kern, i + dim, shape_get_vertex_orientation(shape2, i));

        for (size_t j = i + 1; j < dim; j++)
          {
            kernels_qset_nu(kern, i, j, shape_get_vertex_angle(shape1, i, j));
            kernels_qset_nu(kern, i + dim, j + dim, shape_get_vertex_angle(shape2, i, j));
          }
      }

    /* Free memory */
    free(degen1);
    free(degen2);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int cov_bb_full(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the full BB covariance matrix element

    */

    /* Gaussian pp covariance matrix */
    double covbbGauss[3];
    cov_bb_gauss(covArg, covbbGauss);

    /* Non-Gaussian bb covariance matrix */
    double covbbNGauss[3];
    cov_bb_ngauss(covArg, covbbNGauss);

    /* Final result */
    result[0] = covbbGauss[0] + covbbNGauss[0];
    result[1] = sqrt(covbbGauss[1]*covbbGauss[1] + covbbNGauss[1]*covbbNGauss[1]);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*cov_bb(const char *label))(cov_arg_t*, double*)
{
    /*

        Return the requested function that calculates part (or full) of the bb covariance matrix

    */

    /* Gaussian PP Covariance Matrix */
    if (!strcmp(label, _covbbLabelGauss))
      {
        return cov_bb_gauss;
      }

    /* Non-Gaussian PP Covariance Matrix */
    if (!strcmp(label, _covbbLabelNGauss))
      {
        return cov_bb_ngauss;
      }

    /* Full PP Covariance Matrix */
    if (!strcmp(label, _covbbLabelFull))
      {
        return cov_bb_full;
      }

    /* Label not found */
    return _cov_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   (Indirect) BB Covariance Matrix   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _cov_bb(void *var, void *params)
{
    /*

        Calculate the bb covariance matrix

    */

    (void) params;

//    if (_covbbInMat == NULL)
      {
        double result[3];
        cov_bb_full(var, result);

        return result[0];
      }

    /* TODO : Must get sample args from file as well... (to know which index the ki, mui belong to in the matrix) */

    printf("Not yet supported!\n");
    exit(1);

    return NAN;
}


static mat_t *_cov_bb_mat(void *var, void *params)
{
    /*

        TODO: Get the full BB covariance matrix

    */

    (void) params;

    cov_mat_t *covMat = (cov_mat_t*) var;

    /* Must first get the covariance matrix */
    if (_covbbInMat == NULL || (_covbbInMat -> dim[0] != _sampleRedshift_ -> sampleArg -> size || _covbbInMat -> mblock -> bdim[0][0] != _sampleShapeTri_ -> sizeFull))
      {
       mats_t *cov = cov_mat_poly(covMat);

        _covbbInMat = mat_free(_covbbInMat);
        _covbbInMat = mat_cp(mats_get_mat(cov, 0));

        cov = mats_free_full(cov);
      }

    return _covbbInMat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PB COVARIANCE MATRIX     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   (Direct) PB Covariance Matrix   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int cov_pb_gauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the gaussian part of the PB covariance matrix element

    */


    (void) covArg;
    (void) result;

    result[0] = 0.;
    result[1] = 0.;

    return 0;
}


int cov_pb_ngauss(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the non-gaussian part of the PB covariance matrix element

    */

    /* specArg struct */
    spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

    /* Kern struct */
    kern_t *kern = spec_arg_get_kern(specArg);

    /* Variables */
    shape_t *shape1 = covArg -> shape1;
    shape_t *shape2 = covArg -> shape2;

    /* Get the degeneracies */
    size_t *degen1 = shape_degen_vertexes(shape1);
    size_t *degen2 = shape_degen_vertexes(shape2);

    /* Initialise the result */
    result[0] = 0.;
    result[1] = 0.;

    /* Additional factor */
    double factor = 1. / (1.e9 * kern -> surv -> V);


    /**  Bispectrum x Power Spectrum Contributions  **/

    double resultBP[3];

    for (size_t i = 0; i < shape1 -> dim; i++)

      { // Start loop over shape1

        /* Skip vertexes that have already been included */
        if (degen1[i] == 0)
            continue;

        for (size_t j = 0; j < shape2 -> dim; j++)

          { // Start loop over shape2

            /* Skip vertexes that have already been included */
            if (degen2[j] == 0)
                continue;

            /* Edges don't match */
            if (_avrShapeTri_ && _avrShapeLine_)
              {
                if (!shape_comp_vertex_bin(shape1, i, shape2, j, (covArg -> dk1 + covArg -> dk2) / 2., (covArg -> dmu1 + covArg -> dmu2) / 2.))
                    continue;
              }

            else
              {
                if (!shape_comp_vertex(shape1, i, shape2, j))
                    continue;
              }

            /* Set the line variables */
            for (size_t n1 = 0; n1 < shape1 -> dim; n1++)
              {
                kernels_qset_k(kern, n1, shape_get_vertex_length(shape1, (i + n1) % shape1 -> dim));
                kernels_qset_mu(kern, n1, shape_get_vertex_orientation(shape1, (i + n1) % shape1 -> dim));

                for (size_t n2 = n1 + 1; n2 < shape1 -> dim; n2++)
                    kernels_qset_nu(kern, n1, n2, shape_get_vertex_angle(shape1, (i + n1) % shape1 -> dim, (i + n2) % shape1 -> dim));

              }

            /* Set the triangle variables */
            for (size_t n1 = 0; n1 < shape2 -> dim; n1++)
              {
                kernels_qset_k(kern, n1 + shape1 -> dim, shape_get_vertex_length(shape2, (j + n1) % shape2 -> dim));
                kernels_qset_mu(kern, n1 + shape1 -> dim, shape_get_vertex_orientation(shape2, (j + n1) % shape2 -> dim));

                for (size_t n2 = n1 + 1; n2 < shape2 -> dim; n2++)
                    kernels_qset_nu(kern, n1 + shape1 -> dim, n2 + shape1 -> dim, shape_get_vertex_angle(shape2, (j + n1) % shape2 -> dim, (j + n2) % shape2 -> dim));
              }

            /* Calculate the result */
            if (_avrShapeTri_ && _avrShapeLine_)
                avr_cov_pb_ngauss_bp(covArg, resultBP);

            else
                avr_cov_pb_ngauss_bp_inf(covArg, resultBP);

            resultBP[0] *= (double) (degen1[i] * degen2[j]) * factor;
            resultBP[1] *= (double) (degen1[i] * degen2[j]) * factor;

            /* Increment the total result */
            result[0] += resultBP[0];
            result[1] = sqrt(result[1]*result[1] + resultBP[1]*resultBP[1]);

          } // End loop over shape2

      } // End loop over shape1

    /* Reset variables */
    for (size_t i = 0; i < shape1 -> dim; i++)
      {
        kernels_qset_k(kern, i, shape_get_vertex_length(shape1, i));
        kernels_qset_mu(kern, i, shape_get_vertex_orientation(shape1, i));

        for (size_t j = i + 1; j < shape1 -> dim; j++)
          {
            kernels_qset_nu(kern, i, j, shape_get_vertex_angle(shape1, i, j));
          }
      }

    for (size_t i = 0; i < shape2 -> dim; i++)
      {
        kernels_qset_k(kern, i + shape1 -> dim, shape_get_vertex_length(shape2, i));
        kernels_qset_mu(kern, i + shape1 -> dim, shape_get_vertex_orientation(shape2, i));

        for (size_t j = i + 1; j < shape2 -> dim; j++)
          {
            kernels_qset_nu(kern, i + shape1 -> dim, j + shape1 -> dim, shape_get_vertex_angle(shape2, i, j));
          }
      }

    /* Free memory */
    free(degen1);
    free(degen2);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int cov_pb_full(cov_arg_t *covArg, double *result)
{
    /*

        Calculate the full PB covariance matrix element

    */

    /* Gaussian pb covariance matrix */
    double covpbGauss[3];
    cov_pb_gauss(covArg, covpbGauss);

    /* Non-Gaussian pb covariance matrix */
    double covpbNGauss[3];
    cov_pb_ngauss(covArg, covpbNGauss);

    /* Final result */
    result[0] = covpbGauss[0] + covpbNGauss[0];
    result[1] = sqrt(covpbGauss[1]*covpbGauss[1] + covpbNGauss[1]*covpbNGauss[1]);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*cov_pb(const char *label))(cov_arg_t*, double*)
{
    /*

        Return the requested function that calculates part (or full) of the pb covariance matrix

    */

    /* Gaussian PB Covariance Matrix */
    if (!strcmp(label, _covpbLabelGauss))
      {
        return cov_pb_gauss;
      }

    /* Non-Gaussian PB Covariance Matrix */
    if (!strcmp(label, _covpbLabelNGauss))
      {
        return cov_pb_ngauss;
      }

    /* Full PB Covariance Matrix */
    if (!strcmp(label, _covpbLabelFull))
      {
        return cov_pb_full;
      }

    /* Label not found */
    return _cov_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   (Indirect) PB Covariance Matrix   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _cov_pb(void *var, void *params)
{
    /*

        Calculate the pb covariance matrix

    */

    (void) params;

//    if (_covpbInMat == NULL)
      {
        double result[3];
        cov_pb_full(var, result);

        return result[0];
      }

    /* TODO : Must get sample args from file as well... (to know which index the ki, mui belong to in the matrix) */

    printf("Not yet supborted!\n");
    exit(1);

    return NAN;
}


static mat_t *_cov_pb_mat(void *var, void *params)
{
    /*

        TODO: Get the full PB covariance matrix

    */

    (void) params;

    cov_mat_t *covMat = (cov_mat_t*) var;

    /* Must first get the covariance matrix */
    if (_covpbInMat == NULL || (_covpbInMat -> dim[0] != _sampleRedshift_ -> sampleArg -> size || _covpbInMat -> mblock -> bdim[0][0] != _sampleShapeLine_ -> sizeFull || _covpbInMat -> mblock -> bdim[1][0] != _sampleShapeTri_ -> sizeFull))
      {
        mats_t *cov = cov_mat_poly(covMat);

        _covpbInMat = mat_free(_covpbInMat);
        _covpbInMat = mat_cp(mats_get_mat(cov, 0));

        cov = mats_free_full(cov);
      }

    return _covpbInMat;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     POLY SPECTRUM COVARIANCE MATRIX     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------   (Direct) Poly Spectrum Covariance Matrix   -----------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int (*cov_poly(const char *covLabel, const char *label))(cov_arg_t*, double*)
{
    /*

        Get the 'label' part of the poly spectrum covariance matrix with label 'covLabel'.

    */

    /* PP Covariance Matrix */
    if (!strcmp(covLabel, _idCovPP_))
      {
        return cov_pp(label);
      }

    /* BB Covariance Matrix */
    if (!strcmp(covLabel, _idCovBB_))
      {
        return cov_bb(label);
      }

    /* PB Covariance Matrix */
    if (!strcmp(covLabel, _idCovPB_))
      {
        return cov_pb(label);
      }


    /**  Label not found  **/

    return _cov_zero;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------   (Indirect) Poly Spectrum Covariance Matrix   ----------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double (*_cov_poly(const char *covLabel, const char *label))(void*, void*)
{
    /*

        Get the 'label' part of the poly spectrum covariance matrix with label 'covLabel'.

    */

    (void) label;


    /* PP Covariance Matrix */
    if (!strcmp(covLabel, _idCovPP_))
      {
        return _covPP_;
      }

    /* BB Covariance Matrix */
    if (!strcmp(covLabel, _idCovBB_))
      {
        return _covBB_;
      }

    /* PB Covariance Matrix */
    if (!strcmp(covLabel, _idCovPB_))
      {
        return _covPB_;
      }


    /**  Label not found  **/

    return _zeroFunc_;
}


static mat_t *(*_cov_poly_mat(const char *covLabel))(void *var, void *params)
{
    /*

        TODO: Get the full covariance matrix of given label

    */

    /* PP Covariance Matrix */
    if (!strcmp(covLabel, _idCovPP_))
      {
        return _covPPMat_;
      }

    /* BB Covariance Matrix */
    if (!strcmp(covLabel, _idCovBB_))
      {
        return _covBBMat_;
      }

    /* PB Covariance Matrix */
    if (!strcmp(covLabel, _idCovPB_))
      {
        return _covPBMat_;
      }

    return NULL;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

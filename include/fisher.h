#ifndef FISHER_H_INCLUDED
#define FISHER_H_INCLUDED

#include "common.h"

#include "misc.h"

#include "dat.h"
#include "mat.h"

#include "shape.h"
#include "sample.h"
#include "interpolate.h"

#include "average.h"

#include "spectra.h"
#include "covariance.h"

#include "print.h"



typedef struct
{
    /*

        Output fisher matrix struct

    */

    /* Labels of the parameters */
    size_t size;
    char **labels;

    /* Logarithmic dependence on the parameters */
    bool *log;

    /* File */
    char *file;

    /* Binary */
    bool binary;

    /* Precision (ignored if binary is true) */
    int precision;

} fish_out_t;


typedef struct
{
    /*

        Fisher matrix struct

    */

    /* Id */
    const char *id;

    /* Output parameters */
    fish_out_t *out;

    /* Print parameters */
    print_flags_t *flags;

} fish_mat_t;



/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%    EXTERN    %%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


extern bool flss_get_avr_shape_flag(
    const char *id);


extern sample_raw_t *flss_get_sample_redshift(void);

extern sample_shape_t *flss_get_sample_shape(
    const char *id);



/*  ----------------------------------------------------  */
/*  -------------   Global Identifiers   ---------------  */
/*  ----------------------------------------------------  */


/**  Fisher Parameter Identifiers  **/

extern const char *_idParamsA2Ga_;
extern const char *_idParamsD2Ga_;
extern const char *_idParamsA3GaA_;
extern const char *_idParamsA3GaB_;
extern const char *_idParamsD3GaA_;
extern const char *_idParamsD3GaB_;

extern const char *_idParamsB1_;
extern const char *_idParamsB2_;
extern const char *_idParamsC2Ga_;
extern const char *_idParamsBGam3_;

extern const char *_idParamsF_;
extern const char *_idParamsSigv_;

extern const char *_idParamsD_;
extern const char *_idParamsH_;

extern const char *_idParamsC0_;
extern const char *_idParamsC2_;
extern const char *_idParamsC4_;

extern const char *_idParamsPsn_;
extern const char *_idParamsBsn1_;
extern const char *_idParamsBsn2_;

extern const char *_idParamsPk_;


/**  Fisher Parameter Multiplicities  **/

extern const int _multParamsA2Ga_[__ID_VAR_SIZE__];
extern const int _multParamsD2Ga_[__ID_VAR_SIZE__];
extern const int _multParamsA3GaA_[__ID_VAR_SIZE__];
extern const int _multParamsA3GaB_[__ID_VAR_SIZE__];
extern const int _multParamsD3GaA_[__ID_VAR_SIZE__];
extern const int _multParamsD3GaB_[__ID_VAR_SIZE__];

extern const int _multParamsB1_[__ID_VAR_SIZE__];
extern const int _multParamsB2_[__ID_VAR_SIZE__];
extern const int _multParamsC2Ga_[__ID_VAR_SIZE__];
extern const int _multParamsBGam3_[__ID_VAR_SIZE__];

extern const int _multParamsF_[__ID_VAR_SIZE__];
extern const int _multParamsSigv_[__ID_VAR_SIZE__];

extern const int _multParamsD_[__ID_VAR_SIZE__];
extern const int _multParamsH_[__ID_VAR_SIZE__];

extern const int _multParamsC0_[__ID_VAR_SIZE__];
extern const int _multParamsC2_[__ID_VAR_SIZE__];
extern const int _multParamsC4_[__ID_VAR_SIZE__];

extern const int _multParamsPsn_[__ID_VAR_SIZE__];
extern const int _multParamsBsn1_[__ID_VAR_SIZE__];
extern const int _multParamsBsn2_[__ID_VAR_SIZE__];

extern const int _multParamsPk_[__ID_VAR_SIZE__];


/**  Spectra Identifiers  **/

extern const char *_idSpecDPnl_;
extern const char *_idSpecDBtr_;


/**  Covariance Matrix Identifiers  **/

extern const char *_idCovPP_;
extern const char *_idCovPB_;
extern const char *_idCovBB_;


/**  Fisher Matrix Identifiers  **/

extern const char *_idFishPP_;
extern const char *_idFishPB_;
extern const char *_idFishBB_;


/**  Samples  **/

extern sample_raw_t *_sampleRedshift_;
extern sample_shape_t *_sampleShapeLine_;
extern sample_shape_t *_sampleShapeTri_;



/*  ----------------------------------------------------  */
/*  ---------------   Initialise/Free   ----------------  */
/*  ----------------------------------------------------  */


int fish_ini(void);
int fish_free(void);



/*  ----------------------------------------------------  */
/*  ---------   Fisher Matrix Info Functions   ---------  */
/*  ----------------------------------------------------  */


size_t fish_info_get_kern_order(
    const char *id);



/*  ----------------------------------------------------  */
/*  ------------   Fisher Matrix Structs   -------------  */
/*  ----------------------------------------------------  */


fish_out_t *fish_out_new(void);

fish_out_t *fish_out_free(
    fish_out_t *out);

fish_out_t *fish_out_cp(
    fish_out_t *out);

/*  ----------------------------------------------------  */

int fish_out_add_label(
    fish_out_t *out,
    const char *label);

int fish_out_rm_label(
    fish_out_t *out,
    const char *label);

int fish_out_set_log(
    fish_out_t *out,
    const char *label,
    bool log);


int fish_out_set_file(
    fish_out_t *out,
    const char *file);

int fish_out_set_precision(
    fish_out_t *out,
    int precision);

int fish_out_set_binary(
    fish_out_t *out,
    bool binary);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fish_mat_t *fish_mat_new(
    const char *id);

fish_mat_t *fish_mat_free(
    fish_mat_t *fishMat);

/*  ----------------------------------------------------  */

int fish_mat_set_out(
    fish_mat_t *fishMat,
    fish_out_t *out);


int fish_mat_set_flags(
    fish_mat_t *fishMat,
    print_flags_t *flags);

/*  ----------------------------------------------------  */

fish_out_t *fish_mat_get_out(
    fish_mat_t *fishMat);

print_flags_t *fish_mat_get_flags(
    fish_mat_t *fishMat);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fish_mat_t *fish_mat_pp_default(void);
fish_mat_t *fish_mat_bb_default(void);
fish_mat_t *fish_mat_pb_default(void);



/*  ----------------------------------------------------  */
/*  ----------------   Fisher Matrix   -----------------  */
/*  ----------------------------------------------------  */


mats_t *fish_mat_poly(
    fish_mat_t *fishMat);

mats_t *fish_mat_poly_eig(
    fish_mat_t *fishMat);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */



#endif // FISHER_H_INCLUDED

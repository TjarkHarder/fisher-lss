#ifndef COV_H_INCLUDED
#define COV_H_INCLUDED

#include "common.h"


#include "misc.h"

#include "mat.h"

#include "shape.h"
#include "sample.h"
#include "interpolate.h"

#include "kernels.h"

#include "spectra.h"

#include "print.h"



typedef struct
{
    /*

        Output covariance matrix struct

    */

    /* Number of output matrices */
    size_t size;

    /* Labels of the contributions */
    char **labels;

    /* Output the integration errors */
    bool *err;

    /* File */
    char *file;

    /* Binary */
    bool binary;

    /* Precision (ignored if binary is true) */
    int precision;

} cov_out_t;


typedef struct
{
    /*

        Covariance matrix struct

    */

    /* Id */
    const char *id;

    /* Output parameters */
    cov_out_t *out;

    /* Print parameters */
    print_flags_t *flags;

} cov_mat_t;


typedef struct
{
    /*

        Covariance matrix arguments struct

    */

    /* Id */
    const char *id;

    /* SpecArg */
    spec_arg_t *specArg;

    /* Shapes */
    size_t i1;
    size_t i2;

    shape_t *shape1;
    shape_t *shape2;

    /* Bin widths */
    double dk1;
    double dmu1;

    double dk2;
    double dmu2;

} cov_arg_t;




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
/*  -----------------   Average Flag   -----------------  */
/*  ----------------------------------------------------  */


extern bool _avrShapeLine_;
extern bool _avrShapeTri_;



/*  ----------------------------------------------------  */
/*  ----------------   Shape Average   -----------------  */
/*  ----------------------------------------------------  */


extern int avr_shape_line_vol_ana(
    spec_arg_t *specArg,
    double *result);



/*  ----------------------------------------------------  */
/*  -------   Direct Covariance Matrix Average   -------  */
/*  ----------------------------------------------------  */


extern double avr_cov_direct(
    int (*avr)(cov_arg_t*, double*),
    cov_arg_t *covArg);



/*  ----------------------------------------------------  */
/*  ---------   PP Covariance Matrix Average   ---------  */
/*  ----------------------------------------------------  */


extern int avr_cov_pp_gauss(
    cov_arg_t *covArg,
    double *result);

extern int avr_cov_pp_gauss_inf(
    cov_arg_t *covArg,
    double *result);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


extern int avr_cov_pp_ngauss(
    cov_arg_t *covArg,
    double *result);

extern int avr_cov_pp_ngauss_inf(
    cov_arg_t *covArg,
    double *result);



/*  ----------------------------------------------------  */
/*  ---------   BB Covariance Matrix Average   ---------  */
/*  ----------------------------------------------------  */


extern int avr_cov_bb_gauss(
    cov_arg_t *covArg,
    double *result);

extern int avr_cov_bb_gauss_inf(
    cov_arg_t *covArg,
    double *result);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


extern int avr_cov_bb_ngauss_tp(
    cov_arg_t *covArg,
    double *result);

extern int avr_cov_bb_ngauss_tp_inf(
    cov_arg_t *covArg,
    double *result);

/*  ----------------------------------------------------  */

extern int avr_cov_bb_ngauss_bb(
    cov_arg_t *covArg,
    double *result);

extern int avr_cov_bb_ngauss_bb_inf(
    cov_arg_t *covArg,
    double *result);


/*  ----------------------------------------------------  */
/*  ---------   BB Covariance Matrix Average   ---------  */
/*  ----------------------------------------------------  */


extern int avr_cov_pb_ngauss_bp(
    cov_arg_t *covArg,
    double *result);

extern int avr_cov_pb_ngauss_bp_inf(
    cov_arg_t *covArg,
    double *result);



/*  ----------------------------------------------------  */
/*  -------------   Global Identifiers   ---------------  */
/*  ----------------------------------------------------  */


/**  Miscellaneous Identifiers  **/

extern const char *_idIntError_;

extern const char *_idGauss_;
extern const char *_idNGauss_;
extern const char *_idFull_;


/**  Variable Identifiers  **/

extern const char *_idVarZ_;
extern const char *_idVarK_;
extern const char *_idVarMu_;


/**  Covariance Matrix Identifiers  **/

extern const char *_idCovPP_;
extern const char *_idCovBB_;
extern const char *_idCovPB_;


/**  Spectra Identifiers  **/

extern const char *_idSpecPnl_;
extern const char *_idSpecBtr_;
extern const char *_idSpecTtr_;


/**  Samples  **/

extern sample_raw_t *_sampleRedshift_;
extern sample_shape_t *_sampleShapeLine_;
extern sample_shape_t *_sampleShapeTri_;





/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%    INTERN    %%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  ---------------   Initialise/Free   ----------------  */
/*  ----------------------------------------------------  */


int cov_ini(void);
int cov_free(void);



/*  ----------------------------------------------------  */
/*  -------   Covariance Matrix Info Functions   -------  */
/*  ----------------------------------------------------  */


intgrt_t *cov_info_get_average_integrate(
    const char *id);

size_t cov_info_get_kern_order(
    const char *id);



/*  ----------------------------------------------------  */
/*  ---------   Covariance Matrix from File   ----------  */
/*  ----------------------------------------------------  */


int cov_in_pp_set_file(
    const char *file);

int cov_in_pp_set_label(
    const char *label);

/*  ----------------------------------------------------  */

int cov_in_pp_setup(size_t *size);
int cov_in_pp_reset(void);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int cov_in_bb_set_file(
    const char *file);

int cov_in_bb_set_label(
    const char *label);

/*  ----------------------------------------------------  */

int cov_in_bb_setup(size_t *size);
int cov_in_bb_reset(void);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int cov_in_pb_set_file(
    const char *file);

int cov_in_pb_set_label(
    const char *label);

/*  ----------------------------------------------------  */

int cov_in_pb_setup(size_t *size1, size_t *size2);
int cov_in_pb_reset(void);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


char *cov_in_get_label(
    const char *id);



/*  ----------------------------------------------------  */
/*  ----------   Covariance Matrix Structs   -----------  */
/*  ----------------------------------------------------  */



cov_arg_t *cov_arg_new(
    const char *id);

cov_arg_t *cov_arg_free(
    cov_arg_t *covArg);

/*  ----------------------------------------------------  */

int cov_arg_set_shape(
    cov_arg_t *covArg,
    shape_t *shape,
    size_t index);

int cov_arg_set_dk(
    cov_arg_t *covArg,
    double dk,
    size_t index);

int cov_arg_set_dmu(
    cov_arg_t *covArg,
    double dmu,
    size_t index);

/*  ----------------------------------------------------  */

spec_arg_t *cov_arg_get_spec_arg(
    cov_arg_t *covArg);

kern_t *cov_arg_get_kern(
    cov_arg_t *covArg);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


cov_out_t *cov_out_new(void);

cov_out_t *cov_out_free(
    cov_out_t *out);

cov_out_t *cov_out_cp(
    cov_out_t *out);

/*  ----------------------------------------------------  */

int cov_out_add_label(
    cov_out_t *out,
    const char *label);

int cov_out_rm_label(
    cov_out_t *out,
    const char *label);


int cov_out_set_err(
    cov_out_t *out,
    const char *label,
    bool err);


int cov_out_set_file(
    cov_out_t *out,
    const char *file);

int cov_out_set_precision(
    cov_out_t *out,
    int precision);

int cov_out_set_binary(
    cov_out_t *out,
    bool binary);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


cov_mat_t *cov_mat_new(
    const char *id);

cov_mat_t *cov_mat_free(
    cov_mat_t *covMat);

/*  ----------------------------------------------------  */

cov_out_t *cov_mat_get_out(
    cov_mat_t *covMat);

print_flags_t *cov_mat_get_flags(
    cov_mat_t *covMat);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


cov_mat_t *cov_mat_pp_default(void);
cov_mat_t *cov_mat_bb_default(void);
cov_mat_t *cov_mat_pb_default(void);



/*  ----------------------------------------------------  */
/*  ------------   PP Covariance Matrix   --------------  */
/*  ----------------------------------------------------  */


int cov_pvol_num(cov_arg_t *covArg, double *result);
int cov_pvol_ana(cov_arg_t *covArg, double *result);


int cov_pp_gauss(cov_arg_t *covArg, double *result);
int cov_pp_ngauss(cov_arg_t *covArg, double *result);
int cov_pp_full(cov_arg_t *covArg, double *result);

int (*cov_pp(const char *label))(cov_arg_t *covArg, double *result);


extern double (*_covPP_)(void*, void*);
extern mat_t *(*_covPPMat_)(void*, void*);



/*  ----------------------------------------------------  */
/*  ------------   BB Covariance Matrix   --------------  */
/*  ----------------------------------------------------  */


int cov_bvol_num(cov_arg_t *covArg, double *result);
int cov_bvol_ana(cov_arg_t *covArg, double *result);


int cov_bb_gauss(cov_arg_t *covArg, double *result);
int cov_bb_ngauss(cov_arg_t *covArg, double *result);
int cov_bb_full(cov_arg_t *covArg, double *result);

int (*cov_bb(const char *label))(cov_arg_t *covArg, double *result);


extern double (*_covBB_)(void*, void*);
extern mat_t *(*_covBBMat_)(void*, void*);



/*  ----------------------------------------------------  */
/*  ------------   PB Covariance Matrix   --------------  */
/*  ----------------------------------------------------  */


int cov_pb_gauss(cov_arg_t *covArg, double *result);
int cov_pb_ngauss(cov_arg_t *covArg, double *result);
int cov_pb_full(cov_arg_t *covArg, double *result);

int (*cov_pb(const char *label))(cov_arg_t *covArg, double *result);


extern double (*_covPB_)(void*, void*);
extern mat_t *(*_covPBMat_)(void*, void*);



/*  ----------------------------------------------------  */
/*  --------   Polyspectrum Covariance Matrix   --------  */
/*  ----------------------------------------------------  */


extern mat_t *(*(*_covPolyMat_)(const char*))(void*, void*);
extern double (*(*_covPoly_)(const char*, const char*))(void*, void*);

int (*cov_poly(const char *covLabel, const char *label))(cov_arg_t*, double*);

/*  ----------------------------------------------------  */

mats_t *cov_mat_poly(
    cov_mat_t *covMat);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // COV_H_INCLUDED

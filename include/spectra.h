#ifndef SPECTRA_H_INCLUDED
#define SPECTRA_H_INCLUDED

#include "common.h"


#include "misc.h"

#include "integrate.h"

#include "shape.h"
#include "sample.h"
#include "interpolate.h"

#include "kernels.h"
#include "integrand.h"

#include "print.h"



typedef struct
{
    /*

        Output spectra struct

    */

    /* Number of output columns */
    size_t size;

    /* Labels of the contributions */
    char **labels;

    /* Output the integration errors */
    bool *err;

    /* Logarithmic derivatives (if derivatives are calculated) */
    bool *log;

    /* File */
    char *file;

    /* Binary */
    bool binary;

    /* Precision (ignored if binary is true) */
    int precision;

} spec_out_t;


typedef struct
{
    /*

        Spectra struct

    */

    /* Id */
    const char *id;

    /* Output parameters */
    spec_out_t *out;

    /* Print parameters */
    print_flags_t *flags;

} spec_dat_t;




typedef struct
{
    /*

        Derivative variables

    */

    /* Id (needed?) */
    char *id;

    /* Dependence of parameter on z, k, mu variables */
    double z;
    double k;
    double mu;

    /* Logarithmic derivative (i.e. d/dlog X = X d/dX) */
    bool log;

} spec_deriv_t;


typedef struct
{
    /*

        Spectrum arguments struct

    */

    /* Id */
    const char *id;

    /* Kernel */
    kern_t *kern;

    /* Deriv */
    spec_deriv_t *deriv;

    /* Bin widths */
    double dz;
    double dk;
    double dmu;

    double binVolume;

} spec_arg_t;



/*  ----------------------------------------------------  */
/*  -------------   Global Identifiers   ---------------  */
/*  ----------------------------------------------------  */


/**  Miscellaneous Identifiers  **/

extern const char *_idIntError_;

extern const char *_idTree_;
extern const char *_idCtr_;
extern const char *_idSn_;

extern const char *_idFull_;


/**  Variable Identifiers  **/

extern const char *_idVarZ_;
extern const char *_idVarK_;
extern const char *_idVarMu_;


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

extern const char *_idSpecPnl_;
extern const char *_idSpecDPnl_;

extern const char *_idSpecBtr_;
extern const char *_idSpecDBtr_;

extern const char *_idSpecTtr_;


/**  Samples  **/

extern sample_raw_t *_sampleRedshift_;
extern sample_shape_t *_sampleShapeLine_;
extern sample_shape_t *_sampleShapeTri_;



/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


extern bool flss_get_avr_shape_flag(
    const char *id);


extern sample_raw_t *flss_get_sample_redshift(void);

extern sample_shape_t *flss_get_sample_shape(
    const char *id);



/*  ----------------------------------------------------  */
/*  ---------------   Initialise/Free   ----------------  */
/*  ----------------------------------------------------  */


int spec_ini(void);
int spec_free(void);



/*  ----------------------------------------------------  */
/*  ------------   Spectra Info Functions   ------------  */
/*  ----------------------------------------------------  */


bool spec_info_isderiv(
    const char *id);

/*  ----------------------------------------------------  */

int spec_info_set_loop_order(
    const char *id,
    size_t loopOrder);

/*  ----------------------------------------------------  */

size_t spec_info_get_order(
    const char *id);


size_t spec_info_get_loop_order(
    const char *id);

intgrt_t *spec_info_get_loop_integrate(
    const char *id,
    size_t loopOrder);


size_t spec_info_get_kern_order(
    const char *id);



/*  ----------------------------------------------------  */
/*  --------------   Spectra from File   ---------------  */
/*  ----------------------------------------------------  */


int spec_in_pnl_set_file(
    const char *file);

int spec_in_pnl_set_label(
    const char *label);

/*  ----------------------------------------------------  */

int spec_in_dpnl_set_file(
    const char *file);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int spec_in_pnl_setup(void);

/*  ----------------------------------------------------  */

int spec_in_dpnl_setup(void);



/*  ----------------------------------------------------  */
/*  ---------------   Spectra Structs   ----------------  */
/*  ----------------------------------------------------  */


spec_deriv_t *spec_deriv_new(void);

spec_deriv_t *spec_deriv_free(
    spec_deriv_t *deriv);

/*  ----------------------------------------------------  */

int spec_deriv_set_var(
    spec_deriv_t *deriv,
    const char *label,
    double value);


int spec_deriv_set_z(
    spec_deriv_t *deriv,
    double z);

int spec_deriv_set_k(
    spec_deriv_t *deriv,
    double k);

int spec_deriv_set_mu(
    spec_deriv_t *deriv,
    double mu);


int spec_deriv_set_log(
    spec_deriv_t *deriv,
    bool log);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


spec_arg_t *spec_arg_new(
    const char *id);

spec_arg_t *spec_arg_free(
    spec_arg_t *specArg);

/*  ----------------------------------------------------  */


int spec_arg_set_dz(
    spec_arg_t *specArg,
    double dz);

int spec_arg_set_dk(
    spec_arg_t *specArg,
    double dk);

int spec_arg_set_dmu(
    spec_arg_t *specArg,
    double dmu);


int spec_arg_set_volume(
    spec_arg_t *specArg,
    double binVolume);


int spec_arg_set_kern(
    spec_arg_t *specArg,
    kern_t *kern);

int spec_arg_set_deriv(
    spec_arg_t *specArg,
    spec_deriv_t *deriv);

/*  ----------------------------------------------------  */

kern_t *spec_arg_get_kern(
    spec_arg_t *specArg);

spec_deriv_t *spec_arg_get_deriv(
    spec_arg_t *specArg);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


spec_out_t *spec_out_new(void);

spec_out_t *spec_out_free(
    spec_out_t *out);

spec_out_t *spec_out_cp(
    spec_out_t *out);

/*  ----------------------------------------------------  */

int spec_out_add_label(
    spec_out_t *out,
    const char *label);

int spec_out_rm_label(
    spec_out_t *out,
    const char *label);


int spec_out_set_err(
    spec_out_t *out,
    const char *label,
    bool err);

int spec_out_set_log(
    spec_out_t *out,
    const char *label,
    bool log);


int spec_out_set_file(
    spec_out_t *out,
    const char *file);

int spec_out_set_precision(
    spec_out_t *out,
    int precision);

int spec_out_set_binary(
    spec_out_t *out,
    bool binary);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


spec_dat_t *spec_dat_new(
    const char *id);

spec_dat_t *spec_dat_free(
    spec_dat_t *specDat);

/*  ----------------------------------------------------  */

int spec_dat_set_out(
    spec_dat_t *specDat,
    spec_out_t *out);


int spec_dat_set_flags(
    spec_dat_t *specDat,
    print_flags_t *flags);

/*  ----------------------------------------------------  */

spec_out_t *spec_dat_get_out(
    spec_dat_t *specDat);

print_flags_t *spec_dat_get_flags(
    spec_dat_t *specDat);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


spec_dat_t *spec_dat_pnl_default(void);
spec_dat_t *spec_dat_dpnl_default(void);

spec_dat_t *spec_dat_btr_default(void);
spec_dat_t *spec_dat_dbtr_default(void);



/*  ----------------------------------------------------  */
/*  ----------------   Power Spectrum   ----------------  */
/*  ----------------------------------------------------  */


/* Calculate the non-linear power spectrum */

extern double (*_specPtr_)(void*, void*);
extern double (*_specPnl_)(void*, void*);


int spec_pnl_tree(spec_arg_t *specArg, double *result);
int spec_pnl_p22(spec_arg_t *specArg, double *result);
int spec_pnl_p13(spec_arg_t *specArg, double *result);
int spec_pnl_ctr(spec_arg_t *specArg, double *result);
int spec_pnl_sn(spec_arg_t *specArg, double *result);
int spec_pnl_full(spec_arg_t *specArg, double *result);

int (*spec_pnl(const char *label))(spec_arg_t *specArg, double *result);


/* Calculate the non-linear power spectrum derivatives */

extern double (*_specDPnlA2Ga_)(void*, void*);
extern double (*_specDPnlD2Ga_)(void*, void*);
extern double (*_specDPnlA3GaA_)(void*, void*);
extern double (*_specDPnlA3GaB_)(void*, void*);
extern double (*_specDPnlD3GaA_)(void*, void*);
extern double (*_specDPnlD3GaB_)(void*, void*);

extern double (*_specDPnlB1_)(void*, void*);
extern double (*_specDPnlB2_)(void*, void*);
extern double (*_specDPnlC2Ga_)(void*, void*);
extern double (*_specDPnlBGam3_)(void*, void*);

extern double (*_specDPnlF_)(void*, void*);
extern double (*_specDPnlSigv_)(void*, void*);

extern double (*_specDPnlD_)(void*, void*);
extern double (*_specDPnlH_)(void*, void*);

extern double (*_specDPnlC0_)(void*, void*);
extern double (*_specDPnlC2_)(void*, void*);
extern double (*_specDPnlC4_)(void*, void*);

extern double (*_specDPnlPsn_)(void*, void*);
extern double (*_specDPnlBsn1_)(void*, void*);
extern double (*_specDPnlBsn2_)(void*, void*);

extern double (*_specDPnlPk_)(void*, void*);

extern double (*(*_specDPnl_)(const char*))(void*, void*);


int spec_dpnl_a2ga(spec_arg_t *specArg, double *result);
int spec_dpnl_d2ga(spec_arg_t *specArg, double *result);
int spec_dpnl_a3gaa(spec_arg_t *specArg, double *result);
int spec_dpnl_a3gab(spec_arg_t *specArg, double *result);
int spec_dpnl_d3gaa(spec_arg_t *specArg, double *result);
int spec_dpnl_d3gab(spec_arg_t *specArg, double *result);

int spec_dpnl_b1(spec_arg_t *specArg, double *result);
int spec_dpnl_b2(spec_arg_t *specArg, double *result);
int spec_dpnl_c2ga(spec_arg_t *specArg, double *result);
int spec_dpnl_bam3(spec_arg_t *specArg, double *result);

int spec_dpnl_f(spec_arg_t *specArg, double *result);
int spec_dpnl_sigv(spec_arg_t *specArg, double *result);

int spec_dpnl_d(spec_arg_t *specArg, double *result);
int spec_dpnl_h(spec_arg_t *specArg, double *result);

int spec_dpnl_c0(spec_arg_t *specArg, double *result);
int spec_dpnl_c2(spec_arg_t *specArg, double *result);
int spec_dpnl_c4(spec_arg_t *specArg, double *result);

int spec_dpnl_psn(spec_arg_t *specArg, double *result);
int spec_dpnl_bsn1(spec_arg_t *specArg, double *result);
int spec_dpnl_bsn2(spec_arg_t *specArg, double *result);

int spec_dpnl_pk(spec_arg_t *specArg, double *result);

int (*spec_dpnl(const char *label))(spec_arg_t*, double*);



/*  ----------------------------------------------------  */
/*  ------------------   Bispectrum   ------------------  */
/*  ----------------------------------------------------  */


/* Calculate the tree-level bispectrum */

extern double (*_specBtr_)(void*, void*);

int spec_btr_tree(spec_arg_t *specArg, double *result);
int spec_btr_sn(spec_arg_t *specArg, double *result);
int spec_btr_full(spec_arg_t *specArg, double *result);

int (*spec_btr(const char *label))(spec_arg_t *specArg, double *result);


/* Calculate the tree-level bispectrum derivatives */

extern double (*_specDBtrA2Ga_)(void*, void*);
extern double (*_specDBtrD2Ga_)(void*, void*);
extern double (*_specDBtrA3GaA_)(void*, void*);
extern double (*_specDBtrA3GaB_)(void*, void*);
extern double (*_specDBtrD3GaA_)(void*, void*);
extern double (*_specDBtrD3GaB_)(void*, void*);

extern double (*_specDBtrB1_)(void*, void*);
extern double (*_specDBtrB2_)(void*, void*);
extern double (*_specDBtrC2Ga_)(void*, void*);
extern double (*_specDBtrBGam3_)(void*, void*);

extern double (*_specDBtrF_)(void*, void*);
extern double (*_specDBtrSigv_)(void*, void*);

extern double (*_specDBtrD_)(void*, void*);
extern double (*_specDBtrH_)(void*, void*);

extern double (*_specDBtrC0_)(void*, void*);
extern double (*_specDBtrC2_)(void*, void*);
extern double (*_specDBtrC4_)(void*, void*);

extern double (*_specDBtrPsn_)(void*, void*);
extern double (*_specDBtrBsn1_)(void*, void*);
extern double (*_specDBtrBsn2_)(void*, void*);

extern double (*_specDBtrPk_)(void*, void*);

extern double (*(*_specDBtr_)(const char*))(void*, void*);


int spec_dbtr_a2ga(spec_arg_t *specArg, double *result);
int spec_dbtr_d2ga(spec_arg_t *specArg, double *result);
int spec_dbtr_a3gaa(spec_arg_t *specArg, double *result);
int spec_dbtr_a3gab(spec_arg_t *specArg, double *result);
int spec_dbtr_d3gaa(spec_arg_t *specArg, double *result);
int spec_dbtr_d3gab(spec_arg_t *specArg, double *result);

int spec_dbtr_b1(spec_arg_t *specArg, double *result);
int spec_dbtr_b2(spec_arg_t *specArg, double *result);
int spec_dbtr_c2ga(spec_arg_t *specArg, double *result);
int spec_dbtr_bam3(spec_arg_t *specArg, double *result);

int spec_dbtr_f(spec_arg_t *specArg, double *result);
int spec_dbtr_sigv(spec_arg_t *specArg, double *result);

int spec_dbtr_d(spec_arg_t *specArg, double *result);
int spec_dbtr_h(spec_arg_t *specArg, double *result);

int spec_dbtr_c0(spec_arg_t *specArg, double *result);
int spec_dbtr_c2(spec_arg_t *specArg, double *result);
int spec_dbtr_c4(spec_arg_t *specArg, double *result);

int spec_dbtr_psn(spec_arg_t *specArg, double *result);
int spec_dbtr_bsn1(spec_arg_t *specArg, double *result);
int spec_dbtr_bsn2(spec_arg_t *specArg, double *result);

int spec_dbtr_pk(spec_arg_t *specArg, double *result);

int (*spec_dbtr(const char *label))(spec_arg_t*, double*);



/*  ----------------------------------------------------  */
/*  ------------------   Trispectrum   -----------------  */
/*  ----------------------------------------------------  */


/* Calculate the tree-level trispectrum */

extern double (*_specTtr_)(void*, void*);

int spec_ttr_tree(spec_arg_t *specArg, double *result);
int spec_ttr_sn(spec_arg_t *specArg, double *result);
int spec_ttr_full(spec_arg_t *specArg, double *result);

int (*spec_ttr(const char *label))(spec_arg_t *specArg, double *result);



/*  ----------------------------------------------------  */
/*  -----------------   Polyspectrum   -----------------  */
/*  ----------------------------------------------------  */


extern double (*(*_specPoly_)(const char*, const char*))(void*, void*);

int (*spec_poly(const char *specLabel, const char *label))(spec_arg_t*, double*);

/*  ----------------------------------------------------  */

dat_t *spec_dat_poly(
    spec_dat_t *spec);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */



#endif // SPECTRA_H_INCLUDED

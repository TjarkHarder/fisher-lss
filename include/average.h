#ifndef AVERAGE_H_INCLUDED
#define AVERAGE_H_INCLUDED

#include "common.h"


#include "misc.h"

#include "integrate.h"

#include "shape.h"
#include "sample.h"
#include "interpolate.h"

#include "integrand-average.h"

#include "kernels.h"

#include "spectra.h"
#include "covariance.h"




/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%    EXTERN    %%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  -----------------   Identifiers   ------------------  */
/*  ----------------------------------------------------  */


/**  Miscellaneuos Identifiers  **/

extern const char *_idInf_;


/**  Shape Identifiers  **/

extern const char *_idShapeLine_;
extern const char *_idShapeTri_;


/**  Average Identifiers  **/

extern const char *_idAvrLine_;
extern const char *_idAvrTri_;

extern const char *_idAvrLineVol_;
extern const char *_idAvrTriVol_;

extern const char *_idAvrLineVolFuncNum_;
extern const char *_idAvrLineVolFuncAna_;

extern const char *_idAvrTriVolFuncNum_;
extern const char *_idAvrTriVolFuncAna_;

extern const char *_idAvrCovPPGauss_;
extern const char *_idAvrCovPPNGaussInf_;
extern const char *_idAvrCovPPNGauss_;

extern const char *_idAvrCovBBGauss_;
extern const char *_idAvrCovBBNGauss_;

extern const char *_idAvrCovPBNGauss_;





/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%    INTERN    %%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */



/*  ----------------------------------------------------  */
/*  -------------------   Setters   --------------------  */
/*  ----------------------------------------------------  */


int avr_set_shape_vol_func(const char *idShape, const char *idFunc);



/*  ----------------------------------------------------  */
/*  -------------------   Getters   --------------------  */
/*  ----------------------------------------------------  */


intgrt_t *avr_get_integrate(const char *id);



/*  ----------------------------------------------------  */
/*  ---------------   Initialise/Free   ----------------  */
/*  ----------------------------------------------------  */


int avr_ini(void);
int avr_free(void);





/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%    SHAPE AVERAGE    %%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  ---------------   Direct Function   ----------------  */
/*  ----------------------------------------------------  */


double avr_shape_direct(
    int (*avr)(double (*)(void*, void*), spec_arg_t*, double*),
    double (*specFunc)(void*, void*),
    spec_arg_t *specArg);

double avr_shape_vol_direct(
    int (*vol)(spec_arg_t*, double*),
    spec_arg_t *specArg);



/*  ----------------------------------------------------  */
/*  --------------   Line-Bin Average   ----------------  */
/*  ----------------------------------------------------  */


int avr_shape_line(
    double (*specFunc)(void*, void*),
    spec_arg_t *specArg,
    double *result);

/*  ----------------------------------------------------  */

int avr_shape_line_vol(
    spec_arg_t *specArg,
    double *result);

int avr_shape_line_vol_ana(
    spec_arg_t *specArg,
    double *result);



/*  ----------------------------------------------------  */
/*  ------------   Triangle-Bin Average   --------------  */
/*  ----------------------------------------------------  */


int avr_shape_tri(
    double (*specFunc)(void*, void*),
    spec_arg_t *specArg,
    double *result);

/*  ----------------------------------------------------  */

int avr_shape_tri_vol(
    spec_arg_t *specArg,
    double *result);

int avr_shape_tri_vol_ana(
    spec_arg_t *specArg,
    double *result);



/*  ----------------------------------------------------  */
/*  ---   Average over Infinitesimally Small Bins   ----  */
/*  ----------------------------------------------------  */


int avr_shape_inf(
    double (*specFunc)(void*, void*),
    spec_arg_t *specArg,
    double *result);



/*  ----------------------------------------------------  */
/*  -------   Get Average and Volume Functions   -------  */
/*  ----------------------------------------------------  */


int (*avr_shape_get_func(size_t order))(double(*)(void*, void*), spec_arg_t*, double*);
int (*avr_shape_vol_get_func(size_t order))(spec_arg_t*, double*);



/*  ----------------------------------------------------  */
/*  -----   Shape Bin-Volume for Sample of Shapes   ----  */
/*  ----------------------------------------------------  */


int avr_shape_sample_vol(
    sample_shape_t *sampleShape);





/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%    COVARIANCE MATRIX AVERAGE    %%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  -------   Direct Covariance Matrix Average   -------  */
/*  ----------------------------------------------------  */


double avr_cov_direct(
    int (*avr)(cov_arg_t*, double*),
    cov_arg_t *covArg);



/*  ----------------------------------------------------  */
/*  ---------   PP Covariance Matrix Average   ---------  */
/*  ----------------------------------------------------  */


int avr_cov_pp_gauss(
    cov_arg_t *covArg,
    double *result);

int avr_cov_pp_gauss_inf(
    cov_arg_t *covArg,
    double *result);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int avr_cov_pp_ngauss(
    cov_arg_t *covArg,
    double *result);

int avr_cov_pp_ngauss_inf(
    cov_arg_t *covArg,
    double *result);



/*  ----------------------------------------------------  */
/*  ---------   BB Covariance Matrix Average   ---------  */
/*  ----------------------------------------------------  */


int avr_cov_bb_gauss(
    cov_arg_t *covArg,
    double *result);

int avr_cov_bb_gauss_inf(
    cov_arg_t *covArg,
    double *result);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int avr_cov_bb_ngauss_tp(
    cov_arg_t *covArg,
    double *result);

int avr_cov_bb_ngauss_tp_inf(
    cov_arg_t *covArg,
    double *result);

/*  ----------------------------------------------------  */

int avr_cov_bb_ngauss_bb(
    cov_arg_t *covArg,
    double *result);

int avr_cov_bb_ngauss_bb_inf(
    cov_arg_t *covArg,
    double *result);



/*  ----------------------------------------------------  */
/*  ---------   PB Covariance Matrix Average   ---------  */
/*  ----------------------------------------------------  */


int avr_cov_pb_ngauss_bp(
    cov_arg_t *covArg,
    double *result);

int avr_cov_pb_ngauss_bp_inf(
    cov_arg_t *covArg,
    double *result);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // AVERAGE_H_INCLUDED

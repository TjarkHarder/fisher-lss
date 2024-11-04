#ifndef FLSS_H_INCLUDED
#define FLSS_H_INCLUDED

#include "common.h"

#include "fiducials.h"

#include "misc.h"

#include "dat.h"
#include "mat.h"

#include "kernels.h"

#include "interpolate.h"

#include "integrand.h"
#include "integrand-average.h"
#include "integrate.h"

#include "average.h"

#include "spectra.h"
#include "covariance.h"
#include "fisher.h"

#include "print.h"



/*  ----------------------------------------------------  */
/*  ----------------   Miscellaneous   -----------------  */
/*  ----------------------------------------------------  */


extern const bool _true_;
extern const bool _false_;

extern double _zeroFunc_();
extern double _oneFunc_();




/*  ----------------------------------------------------  */
/*  ---------------   Random Numbers   -----------------  */
/*  ----------------------------------------------------  */


extern gsl_rng *_randomNumberGen_;



/*  ----------------------------------------------------  */
/*  ----------------   Global Labels   -----------------  */
/*  ----------------------------------------------------  */


/**  Miscellaneous Identifiers  **/

extern const char *_idIntError_;

extern const char *_idTree_;
extern const char *_idCtr_;
extern const char *_idSn_;

extern const char *_idGauss_;
extern const char *_idNGauss_;

extern const char *_idFull_;

extern const char *_idInf_;


/**  Integration Routines  **/

extern const char *_idIntgrtCQUAD_;
extern const char *_idIntgrtVegas_;
extern const char *_idIntgrtDivonne_;


/**  Shape Identifiers  **/

extern const char *_idShapeLine_;
extern const char *_idShapeTri_;


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


/**  Covariance Matrix Identifiers  **/

extern const char *_idCovPP_;
extern const char *_idCovBB_;
extern const char *_idCovPB_;


/**  Fisher Matrix Identifiers  **/

extern const char *_idFishPP_;
extern const char *_idFishBB_;
extern const char *_idFishPB_;


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
/*  -----------------   Global Flags   -----------------  */
/*  ----------------------------------------------------  */


/* Bias */
extern bool _fidInclBias_;

/* Redshift Space Distortions */
extern bool _fidInclRSD_;

/* Counter terms */
extern bool _fidInclCtr_;

/* Shot noise */
extern bool _fidInclSn_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Line Bin-Average */
extern bool _avrShapeLine_;

/* Triangle Bin-Average */
extern bool _avrShapeTri_;



/*  ----------------------------------------------------  */
/*  -------------------   Samples   --------------------  */
/*  ----------------------------------------------------  */


/* Redshift (Raw) Sample */
extern sample_raw_t *_sampleRedshift_;

/* Line (Shape) Sample */
extern sample_shape_t *_sampleShapeLine_;

/* Triangle (Shape) Sample */
extern sample_shape_t *_sampleShapeTri_;



/*  ----------------------------------------------------  */
/*  -------------   Calculate Fiducials   --------------  */
/*  ----------------------------------------------------  */


/* LCDM */

extern fid_lcdm_t *(*_fidLCDM_)(void*, void*);
extern void *_fidParamsLCDM_;

/*  ----------------------------------------------------  */

extern double (*_fidLCDMxOmegaM0_)(void*, void*);
extern void *_fidParamsLCDMxOmegaM0_;

extern double (*_fidLCDMxNs_)(void*, void*);
extern void *_fidParamsLCDMxNs_;

extern double (*_fidLCDMxGrowthIndex_)(void*, void*);
extern void *_fidParamsLCDMxGrowthIndex_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Bootstrap */

extern fid_btst_t *(*_fidBTST_)(void*, void*);
extern void *_fidParamsBTST_;

/*  ----------------------------------------------------  */

extern double (*_fidBTSTxA2Ga_)(void*, void*);
extern void *_fidParamsBTSTxA2Ga_;

extern double (*_fidBTSTxD2Ga_)(void*, void*);
extern void *_fidParamsBTSTxD2Ga_;


extern double (*_fidBTSTxA3GaA_)(void*, void*);
extern void *_fidParamsBTSTxA3GaA_;

extern double (*_fidBTSTxA3GaB_)(void*, void*);
extern void *_fidParamsBTSTxA3GaB_;

extern double (*_fidBTSTxD3GaA_)(void*, void*);
extern void *_fidParamsBTSTxD3GaA_;

extern double (*_fidBTSTxD3GaB_)(void*, void*);
extern void *_fidParamsBTSTxD3GaB_;


extern double (*_fidBTSTxH_)(void*, void*);
extern void *_fidParamsBTSTxH_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Bias */

extern fid_bias_t *(*_fidBias_)(void*, void*);
extern void *_fidParamsBias_;

/*  ----------------------------------------------------  */

extern double (*_fidBiasxB1_)(void*, void*);
extern void *_fidParamsBiasxB1_;

extern double (*_fidBiasxB2_)(void*, void*);
extern void *_fidParamsBiasxB2_;


extern double (*_fidBiasxBG2_)(void*, void*);
extern void *_fidParamsBiasxBG2_;

extern double (*_fidBiasxC2Ga_)(void*, void*);
extern void *_fidParamsBiasxC2Ga_;


extern double (*_fidBiasxBGam3_)(void*, void*);
extern void *_fidParamsBiasxBGam3_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* RSD */

extern fid_rsd_t *(*_fidRSD_)(void*, void*);
extern void *_fidParamsRSD_;

/*  ----------------------------------------------------  */

extern double (*_fidRSDxF_)(void*, void*);
extern void *_fidParamsRSDxF_;


extern double (*_fidRSDxSigv_)(void*, void*);
extern void *_fidParamsRSDxSigv_;

extern double (*_fidRSDxSigs_)(void*, void*);
extern void *_fidParamsRSDxSigs_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Ctr */

extern fid_ctr_t *(*_fidCtr_)(void*, void*);
extern void *_fidParamsCtr_;

/*  ----------------------------------------------------  */

extern double (*_fidCtrxC0_)(void*, void*);
extern void *_fidParamsCtrxC0_;

extern double (*_fidCtrxC2_)(void*, void*);
extern void *_fidParamsCtrxC2_;

extern double (*_fidCtrxC4_)(void*, void*);
extern void *_fidParamsCtrxC4_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Survey */

extern fid_surv_t *(*_fidSurv_)(void*, void*);
extern void *_fidParamsSurv_;

/*  ----------------------------------------------------  */

extern double (*_fidSurvxV_)(void*, void*);
extern void *_fidParamsSurvxV_;

extern double (*_fidSurvxN_)(void*, void*);
extern void *_fidParamsSurvxN_;

extern double (*_fidSurvxSn_)(void*, void*);
extern void *_fidParamsSurvxSn_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Growth */

extern double (*_fidGrowth_)(void*, void*);
extern void *_fidParamsGrowth_;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


/* Linear Power Spectrum */

extern double (*_fidPk_)(void*, void*);
extern double (*_fidExtrapPk_)(void*, void*, void*);
extern void *_fidParamsPk_;

extern double (*_fidDPk_)(void*, void*);
extern double (*_fidExtrapDPk_)(void*, void*, void*);
extern void *_fidParamsDPk_;



/*  ----------------------------------------------------  */
/*  --------------   Calculate Spectra   ---------------  */
/*  ----------------------------------------------------  */


/* Tree-Level Power Spectrum */

extern double (*_specPtr_)(void*, void*);

/* Non-Linear Power Spectrum */

extern double (*_specPnl_)(void*, void*);


/* Non-Linear Power Spectrum Derivatives */

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


/* Tree-Level Bispectrum */

extern double (*_specBtr_)(void*, void*);


/* Tree-Level Bispectrum Derivatives */

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


/* Tree-Level Trispectrum */

extern double (*_specTtr_)(void*, void*);


/* Poly Spectrum */

extern double (*(*_specPoly_)(const char*, const char*))(void*, void*);



/*  ----------------------------------------------------  */
/*  ---------   Calculate Covariance Matrix   ----------  */
/*  ----------------------------------------------------  */


/* PP Covariance Matrix */

extern mat_t *(*_covPPMat_)(void*, void*);

extern double (*_covPP_)(void*, void*);


/* BB Covariance Matrix */

extern mat_t *(*_covBBMat_)(void*, void*);

extern double (*_covBB_)(void*, void*);


/* PB Covariance Matrix */

extern mat_t *(*_covPBMat_)(void*, void*);

extern double (*_covPB_)(void*, void*);


/* Poly Spectrum Covariance Matrix */

extern mat_t *(*(*_covPolyMat_)(const char*))(void*, void*);

extern double (*(*_covPoly_)(const char*, const char*))(void*, void*);



/*  ----------------------------------------------------  */
/*  ---------------   Global Functions   ---------------  */
/*  ----------------------------------------------------  */


/* Initialise everything */

int flss_ini(void);


/* Getters */

char *flss_get_dir(void);

bool flss_get_avr_shape_flag(const char *id);

sample_raw_t *flss_get_sample_redshift(void);
sample_shape_t *flss_get_sample_shape(const char *id);


/* Setters */

int flss_set_threads(unsigned int threads);

int flss_set_avr_shape_flag(const char *id, bool avr);

int flss_set_sample_redshift(sample_raw_t *sampleRedshift);
int flss_set_sample_shape(const char *id, sample_shape_t *sampleShape);

int flss_set_random_seed(size_t seed);


/* Free everything */

int flss_free(void);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */



#endif // FLSS_H_INCLUDED

#ifndef FIDUCIALS_H_INCLUDED
#define FIDUCIALS_H_INCLUDED


#include "common.h"

#include "interpolate.h"
#include "integrate.h"



/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%  Fiducial Structures  %%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */


typedef struct
{
    /*

        LCDM parameters

    */

    double omegaM0;
    double ns;
    double growthIndex;

} fid_lcdm_t;


typedef struct
{
    /*

        Beyond EDS parameters

    */

    double a2Ga;
    double d2Ga;

    double a3GaA;
    double a3GaB;

    double d3GaA;
    double d3GaB;

    double h;

} fid_btst_t;


typedef struct
{
    /*

        Bias parameters

    */

    double b1;
    double b2;

    double bG2;
    double c2Ga;

    double bGam3;

} fid_bias_t;


typedef struct
{
    /*

        Redshift Space Distortion parameters

    */

    double f;

    double sigs;
    double sigv;

} fid_rsd_t;


typedef struct
{
    /*

        Counter terms parameters

    */

    double c0;
    double c2;
    double c4;

} fid_ctr_t;


typedef struct
{
    /*

        Survey parameters

    */

    double n;
    double sn;

    double V;

} fid_surv_t;



/*  ----------------------------------------------------  */
/*  ---------------   Module Functions   ---------------  */
/*  ----------------------------------------------------  */


int fiducials_ini(void);

/*  ----------------------------------------------------  */

int fiducials_setup(void);

/*  ----------------------------------------------------  */

int fiducials_set_file(char type, const char *file);

/*  ----------------------------------------------------  */

int fiducials_free(void);



/*  ----------------------------------------------------  */
/*  --------------   Include Fiducials   ---------------  */
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
/*  ----------   New / Free / Cp Functions   -----------  */
/*  ----------------------------------------------------  */


fid_lcdm_t *fid_lcdm_new(void);

fid_lcdm_t *fid_lcdm_free(
    fid_lcdm_t *lcdm);

/*  ----------------------------------------------------  */

fid_lcdm_t *fid_lcdm_cp(
    fid_lcdm_t *lcdm);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fid_btst_t *fid_btst_new(void);

fid_btst_t *fid_btst_free(
    fid_btst_t *btst);

/*  ----------------------------------------------------  */

fid_btst_t *fid_btst_cp(
    fid_btst_t *btst);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fid_bias_t *fid_bias_new(void);
fid_bias_t *fid_bias_free(
    fid_bias_t *bias);

/*  ----------------------------------------------------  */

fid_bias_t *fid_bias_cp(
    fid_bias_t *bias);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fid_rsd_t *fid_rsd_new(void);
fid_rsd_t *fid_rsd_free(
    fid_rsd_t *rsd);

/*  ----------------------------------------------------  */

fid_rsd_t *fid_rsd_cp(
    fid_rsd_t *rsd);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fid_ctr_t *fid_ctr_new(void);

fid_ctr_t *fid_ctr_free(
    fid_ctr_t *ctr);

/*  ----------------------------------------------------  */

fid_ctr_t *fid_ctr_cp(
    fid_ctr_t *ctr);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


fid_surv_t *fid_surv_new(void);

fid_surv_t *fid_surv_free(
    fid_surv_t *surv);

/*  ----------------------------------------------------  */

fid_surv_t *fid_surv_cp(
    fid_surv_t *surv);



/*  ----------------------------------------------------  */
/*  -------   Calculate Fiducials (to be set)   --------  */
/*  ----------------------------------------------------  */


/*    These functions calculate the fiducials for some    */
/* redshift (first void* pointer). The functions will be  */
/* set by input.c but may be replaced without the need to */
/*             change anything in fiducials.c.            */
/*  The second void *pointer is a free parameter pointer  */
/*                free to be used if needed.              */


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
/*  ------------------   Cosmology   -------------------  */
/*  ----------------------------------------------------  */


double fid_omegam(double z, fid_lcdm_t *lcdm);
double fid_hubble(double z, fid_lcdm_t *lcdm);

double fid_growth_rate(double z, fid_lcdm_t *lcdm);
double fid_growth_fct(double z, fid_lcdm_t *lcdm);





#endif // FIDUCIALS_H_INCLUDED

#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <quadmath.h>

#include "omp.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>


#include "cuba.h"



/*  ----------------------------------------------------  */
/*  ---------------   Macro Definitions   --------------  */
/*  ----------------------------------------------------  */


/* Directory */
#ifndef __FISHERLSSDIR__
#define __FISHERLSSDIR__ "."
#endif

/* Sizes */
#define __ID_VAR_SIZE__ 3

#define __ID_PARAMS_SIZE__ 21

#define __ID_SPEC_SIZE__ 5
#define __ID_COV_SIZE__ 3
#define __ID_FISH_SIZE__ 3

/* Speed of light */
#ifndef __C_LIGHT__
#define __C_LIGHT__ 299792458.
#endif

#ifndef __ABSTOL__
#define __ABSTOL__ 1.e-12
#endif

#ifndef __DENSLIMIT__
#define __DENSLIMIT__ 20.
#endif

#ifndef __XSBUFFER__
#define __XSBUFFER__ 8
#endif

#ifndef __SBUFFER__
#define __SBUFFER__ 32
#endif



/*  ----------------------------------------------------  */
/*  -----------------   Miscellaneous   ----------------  */
/*  ----------------------------------------------------  */


extern const bool _true_;
extern const bool _false_;

extern double _zeroFunc_();
extern double _oneFunc_();



/*  ----------------------------------------------------  */
/*  -----------------   Random Numbers   ---------------  */
/*  ----------------------------------------------------  */


extern gsl_rng *_randomNumberGen_;



#endif // COMMON_H_INCLUDED

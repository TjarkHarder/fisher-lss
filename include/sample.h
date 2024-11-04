#ifndef SAMPLE_H_INCLUDED
#define SAMPLE_H_INCLUDED


#include "common.h"

#include "shape.h"
#include "misc.h"

#include "dat.h"



/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%    STRUCTS    %%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Arguments for the sample functions

    */

    double min;
    double max;

    size_t size;

    double step;

} sample_arg_t;


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Raw sample struct

    */

    sample_arg_t *sampleArg;
    double *array;

} sample_raw_t;


typedef struct
{
    /*

        Shape sample struct (eg line, triangle, ...)

    */

    size_t dim;

    size_t dimLength;
    size_t dimOrientation;

    size_t size; // Size of arrayShape

    size_t sizeParity; // Number of parity invariant shapes
    size_t sizeNoParity; // Number of not parity invariant shapes
    size_t sizeFull; // sizeParity + sizeNoParity

    shape_t **arrayShape;

    sample_raw_t *sampleRawLength;
    sample_raw_t *sampleRawOrientation;

} sample_shape_t;



/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%    EXTERN    %%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */

/*  ----------------------------------------------------  */
/*  ---------------   Extern Functions   ---------------  */
/*  ----------------------------------------------------  */


extern int avr_shape_sample_vol(sample_shape_t *sampleShape);



/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%    SAMPLE ARGUMENTS    %%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  ----------------   New / Free / Cp   ---------------  */
/*  ----------------------------------------------------  */


sample_arg_t *sample_arg_new();

sample_arg_t *sample_arg_free(
    sample_arg_t *sampleArg);

/*  ----------------------------------------------------  */

sample_arg_t *sample_arg_cp(
    sample_arg_t *sampleArg);



/*  ----------------------------------------------------  */
/*  --------------------   Setters   -------------------  */
/*  ----------------------------------------------------  */


int sample_arg_set_max(
    sample_arg_t *sampleArg,
    double max);

int sample_arg_set_min(
    sample_arg_t *sampleArg,
    double min);


int sample_arg_set_size(
    sample_arg_t *sampleArg,
    size_t size);


int sample_arg_set_step(
    sample_arg_t *sampleArg,
    double step);



/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


double sample_arg_get_max(
    sample_arg_t *sampleArg);

double sample_arg_get_min(
    sample_arg_t *sampleArg);


size_t sample_arg_get_size(
    sample_arg_t *sampleArg);


double sample_arg_get_step(
    sample_arg_t *sampleArg);



/*  ----------------------------------------------------  */
/*  --------------------   Compare   -------------------  */
/*  ----------------------------------------------------  */


bool sample_arg_compare_sub(
    sample_arg_t *sampleArg1,
    sample_arg_t *sampleArg2);




/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%    RAW SAMPLING    %%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  ----------------   New / Free / Cp   ---------------  */
/*  ----------------------------------------------------  */


sample_raw_t *sample_raw_new();

sample_raw_t *sample_raw_free(
    sample_raw_t *sampleRaw);

/*  ----------------------------------------------------  */

sample_raw_t *sample_raw_cp(
    sample_raw_t *sampleRaw);


/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


sample_arg_t *sample_raw_get_sample_arg(
    sample_raw_t *sampleRaw);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double sample_raw_get_value(
    sample_raw_t *sampleRaw,
    size_t index);



/*  ----------------------------------------------------  */
/*  ----------------   Linear Sampling   ---------------  */
/*  ----------------------------------------------------  */


sample_raw_t *sample_raw_lin(
    sample_arg_t *sampleArg);




/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%    SHAPE SAMPLING    %%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  ----------------   New / Free / Cp   ---------------  */
/*  ----------------------------------------------------  */


sample_shape_t *sample_shape_new(
    size_t dim,
    size_t size);

sample_shape_t *sample_shape_free(
    sample_shape_t *sampleShape);

/*  ----------------------------------------------------  */

sample_shape_t *sample_shape_cp(
    sample_shape_t *sampleShape);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


bool sample_shape_uniq(
    sample_shape_t *sampleShape,
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


size_t sample_shape_get_size(
    sample_shape_t *sampleShape);

/*  ----------------------------------------------------  */

size_t sample_shape_get_size_parity(
    sample_shape_t *sampleShape);

size_t sample_shape_get_size_nparity(
    sample_shape_t *sampleShape);

size_t sample_shape_get_size_full(
    sample_shape_t *sampleShape);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


shape_t *sample_shape_get_shape(
    sample_shape_t *sampleShape,
    size_t index);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


sample_raw_t *sample_shape_get_sample_raw_length(
    sample_shape_t *sampleShape);

sample_raw_t *sample_shape_get_sample_raw_orientation(
    sample_shape_t *sampleShape);



/*  ----------------------------------------------------  */
/*  ---------------------   Sort   ---------------------  */
/*  ----------------------------------------------------  */


int sample_shape_sort(
    sample_shape_t *sampleShape);

int sample_shape_insort(
    sample_shape_t *sampleShape,
    size_t index,
    bool *success);



/*  ----------------------------------------------------  */
/*  --------------------   Search   --------------------  */
/*  ----------------------------------------------------  */


size_t sample_shape_bsearch(
    sample_shape_t *sampleShape,
    shape_t *shape,
    bool *success);



/*  ----------------------------------------------------  */
/*  -----------------   Line Sampling   ----------------  */
/*  ----------------------------------------------------  */


sample_shape_t *sample_shape_line(
    sample_arg_t *sampleArgLength,
    sample_arg_t *sampleArgOrientation,
    const char *fileName);

sample_shape_t *sample_shape_line_raw(
    sample_raw_t *sampleRawLength,
    sample_raw_t *sampleRawOrientation,
    const char *fileName);



/*  ----------------------------------------------------  */
/*  ---------------   Triangle Sampling   --------------  */
/*  ----------------------------------------------------  */


sample_shape_t *sample_shape_tri(
    sample_arg_t *sampleArgLength,
    sample_arg_t *sampleArgOrientation,
    const char *fileName);

sample_shape_t *sample_shape_tri_raw(
    sample_raw_t *sampleRawLength,
    sample_raw_t *sampleRawOrientation,
    const char *fileName);



/*  ----------------------------------------------------  */
/*  ----------------   Output and Input   --------------  */
/*  ----------------------------------------------------  */


int sample_shape_output(
    const char *fileName,
    sample_shape_t *sampleShape,
    int *prec);

sample_shape_t *sample_shape_input(
    const char *fileName,
    sample_arg_t *sampleArgLength,
    sample_arg_t *sampleArgOrientation,
    size_t *size);




/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


#endif // SAMPLE_H_INCLUDED

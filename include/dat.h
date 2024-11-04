#ifndef DAT_H_INCLUDED
#define DAT_H_INCLUDED

#include "common.h"

#include "misc.h"



/*  ----------------------------------------------------  */
/*  ----------------  Macro Definitions  ---------------  */
/*  ----------------------------------------------------  */


#ifndef HEADER_LEN
#define HEADER_LEN 16384
#endif


#ifndef HEADER_LINE_LEN
#define HEADER_LINE_LEN 128
#endif



/*  ----------------------------------------------------  */
/*  --------------------  Structures  ------------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Data struct

    */

    /* Labels */
    char **xLabels;
    char **yLabels;

    /* Dimensions of the x and y data */
    size_t xDim;
    size_t yDim;

    /* Number of data points */
    size_t size;

    /* Data arrays */
    double *xData;
    double *yData;

} dat_t;



/*  ----------------------------------------------------  */
/*  ------------------   New / Free   ------------------  */
/*  ----------------------------------------------------  */


dat_t *dat_new(
    size_t xDim,
    size_t yDim,
    size_t size);

/*  ----------------------------------------------------  */

dat_t *dat_free(
    dat_t *dat);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


dat_t *dat_cp(
    dat_t *dat);



/*  ----------------------------------------------------  */
/*  ------------   Manipulate the Struct   -------------  */
/*  ----------------------------------------------------  */


dat_t *dat_cat(
    size_t num,
    ...);

/*  ----------------------------------------------------  */

int dat_append_y(
    dat_t *dat,
    char *yLabel);



/*  ----------------------------------------------------  */
/*  --------------------   Setters   -------------------  */
/*  ----------------------------------------------------  */


int dat_set_label(
    dat_t *dat,
    size_t index,
    char type,
    char *label);

/*  ----------------------------------------------------  */

int dat_set_xlabel(
    dat_t *dat,
    size_t index,
    char *label);

int dat_set_ylabel(
    dat_t *dat,
    size_t index,
    char *label);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int dat_set_value(
    dat_t *dat,
    size_t col,
    size_t row,
    char type,
    double value);

/*  ----------------------------------------------------  */

int dat_set_xvalue(
    dat_t *dat,
    size_t col,
    size_t row,
    double value);

int dat_set_yvalue(
    dat_t *dat,
    size_t col,
    size_t row,
    double value);



/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


char *dat_get_label(
    dat_t *dat,
    size_t index,
    char type);

/*  ----------------------------------------------------  */

char *dat_get_xlabel(
    dat_t *dat,
    size_t index);

char *dat_get_ylabel(
    dat_t *dat,
    size_t index);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double dat_get_value(
    dat_t *dat,
    size_t col,
    size_t row,
    char type);

/*  ----------------------------------------------------  */

double dat_get_xvalue(
    dat_t *dat,
    size_t col,
    size_t row);

double dat_get_yvalue(
    dat_t *dat,
    size_t col,
    size_t row);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double *dat_get_array(
    dat_t *dat,
    size_t col,
    char type);

/*  ----------------------------------------------------  */

double *dat_get_xarray(
    dat_t *dat,
    size_t col);

double *dat_get_yarray(
    dat_t *dat,
    size_t col);



/*  ----------------------------------------------------  */
/*  ------------    Input / Output Data    -------------  */
/*  ----------------------------------------------------  */


int dat_output(
    char *fileName,
    dat_t *dat,
    int *precision);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


dat_t *dat_input(
    char *fileName,
    char *yLabel,
    bool fill);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // DAT_H_INCLUDED

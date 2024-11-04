/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTERPOLATE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "interpolate.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INTERPOLATION     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   Initialise Interpolation Functions   ------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


interp_t *interpolate_interp_init(interp_t *(*interpInit)(dat_t *, double (*)(void*, void*, void*)), dat_t *dat, double (*extrapolate)(void*, void*, void*))
{
    /*

        Master function to initialise the interpolation functions for given data.

    */

    interp_t *interp = interpInit(dat, extrapolate);

    return interp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


interp_t *interpolate_interp_init_splinter(dat_t *dat, double (*extrapolate)(void*, void*, void*))
{
    /*

        Setup the n-spline interpolation for input data "dat" using the SPLINTER library by Bjarne Grimstad
        (see: https://github.com/bgrimstad/splinter).

    */


    /* Max degree of B-Spline */

    size_t degree = 3;


    /* Need at least degree + 1 unique y-data points */

    double *yDataUniqBin = calloc(dat -> size, sizeof(double));
    size_t yDataUniqCnt = 0;
    double yDataPnt;
    bool yDataUniq;

    /* Search for unique elements */
    for (size_t n = 0; n < dat -> size; n++)
      {
        yDataUniq = true;
        yDataPnt = dat_get_value(dat, 0, n, 'y');

        for (size_t i = 0; i < yDataUniqCnt; i++)
          {
            /* Not unique if the element is already inside the bin */
            if (yDataUniqBin[i] == yDataPnt)
              {
                yDataUniq = false;
                break;
              }
          }

        /* Add the data point and increment the counter */
        if (yDataUniq)
          {
            yDataUniqBin[yDataUniqCnt] = yDataPnt;

            yDataUniqCnt ++;
          }

        /* Only have to collect degree + 1 unique elements */
        if (yDataUniqCnt == degree + 1)
            break;
      }


    /* Free memory */
    free(yDataUniqBin);

    /* Only one unique element */
    if (yDataUniqCnt == 1)
      {
        /* Spline and acc (both not used) */
        interp_t *interp = malloc(sizeof(interp_t));

        interp -> acc = NULL;
        interp -> spline = NULL;

        /* No extrapolation for a constant */
        interp -> extrapolate = NULL;

        /* Bounds of interpolation (only keep the unique y-data point) */
        interp -> xBounds = NULL;
        interp -> yBounds = malloc(sizeof(double));

        interp -> yBounds[0] = dat_get_value(dat, 0, 0, 'y');

        /* xDimUniq */
        interp -> xDimUniq = dat -> xDim;

        /* xIndUniq */
        interp -> xIndUniq = NULL;

        return interp;
      }

    /* Not enough data points to have spline of degree 3 */
    else if (yDataUniqCnt <= degree)
        degree = yDataUniqCnt - 1;


    /* Need at least degree + 1 unique x-data points in each dimension */

    size_t *degrees = malloc(sizeof(size_t) * dat -> xDim);

    double *xDataUniqBin = calloc(dat -> size, sizeof(double));
    size_t xDataUniqCnt;
    double xDataPnt;
    bool xDataUniq;

    /* Search for unique elements in each dimension */
    for (size_t i = 0; i < dat -> xDim; i++)
      {
        xDataUniqCnt = 0;

        /* Search for unique elements */
        for (size_t n = 0; n < dat -> size; n++)
          {
            xDataUniq = true;
            xDataPnt = dat_get_value(dat, i, n, 'x');

            for (size_t j = 0; j < xDataUniqCnt; j++)
              {
                /* Not unique if the element is already inside the bin */
                if (xDataUniqBin[j] == xDataPnt)
                  {
                    xDataUniq = false;
                    break;
                  }
              }

            /* Add the data point and increment the counter */
            if (xDataUniq)
              {
                xDataUniqBin[xDataUniqCnt] = xDataPnt;

                xDataUniqCnt ++;
              }

            /* Only have to collect degree + 1 unique elements */
            if (xDataUniqCnt == degree + 1)
                break;
          }

        /* Not enough data points -> lower degree */
        if (xDataUniqCnt <= degree)
          {
            degrees[i] = xDataUniqCnt - 1;
          }

        /* Can use degree */
        else
          {
            degrees[i] = degree;
          }
      }

    /* Free memory */
    free(xDataUniqBin);

    /* Remove x dimensions if they only contain one unique data point */
    size_t *xDataIndUniq = malloc(sizeof(size_t) * dat -> xDim);
    size_t xDataDimUniq = 0;

    for (size_t i = 0; i < dat -> xDim; i++)
      {
        /* x-dimension has unique points */
        if (degrees[i])
          {
            xDataIndUniq[xDataDimUniq] = i;
            xDataDimUniq ++;

            continue;
          }

        /* No unique points -> remove the entry */
        for (size_t j = i; j < dat -> xDim - 1; j++)
          {
            degrees[j] = degrees[j + 1];
          }
      }

    /* Reallocate memory */
    xDataIndUniq = realloc(xDataIndUniq, sizeof(size_t) * xDataDimUniq);
    degrees = realloc(degrees, sizeof(size_t) * xDataDimUniq);

    /* No unique data points */
    if (xDataDimUniq == 0)
      {
        printf("Need atleast 2 unique x-data points (in atleast one dimension) for interpolation!\n");
        exit(1);
      }


    /* Create the splinter dataTable */

    /* data as a 1-d array */
    double *totData = malloc(sizeof(double) * (dat -> size * (xDataDimUniq + 1)));
    double *xData = malloc(sizeof(double) * (dat -> size * xDataDimUniq));
    double *yData = malloc(sizeof(double) * dat -> size);

    for (size_t n = 0; n < dat -> size; n++)
      {
        for (size_t i = 0; i < xDataDimUniq; i++)
          {
            totData[(xDataDimUniq + 1) * n + i] = dat_get_value(dat, xDataIndUniq[i], n, 'x');
            xData[xDataDimUniq * n + i] = dat_get_value(dat, xDataIndUniq[i], n, 'x');
          }

        totData[(xDataDimUniq + 1) * n + xDataDimUniq] = dat_get_value(dat, 0, n, 'y');
        yData[n] = dat_get_value(dat, 0, n, 'y');
      }

    /* dataTable */
    splinter_obj_ptr dataTable = splinter_datatable_init();

    /* Add the data array to the dataTable (Splinter version 3.0) */
//    splinter_datatable_add_samples_row_major(dataTable, totData, (int) dat -> size, (int) xDataDimUniq);

    /* Add the data array to the dataTable (Splinter version 4.0 */
    splinter_datatable_add_samples_row_major(dataTable, xData, (int) xDataDimUniq, yData, 1, (int) dat -> size);


    /* Build the B-Spline (Splinter version 3.0) */

    /* bsplineBuilder init */
//    splinter_obj_ptr bsplineBuilder = splinter_bspline_builder_init(dataTable);

    /* Set degrees */
//    splinter_bspline_builder_set_degree(bsplineBuilder, (unsigned int*) degrees, (int) xDataDimUniq);

    /* Get the bspline */
//    splinter_obj_ptr bspline = splinter_bspline_builder_build(bsplineBuilder);


    /* Build the B-Spline (Splinter version 4.0) */

    /* Get the bspline */
    splinter_obj_ptr bspline = splinter_bspline_interpolator(dataTable, (int) degree);


    /* Store the result in the interp_t struct */

    /* Spline and acc (latter only for gsl) */
    interp_t *interp = malloc(sizeof(interp_t));

    interp -> acc = NULL;
    interp -> spline = bspline;
    interp -> extrapolate = extrapolate;

    /* Bounds of interpolation */
    interp -> xBounds = malloc(sizeof(double*) * dat -> xDim);
    interp -> yBounds = malloc(sizeof(double) * 2);

    for (size_t i = 0; i < dat -> xDim; i++)
      {
        interp -> xBounds[i] = malloc(sizeof(double) * 2);

        interp -> xBounds[i][0] = dat_get_value(dat, i, 0, 'x');
        interp -> xBounds[i][1] = dat_get_value(dat, i, 0, 'x');

        for (size_t n = 1; n < dat -> size; n++)
          {
            double xval = dat_get_value(dat, i, n, 'x');

            if (xval < interp -> xBounds[i][0])
                interp -> xBounds[i][0] = xval;

            if (xval > interp -> xBounds[i][1])
                interp -> xBounds[i][1] = xval;
          }
      }

    interp -> yBounds[0] = dat_get_value(dat, 0, 0, 'y');
    interp -> yBounds[1] = dat_get_value(dat, 0, dat -> size - 1, 'y');

    /* xDim + xDimUniq + xIndUniq */
    interp -> xDim = dat -> xDim;
    interp -> xDimUniq = xDataDimUniq;
    interp -> xIndUniq = xDataIndUniq;


    /* Free memory */

    free(degrees);
    free(totData);
    free(xData);
    free(yData);

    splinter_datatable_delete(dataTable);
//    splinter_bspline_builder_delete(bsplineBuilder);


    return interp;
}


/*  ------------------------------------------------------------------------------------------------------  */


interp_t *interpolate_interp_init_gsl(dat_t *dat, double (*extrapolate)(void*, void*, void*))
{
    /*

        Setup the 1-d spline interpolation with gsl.

    */

    /* kDat and pkDat as 1-d arrays */

    double *xData = dat_get_array(dat, 0, 'x');
    double *yData = dat_get_array(dat, 0, 'y');


    /* Interpolate the power spectrum */

    interp_t *interp = malloc(sizeof(interp_t));

    interp -> acc = gsl_interp_accel_alloc();
    interp -> spline = gsl_spline_alloc(gsl_interp_steffen, dat -> size);
    interp -> extrapolate = extrapolate;

    gsl_spline_init(interp -> spline, xData, yData, dat -> size);


    /* Bounds of interpolation */

    interp -> xBounds = malloc(sizeof(double*));
    interp -> xBounds[0] = malloc(sizeof(double) * 2);
    interp -> yBounds = malloc(sizeof(double) * 2);

    interp -> xBounds[0][0] = xData[0];
    interp -> yBounds[0] = yData[0];

    interp -> xBounds[0][1] = xData[dat -> size - 1];
    interp -> yBounds[1] = yData[dat -> size - 1];


    /* xDim */

    interp -> xDimUniq = 1;


    /* Free allocated memory */

    free(xData);
    free(yData);


    return interp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   Evaluate Interpolation Functions   -------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double interpolate_interp_eval(double (*interpEval)(void*, interp_t*, void*), void *values, interp_t *interp, void *params)
{
    /*

        Master interpolation function calls the interpolation function "interpFunc"
        for the input values "values" with additional parameters "params" for the extrapolation function.

    */

    return interpEval(values, interp, params);
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double interpolate_interp_eval_splinter(void *values, interp_t *interp, void *params)
{
    /*

        Evaluate the spline interpolation by splinter (see interpolate_interp_init_splinter).

    */

    double *xValues = (double*) values;

    /* If spline is NULL the interpolation function was fed a constant */
    if (interp -> spline == NULL)
      {
        return interp -> yBounds[0];
      }

    /* Only keep the unique x-dimensions */
    double *xValuesUniq = malloc(sizeof(double) * interp -> xDimUniq);
    bool inBounds = true;

    for (size_t i = 0; i < interp -> xDimUniq; i++)
      {
        xValuesUniq[i] = xValues[interp -> xIndUniq[i]];

        /* Check the bounds */
        if (xValuesUniq[i] < interp -> xBounds[interp -> xIndUniq[i]][0] || xValuesUniq[i] > interp -> xBounds[interp -> xIndUniq[i]][1])
            inBounds = false;
      }

    /* Result of the interpolation */
    double *res, resReturn;

    /* In bounds, no need for extrapolation */
    if (inBounds)
      {
        res = splinter_bspline_eval_row_major((splinter_obj_ptr) interp -> spline, xValuesUniq, (int) interp -> xDimUniq);
        resReturn = *res;
        free(res);
      }

    /* Out of bounds, need to extrapolate */
    else
      {
        if (interp -> extrapolate == NULL)
          {
            for (size_t i = 0; i < interp -> xDimUniq; i++)
              {
                xValuesUniq[i] = xValues[interp -> xIndUniq[i]];

                printf("%e < %e < %e\n", interp -> xBounds[interp -> xIndUniq[i]][0], xValuesUniq[i], interp -> xBounds[interp -> xIndUniq[i]][1]);
              }

            /* Cannot extrapolate */
            printf("Cannot extrapolate if no extrapolation function has been provided.\n");
            exit(1);

            return NAN;
          }

        resReturn = interp -> extrapolate(values, interp, params);
      }

    /* Free memories */
    free(xValuesUniq);

    return resReturn;
}


/*  ------------------------------------------------------------------------------------------------------  */


double interpolate_interp_eval_gsl(void *values, interp_t *interp, void *params)
{
    /*

        Evaluate the interpolation function by gsl (1-d).

    */

    /* Convert input variables */
    double *xValues = (double*) values;

    double res;

    /* x is inside interpolated region */
    if (*xValues >= interp -> xBounds[0][0] && *xValues <= interp -> xBounds[0][1])
        res = gsl_spline_eval(interp -> spline, *xValues, interp -> acc);

    /* Out of bounds */
    else
      {
        if (interp -> extrapolate == NULL)
          {
            printf("Cannot extrapolate if no extrapolation function has been provided.\n");
            exit(1);

            return NAN;
          }

        res = interp -> extrapolate(values, interp, params);
      }

    return res;
}


double interpolate_interp_eval_deriv_gsl(void *values, interp_t *interp, void *params)
{
    /*

        Evaluate the derivative of the interpolation function by gsl (1-d)

    */

    /* Convert input variables */
    double *xValues = (double*) values;

    double res;

    /* x in interpolated region */
    if (*xValues >= interp -> xBounds[0][0] && *xValues <= interp -> xBounds[0][1])
        res = gsl_spline_eval_deriv(interp -> spline, *xValues, interp -> acc);

    /* Out of bounds */
    else
      {
        if (interp -> extrapolate == NULL)
          {
            goto extrapolateErr;
          }

        res = interp -> extrapolate(values, interp, params);
      }

    return res;

    /* Cannot extrapolate */
    extrapolateErr:
    printf("\nCannot extrapolate if no extrapolation function has been provided.\n\n");
    exit(1);

    return 0.;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Free Interpolation Functions   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


interp_t *interpolate_interp_free(interp_t *(*interpFree)(interp_t*), interp_t *interp)
{
    /*

        Master function to free the allocated memory of the interp_t struct.

    */

    interp = interpFree(interp);

    return interp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


interp_t *interpolate_interp_free_splinter(interp_t *interp)
{
    /*

        Free the allocated memory of an interp_t struct for the splinter library (see interpolate_interp_init_splinter).

    */


    /* If interp is NULL simply return NULL as well */
    if (interp == NULL)
        return NULL;

    /* Free interp's contents */
    if (interp -> spline != NULL) splinter_bspline_delete((splinter_obj_ptr) interp -> spline);

    if (interp -> xBounds != NULL)
      {
        for (size_t i = 0; i < interp -> xDim; i++)
          {
            if (interp -> xBounds[i] != NULL)  free(interp -> xBounds[i]);
          }

        free(interp -> xBounds);
      }

    free(interp -> yBounds);
    free(interp -> xIndUniq);

    /* Free interp as well */
    free(interp);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


interp_t *interpolate_interp_free_gsl(interp_t *interp)
{
    /*

        Free the allocated memory of an interp_t struct.

    */


    /* If interp is NULL simply return NULL as well */
    if (interp == NULL)
        return NULL;

    /* Free interp's contents */
    if (interp -> acc != NULL) gsl_interp_accel_free(interp -> acc);
    if (interp -> spline != NULL) gsl_spline_free(interp -> spline);

    if (interp -> xBounds != NULL)
      {
        if (interp -> xBounds[0] != NULL) free(interp -> xBounds[0]);

        free(interp -> xBounds);
      }

    free(interp -> yBounds);

    /* Free interp as well */
    free(interp);

    return NULL;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

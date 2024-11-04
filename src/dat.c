/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     DAT.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "dat.h"


// TODO: Maybe add an additional header to the output function...


/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     DATA STRUCT     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Create / Free Struct   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


dat_t *dat_new(size_t xDim, size_t yDim, size_t size)
{
    /*

        Create a new dat_t struct and allocate memory to everything

    */

    /* Struct */
    dat_t *dat = malloc(sizeof(dat_t));

    /* Dimensions */
    dat -> xDim = xDim;
    dat -> yDim = yDim;
    dat -> size = size;

    /* x */
    if (dat -> xDim == 0)
      {
        dat -> xLabels = NULL;
        dat -> xData = NULL;
      }

    else
      {
        dat -> xLabels = malloc(sizeof(char*) * dat -> xDim);
        dat -> xData = calloc(dat -> xDim * dat -> size, sizeof(double));

        for (size_t i = 0; i < dat -> xDim; i++)
          {
            dat -> xLabels[i] = NULL;
          }
      }

    /* y */
    if (dat -> yDim == 0)
      {
        dat -> yLabels = NULL;
        dat -> yData = NULL;
      }

    else
      {
        dat -> yLabels = malloc(sizeof(char*) * dat -> yDim);
        dat -> yData = calloc(dat -> yDim * dat -> size, sizeof(double));

        for (size_t i = 0; i < dat -> yDim; i++)
          {
            dat -> yLabels[i] = NULL;
          }
      }

    return dat;
}


/*  ------------------------------------------------------------------------------------------------------  */


dat_t *dat_free(dat_t *dat)
{
    /*

        Free a dat_t structure.

    */

    /* If dat is NULL simply return NULL */
    if (dat == NULL)
        return NULL;

    /* Labels */
    for (size_t i = 0; i < dat -> xDim; i++)
      {
        if (dat -> xLabels[i] != NULL) free(dat -> xLabels[i]);
      }

    free(dat -> xLabels);

    for (size_t i = 0; i < dat -> yDim; i++)
      {
        if (dat -> yLabels[i] != NULL) free(dat -> yLabels[i]);
      }

    free(dat -> yLabels);

    /* Data */
    free(dat -> xData);
    free(dat -> yData);

    /* Free dat itself */
    free(dat);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


dat_t *dat_cp(dat_t *dat)
{
    /*

        Copy a dat_t struct

    */

    /* If dat is NULL return NULL as well */
    if (dat == NULL)
        return NULL;

    /* Copy dat */
    dat_t *datCp = dat_new(dat -> xDim, dat -> yDim, dat -> size);

    /* Copy the data arrays */
    for (size_t n = 0; n < datCp -> size; n++)
      {
        for (size_t i = 0; i < datCp -> xDim; i++)
            dat_set_value(datCp, i, n, 'x', dat_get_value(dat, i, n, 'x'));

        for (size_t i = 0; i < datCp -> yDim; i++)
            dat_set_value(datCp, i, n, 'y', dat_get_value(dat, i, n, 'y'));
      }

    /* Copy the labels */
    for (size_t i = 0; i < datCp -> xDim; i++)
        dat_set_label(datCp, i, 'x', dat_get_label(dat, i, 'x'));

    for (size_t i = 0; i < datCp -> yDim; i++)
        dat_set_label(datCp, i, 'y', dat_get_label(dat, i, 'y'));

    return datCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Manipulate the Struct   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


dat_t *dat_cat(size_t num, ...)
{
    /*

        Concatenate several dat structs.

        This function only concatenates the x data ignoring the y data.

    */

    if (num == 0)
        return NULL;

    dat_t *datCat = NULL;

    size_t xDim = 0;
    size_t size = 1;

    va_list ap;

    /* Get the total x dimension and all xSizes */
    va_start(ap, num);

    for (size_t i = 0; i < num; i++)
      {
        dat_t *dat = va_arg(ap, dat_t*);

        if (dat == NULL)
            continue;

        xDim += dat -> xDim;
        size *= dat -> size;
      }


    va_end(ap);

    /* Get the data */
    va_start(ap, num);

    datCat = dat_new(xDim, 0, size);
    xDim = 0;

    size_t subSize1 = datCat -> size;
    size_t subSize2 = subSize1;

    for (size_t i = 0; i < num; i++)
      {
        dat_t *dat = va_arg(ap, dat_t*);

        if (dat == NULL)
            continue;

        subSize2 /= dat -> size;

        for (size_t j = 0; j < dat -> xDim; j++)
          {
            /* Set the label */
            dat_set_label(datCat, xDim, 'x', dat -> xLabels[j]);

            /* Add all the data points to the right location */
            for (size_t n = 0; n < datCat -> size; n++)
                dat_set_value(datCat, xDim, n, 'x', dat_get_value(dat, j, (n % subSize1 - n % subSize2) / subSize2, 'x'));

            xDim += 1;
          }

        subSize1 = subSize2;
      }

    va_end(ap);

    return datCat;
}


/*  ------------------------------------------------------------------------------------------------------  */


int dat_append_y(dat_t *dat, char *yLabel)
{
    /*

        Append a new y column with some label to dat

        Note: No attempt at reordering the data inside the array is made. This should only be used for appending columns
        to a new struct!

    */

    /* Increment the y dimension */
    dat -> yDim += 1;

    /* Reallocate memory to the yLabels and set the label */
    dat -> yLabels = realloc(dat -> yLabels, sizeof(char*) * dat -> yDim);
    dat -> yLabels[dat -> yDim - 1] = NULL;
    dat_set_label(dat, dat -> yDim - 1, 'y', yLabel);

    /* Reallocate memory to the data array */
    dat -> yData = realloc(dat -> yData, sizeof(double) * (dat -> yDim * dat -> size));

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Setters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int dat_set_label(dat_t *dat, size_t index, char type, char *label)
{
    /*

        Set the index' label of given type for dat

    */

    /* Set an x label */
    if (type == 'x')
      {
        dat_set_xlabel(dat, index, label);

        return 0;
      }

    /* Set an y label */
    if (type == 'y')
      {
        dat_set_ylabel(dat, index, label);

        return 0;
      }

    /* Wrong type */
    printf("Expected 'x' or 'y' for type but instead got '%c'!\n", type);
    exit(1);

    return 1;
}


/*  ------------------------------------------------------------------------------------------------------  */


int dat_set_xlabel(dat_t *dat, size_t index, char *label)
{
    /*

        Set the index' x label

    */

    /* Out of bounds */
    if (index >= dat -> xDim)
      {
        printf("Cannot set the x label at index %ld if the dimension is %ld.\n", index, dat -> xDim);
        exit(1);

        return 1;
      }

    /* Unset the label */
    if (label == NULL || strlen(label) == 0)
      {
        free(dat -> xLabels[index]);
        dat -> xLabels[index] = NULL;

        return 0;
      }

    /* Set the label */
    misc_scp(&dat -> xLabels[index], label);

    return 0;
}


int dat_set_ylabel(dat_t *dat, size_t index, char *label)
{
    /*

        Set the index' y label

    */

    /* Out of bounds */
    if (index >= dat -> yDim)
      {
        printf("Cannot set the y label at index %ld if the dimension is %ld.\n", index, dat -> yDim);
        exit(1);

        return 1;
      }

    /* Unset the label */
    if (label == NULL || strlen(label) == 0)
      {
        free(dat -> yLabels[index]);
        dat -> yLabels[index] = NULL;

        return 0;
      }

    /* Set the label */
    misc_scp(&dat -> yLabels[index], label);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int dat_set_value(dat_t *dat, size_t col, size_t row, char type, double value)
{
    /*

        Set a value in dat

    */

    /* Set an x value */
    if (type == 'x')
      {
        dat_set_xvalue(dat, col, row, value);

        return 0;
      }

    /* Set an y value */
    if (type == 'y')
      {
        dat_set_yvalue(dat, col, row, value);

        return 0;
      }

    /* Incorrect type */
    printf("Expected 'x' or 'y' for type but instead got '%c'!\n", type);
    exit(1);

    return 1;
}


/*  ------------------------------------------------------------------------------------------------------  */


int dat_set_xvalue(dat_t *dat, size_t col, size_t row, double value)
{
    /*

        Set an x value in dat

    */

    if (col >= dat -> xDim || row >= dat -> size)
      {
        printf("Expected xDim = %ld > %ld = col and size = %ld > %ld = row!\n", dat -> xDim, col, dat -> size, row);
        exit(1);

        return 1;
      }

    dat -> xData[dat -> size * col + row] = value;

    return 0;
}


int dat_set_yvalue(dat_t *dat, size_t col, size_t row, double value)
{
    /*

        Set an y value in dat

    */

    if (col >= dat -> yDim || row >= dat -> size)
      {
        printf("Expected yDim = %ld > %ld = col and size = %ld > %ld = row!\n", dat -> yDim, col, dat -> size, row);
        exit(1);

        return 1;
      }

    dat -> yData[dat -> size * col + row] = value;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Getters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *dat_get_label(dat_t *dat, size_t index, char type)
{
    /*

        Get the index' label of given type for dat

    */

    /* Get an x label */
    if (type == 'x')
      {
        return dat_get_xlabel(dat, index);
      }

    /* Get an y label */
    if (type == 'y')
      {
        return dat_get_ylabel(dat, index);
      }

    /* Wrong type */
    printf("Expected 'x' or 'y' for type but instead got '%c'!\n", type);
    exit(1);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


char *dat_get_xlabel(dat_t *dat, size_t index)
{
    /*

        Get the index' x label

    */

    /* Out of bounds */
    if (index >= dat -> xDim)
      {
        printf("Cannot get the x label at index %ld if the dimension is %ld.\n", index, dat -> xDim);
        exit(1);

        return NULL;
      }

    /* Return the label */
    return dat -> xLabels[index];
}


char *dat_get_ylabel(dat_t *dat, size_t index)
{
    /*

        Get the index' y label

    */

    /* Out of bounds */
    if (index >= dat -> yDim)
      {
        printf("Cannot get the y label at index %ld if the dimension is %ld.\n", index, dat -> yDim);
        exit(1);

        return NULL;
      }

    /* Return the label */
    return dat -> yLabels[index];
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double dat_get_value(dat_t *dat, size_t col, size_t row, char type)
{
    /*

        Get a value from dat

    */

    /* Get an x value */
    if (type == 'x')
      {
        return dat_get_xvalue(dat, col, row);
      }

    /* Get an y value */
    if (type == 'y')
      {
        return dat_get_yvalue(dat, col, row);
      }

    /* Incorrect type */
    printf("Expected 'x' or 'y' for type but instead got '%c'!\n", type);
    exit(1);

    return NAN;
}


/*  ------------------------------------------------------------------------------------------------------  */


double dat_get_xvalue(dat_t *dat, size_t col, size_t row)
{
    /*

        Get an x value from dat

    */

    if (col >= dat -> xDim || row >= dat -> size)
      {
        printf("Expected xDim = %ld > %ld = col and size = %ld > %ld = row!\n", dat -> xDim, col, dat -> size, row);
        exit(1);

        return NAN;
      }

    return dat -> xData[dat -> size * col + row];
}


double dat_get_yvalue(dat_t *dat, size_t col, size_t row)
{
    /*

        Get an y value from dat

    */

    if (col >= dat -> yDim || row >= dat -> size)
      {
        printf("Expected yDim = %ld > %ld = col and size = %ld > %ld = row!\n", dat -> yDim, col, dat -> size, row);
        exit(1);

        return NAN;
      }

    return dat -> yData[dat -> size * col + row];
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double *dat_get_array(dat_t *dat, size_t col, char type)
{
    /*

        Get a full column from dat

    */

    /* Get an x value */
    if (type == 'x')
      {
        return dat_get_xarray(dat, col);
      }

    /* Get an y value */
    if (type == 'y')
      {
        return dat_get_yarray(dat, col);
      }

    /* Incorrect type */
    printf("Expected 'x' or 'y' for type but instead got '%c'!\n", type);
    exit(1);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


double *dat_get_xarray(dat_t *dat, size_t col)
{
    /*

        Get an x array from dat

    */

    if (col >= dat -> xDim)
      {
        printf("Expected xDim = %ld > %ld = col!\n", dat -> xDim, col);
        exit(1);

        return NULL;
      }

    double *array = malloc(sizeof(double) * dat -> size);

    for (size_t n = 0; n < dat -> size; n++)
        array[n] = dat_get_value(dat, col, n, 'x');

    return array;
}


double *dat_get_yarray(dat_t *dat, size_t col)
{
    /*

        Get an y array from dat

    */

    if (col >= dat -> yDim)
      {
        printf("Expected yDim = %ld > %ld = col!\n", dat -> yDim, col);
        exit(1);

        return NULL;
      }

    double *array = malloc(sizeof(double) * dat -> size);

    for (size_t n = 0; n < dat -> size; n++)
        array[n] = dat_get_value(dat, col, n, 'y');

    return array;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     OUTPUT DATA     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Output Data   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static char *_dat_output_fprintf_get_header(dat_t *dat, size_t offset);
static char *_dat_output_fwrite_get_header(dat_t *dat, size_t offset);

static int _dat_output_fprintf_header(FILE *stream, char *header);
static int _dat_output_fwrite_header(FILE *stream, char *header);

static int _dat_output_fprintf(FILE *stream, dat_t *dat, int precision);
static int _dat_output_fwrite(FILE *stream, dat_t *dat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int dat_output(char *fileName, dat_t *dat, int *precision)
{
    /*

        Output dat to a file.

        If precision is NULL write in binary mode else in text mode with 'precision' number of decimal places.

    */

    /* Open the file */
    FILE *stream = fopen(fileName, "w");

    if (stream == NULL)
      {
        printf("Could not write data to the file '%s'. Make sure the location exists.\n", fileName);

        return 0;
      }

    /* Get the header */
    char *header;

    /* Write the data to file */
    if (precision == NULL)
      {
        header = _dat_output_fwrite_get_header(dat, 0);

        _dat_output_fwrite_header(stream, header);
        _dat_output_fwrite(stream, dat);
      }

    else
      {
        header = _dat_output_fprintf_get_header(dat, 0);

        _dat_output_fprintf_header(stream, header);
        _dat_output_fprintf(stream, dat, *precision);
      }

    /* Free memory */
    free(header);

    /* Close the file */
    fclose(stream);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_output_fprintf_header(FILE *stream, char *header)
{
    /*

        Write the total header to stream in text mode

    */

    /* Write header to file */
    fprintf(stream, "%s", header);

    return 0;
}


static int _dat_output_fwrite_header(FILE *stream, char *header)
{
    /*

        Write the total header to stream in binary mode

    */

    /* Write header to file */
    fwrite(header, strlen(header) + 1, 1, stream);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_output_fprintf(FILE *stream, dat_t *dat, int precision)
{
    /*

        Write dat to stream in text mode

    */

    for (size_t n = 0; n < dat -> size; n++)
      {
        /* x data */
        for (size_t i = 0; i < dat -> xDim; i++)
          {
            fprintf(stream, "%.*e", precision, dat_get_value(dat, i, n, 'x'));

            /* Seperate the values by white space */
            if (i < dat -> xDim - 1)
                fprintf(stream, " ");
          }

        /* Seperate x data and y data by white space */
        fprintf(stream, " ");

        /* y data */
        for (size_t i = 0; i < dat -> yDim; i++)
          {
            fprintf(stream, "%.*e", precision, dat_get_value(dat, i, n, 'y'));

            /* Seperate the values by white space */
            if (i < dat -> yDim - 1)
                fprintf(stream, " ");
          }

        /* new line */
        fprintf(stream, "\n");
      }

    return 0;
}


static int _dat_output_fwrite(FILE *stream, dat_t *dat)
{
    /*

        Write dat to stream in binary mode

    */

    for (size_t n = 0; n < dat -> size; n++)
      {
        /* x data */
        for (size_t i = 0; i < dat -> xDim; i++)
          {
            double val = dat_get_value(dat, i, n, 'x');
            fwrite(&val, sizeof(double), 1, stream);
          }

        /* y data */
        for (size_t i = 0; i < dat -> yDim; i++)
          {
            double val = dat_get_value(dat, i, n, 'y');
            fwrite(&val, sizeof(double), 1, stream);
          }
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------   Construct a Header for a Text File   --------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _dat_output_fprintf_get_header_guide(char *header);
static int _dat_output_fprintf_get_header_line(dat_t *dat, char type, size_t index, char *header);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static char *_dat_output_fprintf_get_header(dat_t *dat, size_t offset)
{
    /*

        Get the header for dat for a text file

    */

    char *header = malloc(HEADER_LEN);
    header[0] = '\0';

    /* Add guide about displayed information */
    _dat_output_fprintf_get_header_guide(header);

    /* Line break */
    strcat(header, "\n");

    /* Data size */
    _dat_output_fprintf_get_header_line(dat, 's', 0, header);

    /* Line break */
    strcat(header, "\n");

    /* Line offset */
    _dat_output_fprintf_get_header_line(dat, 'l', offset, header);

    /* Line break */
    strcat(header, "\n");

    if (dat -> xDim != 0)
      {
        /* x data */
        for (size_t i = 0; i < dat -> xDim; i++)
          {
            _dat_output_fprintf_get_header_line(dat, 'x', i, header);
          }

        /* Line break */
        strcat(header, "\n");
      }

    if (dat -> yDim != 0)
      {
        /* y data */
        for (size_t i = 0; i < dat -> yDim; i++)
          {
            _dat_output_fprintf_get_header_line(dat, 'y', i, header);
          }

        /* Line break */
        strcat(header, "\n");
      }

    /* Seperate the header from the data */
    char *line = malloc(HEADER_LINE_LEN + 1);

    for (size_t i = 0; i < HEADER_LINE_LEN; i++)
        line[i] = '=';

    line[HEADER_LINE_LEN] = '\0';
    strcat(header, line);

    strcat(header, "\n\n");

    /* Free memory */
    free(line);

    /* Offset must equal number of '\n' characters in the header */
    size_t count = misc_scnch(header, '\n');

    if (offset < count)
      {
        char *headerNew = _dat_output_fprintf_get_header(dat, count);

        free(header);
        header = headerNew;
      }

    return header;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_output_fprintf_get_header_guide(char *header)
{
    /*

        Add a guide about the displayed information

    */

    /* Initialise header */
    header[0] = '\0';

    /* Long line to enclose guide (+2 for \n and \0) */
    char *encloseLine = malloc(HEADER_LINE_LEN + 2);

    for (size_t i = 0; i < HEADER_LINE_LEN; i++)
        encloseLine[i] = '#';

    encloseLine[HEADER_LINE_LEN] = '\n';
    encloseLine[HEADER_LINE_LEN + 1] = '\0';

    /* Add the enclosure line to the header */
    strcat(header, encloseLine);


    /* Add the guide */

    char *line = malloc(HEADER_LINE_LEN + 2);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*  Information about the data is displayed as:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     data size : <size>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     line offset : <offset>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     x data in column <column> : <label>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     y data in column <column> : <label>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    /* Add the enclosure line to the header */
    strcat(header, encloseLine);


    /* Free memory */

    free(line);
    free(encloseLine);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_output_fprintf_get_header_line(dat_t *dat, char type, size_t index, char *header)
{
    /*

        Get the header line for dat located in some column.

    */

    char *headerLine = malloc(HEADER_LINE_LEN);

    /* data size */
    if (type == 's')
      {
        snprintf(headerLine, HEADER_LINE_LEN, "data size : %ld\n", dat -> size);
      }

    /* data size */
    if (type == 'l')
      {
        snprintf(headerLine, HEADER_LINE_LEN, "line offset : %ld\n", index + 1);
      }

    /* x data */
    if (type == 'x')
      {
        /* Without label */
        if (dat -> xLabels[index] == NULL)
          {
            snprintf(headerLine, HEADER_LINE_LEN, "x data in column %ld : --\n", index);
          }

        /* With label */
        else
          {
            snprintf(headerLine, HEADER_LINE_LEN, "x data in column %ld : %s\n", index, dat -> xLabels[index]);
          }
      }

    /* y data */
    else if (type == 'y')
      {
        if (dat -> yLabels[index] == NULL)
          {
            snprintf(headerLine, HEADER_LINE_LEN, "y data in column %ld : --\n", index + dat -> xDim);
          }

        else
          {
            snprintf(headerLine, HEADER_LINE_LEN, "y data in column %ld : %s\n", index + dat -> xDim, dat -> yLabels[index]);
          }
      }


    /* Add the line to the header */
    strcat(header, headerLine);

    /* Free memory */
    free(headerLine);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------   Construct a Header for a Binary File   -------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _dat_output_fwrite_get_header_guide(char *header);
static int _dat_output_fwrite_get_header_line(dat_t *dat, char type, size_t index, size_t offset, char *header);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static char *_dat_output_fwrite_get_header(dat_t *dat, size_t offset)
{
    /*

        Get the header for dat for a binary file

    */

    char *header = malloc(HEADER_LEN);
    header[0] = '\0';

    /* Add guide about displayed information */
    _dat_output_fwrite_get_header_guide(header);

    /* Line break */
    strcat(header, "\n");

    /* data size */
    _dat_output_fwrite_get_header_line(dat, 's', 0, 0, header);

    /* Line break */
    strcat(header, "\n");

    /* Byte offset */
    _dat_output_fwrite_get_header_line(dat, 'b', 0, offset, header);

    /* Line break */
    strcat(header, "\n");

    if (dat -> xDim != 0)
      {
        /* x data */
        for (size_t i = 0; i < dat -> xDim; i++)
          {
            _dat_output_fwrite_get_header_line(dat, 'x', i, i * sizeof(double), header);
          }

        /* Line break */
        strcat(header, "\n");
      }

    if (dat -> yDim != 0)
      {
        /* y data */
        for (size_t i = 0; i < dat -> yDim; i++)
          {
            _dat_output_fwrite_get_header_line(dat, 'y', i, (i + dat -> xDim) * sizeof(double), header);
          }

        /* Line break */
        strcat(header, "\n");
      }

    /* Seperate the header from the data */
    char *line = malloc(HEADER_LINE_LEN + 1);

    for (size_t i = 0; i < HEADER_LINE_LEN; i++)
        line[i] = '=';

    line[HEADER_LINE_LEN] = '\0';
    strcat(header, line);

    strcat(header, "\n\n");

    /* Free memory */
    free(line);

    /* Size of the header must equal string length of the header */
    if (offset < strlen(header))
      {
        char *headerNew = _dat_output_fwrite_get_header(dat, strlen(header));

        free(header);
        header = headerNew;
      }


    return header;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_output_fwrite_get_header_guide(char *header)
{
    /*

        Add a guide about the displayed information

    */

    /* Initialise header */
    header[0] = '\0';

    /* Long line to enclose guide (+2 for \n and \0) */
    char *encloseLine = malloc(HEADER_LINE_LEN + 2);

    for (size_t i = 0; i < HEADER_LINE_LEN; i++)
        encloseLine[i] = '#';

    encloseLine[HEADER_LINE_LEN] = '\n';
    encloseLine[HEADER_LINE_LEN + 1] = '\0';

    /* Add the enclosure line to the header */
    strcat(header, encloseLine);


    /* Add the guide */

    char *line = malloc(HEADER_LINE_LEN + 2);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*  Information about the data is displayed as:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     data size : <size>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     byte offset : <offset>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     x data at bytes <byte> : <label>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                   :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                   :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     y data at bytes <byte> : <label>\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                   :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                   :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*  where <byte> is to be understood as an offset to the first element in a 'row'. For example:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*   ===========================\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     data size : 250\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     byte offset : 800\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     x data at bytes 0 : x1\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     x data at bytes 8 : x2\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     x data at bytes 16 : x3\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     y data at bytes 24 : y1\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     y data at bytes 32 : y2\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*   ===========================\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*  would have a total of 250 data points (i.e. 250 points of x1, x2, x3, y1 and y2) starting at the offset byte 800.\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*  The x1 data points are at byte positions (800 + ) 0, 40, 80, ..., x2 at (800 + ) 8, 48, 88, ... and so on.\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    /* Add the enclosure line to the header */
    strcat(header, encloseLine);


    /* Free memory */

    free(line);
    free(encloseLine);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_output_fwrite_get_header_line(dat_t *dat, char type, size_t index, size_t offset, char *header)
{
    /*

        Get the header line for dat located in some column.

    */

    char *headerLine = malloc(HEADER_LINE_LEN);

    /* data size */
    if (type == 's')
      {
        snprintf(headerLine, HEADER_LINE_LEN, "data size : %ld\n", dat -> size);
      }

    /* byte offset */
    if (type == 'b')
      {
        snprintf(headerLine, HEADER_LINE_LEN, "byte offset : %ld\n", offset + 1);
      }

    /* x data */
    if (type == 'x')
      {
        /* Without label */
        if (dat -> xLabels[index] == NULL)
          {
            snprintf(headerLine, HEADER_LINE_LEN, "x data at bytes %ld : --\n", offset);
          }

        /* With label */
        else
          {
            snprintf(headerLine, HEADER_LINE_LEN, "x data at bytes %ld : %s\n", offset, dat -> xLabels[index]);
          }
      }

    /* y data */
    else if (type == 'y')
      {
        if (dat -> yLabels[index] == NULL)
          {
            snprintf(headerLine, HEADER_LINE_LEN, "y data at bytes %ld : --\n", offset);
          }

        else
          {
            snprintf(headerLine, HEADER_LINE_LEN, "y data at bytes %ld : %s\n", offset, dat -> yLabels[index]);
          }
      }

    /* Add the line to the header */
    strcat(header, headerLine);

    /* Free memory */
    free(headerLine);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INPUT DATA     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Input: Local Structures   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


typedef struct
{
    /*

        Entries in a header line

    */

    char type;

    char *loc;

    char *label;
    char *size;

    char *offset;

} _dat_header_line_t;


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static _dat_header_line_t *_dat_header_line_free(_dat_header_line_t *headerLine)
{
    /*

        Free a _dat_header_line_t struct

    */

    /* If headerLine is NULL simply return NULL */
    if (headerLine == NULL)
        return NULL;

    /* Free headerLine's contents */
    free(headerLine -> loc);

    free(headerLine -> label);
    free(headerLine -> size);

    free(headerLine -> offset);

    /* Free headerLine itself */
    free(headerLine);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Input Data from File   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static dat_t *_dat_input_fscanf_data(char *fileName, dat_t *headerDat, size_t *yCol, size_t yDimTot, size_t offset);
static dat_t *_dat_input_fread_data(char *fileName, dat_t *headerDat, size_t *yByte, size_t yDimTot, size_t offset);

static dat_t *_dat_input_get_header(char *fileName, char *yLabel, size_t **yLoc, size_t *yDimTot, size_t *offset, bool *binary);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


dat_t *dat_input(char *fileName, char *yLabel, bool fill)
{
    /*

        Read data from file. If yLabel is not NULL get the data with this label.

        Only read the actual data if fill is true otherwise just get the data labels, dimensions and sizes from the header.

    */

    size_t *yLoc = NULL;
    size_t yDimTot = 0;

    size_t offset = 0;
    bool binary = false;

    /* Get the data struct from the header */
    dat_t *headerDat = _dat_input_get_header(fileName, yLabel, &yLoc, &yDimTot, &offset, &binary);

    /* No data in the file */
    if (headerDat == NULL)
      {
        if (yLabel == NULL)
            printf("Could not find any data in the file '%s'.\n", fileName);

        else
            printf("Could not find the specified data with label '%s' in the file '%s'.\n", yLabel, fileName);

        exit(1);

        return NULL;
      }

    /* Fill in the data */
    if (fill)
      {
        dat_t *dat;

        if (binary)
            dat = _dat_input_fread_data(fileName, headerDat, yLoc, yDimTot, offset);

        else
            dat = _dat_input_fscanf_data(fileName, headerDat, yLoc, yDimTot, offset);

        headerDat = dat_free(headerDat);
        headerDat = dat;
      }

    /* Free memory */
    free(yLoc);

    return headerDat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _dat_input_fscanf_skip_header(FILE *stream, size_t offset);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static dat_t *_dat_input_fscanf_data(char *fileName, dat_t *headerDat, size_t *yCol, size_t yDimTot, size_t offset)
{
    /*

        Read data from a text file

    */

    /* Open the file */
    FILE *stream = fopen(fileName, "r");

    /* Skip to the end of the header */
    _dat_input_fscanf_skip_header(stream, offset);

    /* New struct */
    dat_t *dat = dat_new(headerDat -> xDim, headerDat -> yDim, headerDat -> size);

    for (size_t i = 0; i < dat -> xDim; i++)
        dat_set_xlabel(dat, i, headerDat -> xLabels[i]);

    for (size_t i = 0; i < dat -> yDim; i++)
        dat_set_ylabel(dat, i, headerDat -> yLabels[i]);

    /* Read in the data */
    for (size_t n = 0; n < dat -> size; n++)
      {
        double val;
        size_t col = 0;

        /* x data */
        for (size_t i = 0; i < dat -> xDim; i++)
          {
            int numScan = fscanf(stream, "%lf", &val);
            (void) numScan; // Ignore number of variables scanned and assume everything went smoothly :)

            dat_set_value(dat, i, n, 'x', val);
            col++;
          }

        size_t yColLoc = 0;

        /* y data */
        for (size_t i = 0; i < yDimTot; i++)
          {
            int numScan = fscanf(stream, "%lf", &val);
            (void) numScan; // Ignore number of variables scanned and assume everything went smoothly :)

            /* test if the column is in yCol */
            if (yColLoc < dat -> yDim && col == yCol[yColLoc])
              {
                dat_set_value(dat, yColLoc, n, 'y', val);
                yColLoc++;
              }

            col++;
          }
      }

    /* Close the file */
    fclose(stream);

    return dat;
}


static int _dat_input_fscanf_skip_header(FILE *stream, size_t offset)
{
    /*

        Get to the end of the header

    */

    char c;
    size_t lineCount = 0;

    do
      {
        c = (char) fgetc(stream);

        /* Line ends with '\n' */
        if (c == '\n')
          {
            lineCount++;
          }

        /* Header ends (-1 since next line already contains data) */
        if (lineCount == offset - 1)
          {
            break;
          }
      }

    while (c != EOF);

    if (c == EOF)
      {
        printf("File ended before finding the end of the header!\n");
        exit(1);

        return 1;
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static dat_t *_dat_input_fread_data(char *fileName, dat_t *headerDat, size_t *yByte, size_t yDimTot, size_t offset)
{
    /*

        Read data from a binary file

    */

    /* Open the file */
    FILE *stream = fopen(fileName, "r");

    /* Skip to the end of the header */
    fseek(stream, (long int) offset, SEEK_SET);

    /* New struct */
    dat_t *dat = dat_new(headerDat -> xDim, headerDat -> yDim, headerDat -> size);

    for (size_t i = 0; i < dat -> xDim; i++)
        dat_set_xlabel(dat, i, headerDat -> xLabels[i]);

    for (size_t i = 0; i < dat -> yDim; i++)
        dat_set_ylabel(dat, i, headerDat -> yLabels[i]);

    /* Read in the data */
    for (size_t n = 0; n < dat -> size; n++)
      {
        double val;
        size_t byte = 0;

        /* x data */
        for (size_t i = 0; i < dat -> xDim; i++)
          {
            size_t numRead = fread(&val, sizeof(double), 1, stream);
            (void) numRead; // Ignore number of variables read and assume everything went smoothly :)

            dat_set_value(dat, i, n, 'x', val);
            byte += sizeof(double);
          }

        size_t yByteLoc = 0;

        /* y data */
        for (size_t i = 0; i < yDimTot; i++)
          {
            size_t numRead = fread(&val, sizeof(double), 1, stream);
            (void) numRead; // Ignore number of variables read and assume everything went smoothly :)

            /* Test if the byte is in yByte (+ make sure to be in bounds) */
            if (yByteLoc < dat -> yDim && byte == yByte[yByteLoc])
              {
                dat_set_value(dat, yByteLoc, n, 'y', val);
                yByteLoc++;
              }

            byte += sizeof(double);
          }
      }

    /* Close the file */
    fclose(stream);

    return dat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Data Struct from Header   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _dat_input_header_skip_guide(FILE *stream);

static int _dat_input_get_header_dat(FILE *stream, dat_t *headerDat, char *yLabel, size_t **yLoc, size_t *yDimTot, size_t *offset, bool *binary);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static dat_t *_dat_input_get_header(char *fileName, char *yLabel, size_t **yLoc, size_t *yDimTot, size_t *offset, bool *binary)
{
    /*

        Read a file (fileName) and return all of the data or for only one yLabel.

    */

    /* Open the file */
    FILE *stream = fopen(fileName, "r");

    if (stream == NULL)
      {
        printf("Could not open the file '%s'\n", fileName);
        exit(1);

        return NULL;
      }

    /* Skip the guide */
    _dat_input_header_skip_guide(stream);

    /* Read data */
    dat_t *headerDat = dat_new(0, 0, 0);

    int headerStatus;

    do
      {
        headerStatus = _dat_input_get_header_dat(stream, headerDat, yLabel, yLoc, yDimTot, offset, binary);
      }

    while (headerStatus == 0);

    /* Close the file */
    fclose(stream);

    /* Header ended and correct label was not found */
    if (headerDat -> size == 0 || headerDat -> yDim == 0)
      {
        headerDat = dat_free(headerDat);
      }

    return headerDat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_input_header_skip_guide(FILE *stream)
{
    /*

        Get to the end of the guide in the header

    */

    char c;

    size_t guideBarrierLen = 0;
    bool guideStart = false;
    bool guideEnd = false;

    do
      {
        c = (char) fgetc(stream);

        /* Increment barrier length */
        if (c == '#')
          {
            guideBarrierLen ++;
          }

        /* Reset barrier length */
        else
          {
            guideBarrierLen = 0;
          }

        /* Guide starts / ends */
        if (guideBarrierLen == HEADER_LINE_LEN)
          {
            if (!guideStart)
                guideStart = true;

            else
                guideEnd = true;
          }
      }

    while (!guideEnd && c != EOF);

    if (c == EOF)
      {
        printf("File ended before finding the beginning of the data information in the header!\n");
        exit(1);

        return 1;
      }

    /* Need to read in two more characters (\n) */
    c = (char) fgetc(stream);
    c = (char) fgetc(stream);


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static _dat_header_line_t *_dat_input_get_header_line_dissect(FILE *stream, bool *end);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _dat_input_get_header_dat(FILE *stream, dat_t *headerDat, char *yLabel, size_t **yLoc, size_t *yDimTot, size_t *offset, bool *binary)
{
    /*

        Get the data struct of the header.

        Returns 1 if the header ended and 0 else

    */

    /* Header line */
    bool headerEnd = false;
    _dat_header_line_t *headerLine = _dat_input_get_header_line_dissect(stream, &headerEnd);

    /* Header has ended */
    if (headerEnd)
        return 1;

    /* Empty line */
    if (headerLine == NULL)
        return 0;

    /* Data size */
    if (headerLine -> type == 's')
      {
        headerDat -> size = (size_t) atoi(headerLine -> size);
      }

    /* Byte offset */
    else if (headerLine -> type == 'b')
      {
        *binary = true;
        *offset = (size_t) atoi(headerLine -> offset);
      }

    /* Line offset */
    else if (headerLine -> type == 'l')
      {
        *binary = false;
        *offset = (size_t) atoi(headerLine -> offset);
      }

    /* Add another x dimension */
    else if (headerLine -> type == 'x')
      {
        headerDat -> xDim += 1;

        /* Set the label (must first reallocate memory and set the label to NULL) */
        headerDat -> xLabels = realloc(headerDat -> xLabels, sizeof(char*) * headerDat -> xDim);
        headerDat -> xLabels[headerDat -> xDim - 1] = NULL;
        dat_set_xlabel(headerDat, headerDat -> xDim - 1, headerLine -> label);
      }

    /* Add another y dimension if the label is correct */
    else if (headerLine -> type == 'y')
      {
        *yDimTot += 1;

        /* Correct label of the y data */
        bool correctLabel = true;

        /* yLabel is provided but y data does not have the same label */
        if (yLabel != NULL && strcmp(headerLine -> label, yLabel))
          {
            correctLabel = false;
          }

        /* Add another y dimension */
        if (correctLabel)
          {
            headerDat -> yDim += 1;

            /* Set the label (must first reallocate memory and set the label to NULL) */
            headerDat -> yLabels = realloc(headerDat -> yLabels, sizeof(char*) * headerDat -> yDim);
            headerDat -> yLabels[headerDat -> yDim - 1] = NULL;
            dat_set_ylabel(headerDat, headerDat -> yDim - 1, headerLine -> label);

            /* Add the location */
            *yLoc = realloc(*yLoc, sizeof(size_t) * headerDat -> yDim);
            (*yLoc)[headerDat -> yDim - 1] = (size_t) atoi(headerLine -> loc);
          }
      }

    /* Free memory */
    headerLine = _dat_header_line_free(headerLine);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static char *_dat_input_get_header_line(FILE *stream);

static bool _dat_input_get_header_line_empty(char *line);
static bool _dat_input_get_header_line_end(char *line);

static char _dat_input_get_header_line_type(char *line);

static char *_dat_input_get_header_line_size(char *line);
static char *_dat_input_get_header_line_offset(char *line);

static char *_dat_input_get_header_line_loc(char *line);
static char *_dat_input_get_header_line_label(char *line);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

static _dat_header_line_t *_dat_input_get_header_line_dissect(FILE *stream, bool *end)
{
    /*

        Read a line from a dat header and store the relevant information in a _dat_header_line_t struct.

    */

    /* Read the next line */
    char *line = _dat_input_get_header_line(stream);

    /* Empty line (copy and pasting stuff in file automatically adds \r at the end of every line... interesting stuff :)) */
    if (_dat_input_get_header_line_empty(line))
      {
        free(line);

        return NULL;
      }

    /* Header has ended */
    if (_dat_input_get_header_line_end(line))
      {
        *end = true;
        free(line);

        return NULL;
      }

    /* headerLine struct */
    _dat_header_line_t *headerLine = malloc(sizeof(_dat_header_line_t));

    headerLine -> loc = NULL;

    headerLine -> size = NULL;
    headerLine -> label = NULL;

    headerLine -> offset = NULL;

    /* Type */
    headerLine -> type = _dat_input_get_header_line_type(line);

    /* Size */
    if (headerLine -> type == 's')
      {
        headerLine -> size = _dat_input_get_header_line_size(line);
      }

    else if (headerLine -> type == 'b' || headerLine -> type == 'l')
      {
        headerLine -> offset = _dat_input_get_header_line_offset(line);
      }

    /* Column and label */
    else
      {
        headerLine -> loc = _dat_input_get_header_line_loc(line);
        headerLine -> label = _dat_input_get_header_line_label(line);
      }

    free(line);

    return headerLine;
}


/*  ------------------------------------------------------------------------------------------------------  */


static char *_dat_input_get_header_line(FILE *stream)
{
    /*

        Read a line from stream

    */

    size_t lineLoc = 0;
    char *line = malloc(HEADER_LINE_LEN + 2);
    char c;

    do
      {
        c = (char) fgetc(stream);

        line[lineLoc] = c;
        lineLoc++;
      }

    while(c != '\n' && c != EOF && lineLoc < HEADER_LINE_LEN);

    line[lineLoc] = '\0';

    return line;
}


/*  ------------------------------------------------------------------------------------------------------  */


static bool _dat_input_get_header_line_empty(char *line)
{
    /*

        Check if the line is empty (ie only ' ', '\r\n' or '\n')

    */

    char c;
    size_t index = 0;

    bool empty = true;

    do
      {
        c = line[index];
        index++;

        if (!isspace(c) && c != 0)
          {
            empty = false;
            break;
          }
      }

    while (c != '\0');

    return empty;
}


static bool _dat_input_get_header_line_end(char *line)
{
    /*

        Check if the line signals the end of the header (128 x '=')

    */

    size_t headerBarrierLen = 0;
    size_t lineLen = strlen(line);

    /* Line is smaller than max header line length -> can return false */
    if (lineLen < HEADER_LINE_LEN)
        return false;

    /* Check how many '='s are in line */
    for (size_t i = 0; i < lineLen; i++)
      {
        if (line[i] == '=')
            headerBarrierLen ++;
      }

    /* If entire line is made up of '='s return true */
    if (headerBarrierLen == HEADER_LINE_LEN)
        return true;

    return false;
}


/*  ------------------------------------------------------------------------------------------------------  */


static char _dat_input_get_header_line_type(char *line)
{
    /*

        Check which type of information the line contains

    */

    if (misc_sin_rm(line, "data size"))
        return 's';

    if (misc_sin_rm(line, "byte offset"))
        return 'b';

    if (misc_sin_rm(line, "line offset"))
        return 'l';

    if (misc_sin_rm(line, "x data in column") || misc_sin_rm(line, "x data at bytes"))
        return 'x';

    if (misc_sin_rm(line, "y data in column") || misc_sin_rm(line, "y data at bytes"))
        return 'y';

    printf("Expected line to be of shape '%s', '%s', '%s', '%s', %s', '%s' or '%s' but instead got '%s'!\n"
           , "data size : <>", "byte offset : <>", "line offset : <>", "x data at bytes <> : <>", "x data in column <> : <>", "y data at bytes <> : <>","y data in column <> : <>", line);

    exit(1);

    return '\0';
}


static char *_dat_input_get_header_line_size(char *line)
{
    /*

        Get the size of the data

    */

    char c;
    size_t index = 0;

    char *size = NULL;
    size_t sizeLen = 0;

    do
      {
        c = line[index++];

        /* Only need to store the digits */
        if (isdigit(c))
          {
            sizeLen++;

            size = realloc(size, sizeLen + 1);
            size[sizeLen - 1] = c;
            size[sizeLen] = '\0';
          }
      }

    while (c != '\0');

    return size;
}


static char *_dat_input_get_header_line_offset(char *line)
{
    /*

        Get the offset of the data

    */

    char c;
    size_t index = 0;

    char *offset = NULL;
    size_t offsetLen = 0;

    do
      {
        c = line[index++];

        /* Only need to store the digits */
        if (isdigit(c))
          {
            offsetLen++;

            offset = realloc(offset, offsetLen + 1);
            offset[offsetLen - 1] = c;
            offset[offsetLen] = '\0';
          }
      }

    while (c != '\0');

    return offset;
}


static char *_dat_input_get_header_line_loc(char *line)
{
    /*

        Get the column of the data

    */

    char c;
    size_t index = 0;

    char *loc = NULL;
    size_t locLen = 0;

    do
      {
        c = line[index++];

        /* Break before label begins */
        if (c == ':')
          {
            index++; // Skip white space after ':'
            break;
          }

        /* Only need to store the digits */
        if (isdigit(c))
          {
            locLen++;

            loc = realloc(loc, locLen + 1);
            loc[locLen - 1] = c;
            loc[locLen] = '\0';
          }
      }

    while (c != '\0');

    /* Shorten the line */
    size_t lineLen = strlen(line);

    for (size_t i = 0; i < lineLen - index + 1; i++)
      {
        line[i] = line[i + index];
      }

    return loc;
}


static char *_dat_input_get_header_line_label(char *line)
{
    /*

        Get the label of the data

    */

    char c;
    size_t index = 0;

    char *label = NULL;
    size_t labelLen = 0;

    do
      {
        c = line[index++];

        if (c == '\r' || c == '\n')
            break;

        labelLen++;

        label = realloc(label, labelLen + 1);
        label[labelLen - 1] = c;
        label[labelLen] = '\0';
      }

    while (c != '\0');

    return label;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

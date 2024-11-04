/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SAMPLE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "sample.h"


/**  Output directory and file extension  **/

static const char *_outDir = "/output/sample/";





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SAMPLE ARGUMENTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   New / Free / Cp Functions   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_arg_t *sample_arg_new()
{
    /*

        Create a new sample_arg_t struct

    */

    sample_arg_t *sampleArg = malloc(sizeof(sample_arg_t));

    sampleArg -> size = 0;
    sampleArg -> step = 0.;

    sampleArg -> min = 0.;
    sampleArg -> max = 0.;

    return sampleArg;
}


sample_arg_t *sample_arg_free(sample_arg_t *sampleArg)
{
    /*

        Free sampleArg

    */

    /* Check for NULL */
    if (sampleArg == NULL)
        return NULL;

    free(sampleArg);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_arg_t *sample_arg_cp(sample_arg_t *sampleArg)
{
    /*

        Copy sampleArg

    */

    /* Check for NULL */
    if (sampleArg == NULL)
        return NULL;

    sample_arg_t *sampleArgCp = sample_arg_new();

    sampleArgCp -> size = sampleArg -> size;
    sampleArgCp -> step = sampleArg -> step;

    sampleArgCp -> min = sampleArg -> min;
    sampleArgCp -> max = sampleArg -> max;

    return sampleArgCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Setters   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int sample_arg_set_max(sample_arg_t *sampleArg, double max)
{
    /*

        Set the upper bound of the sampling

    */

    /* size and step must be zero */
    if (sampleArg -> size != 0 || sampleArg -> step != 0.)
      {
        printf("Resetting 'sample_arg_t' struct! Must set 'min' and 'max' before setting 'size' or 'step'.\n");

        sampleArg -> size = 0;
        sampleArg -> step = 0.;
      }

    sampleArg -> max = max;

    return 0;
}


int sample_arg_set_min(sample_arg_t *sampleArg, double min)
{
    /*

        Set the lower bound of the sampling

    */

    /* size and step must be zero */
    if (sampleArg -> size != 0 || sampleArg -> step != 0.)
      {
        printf("Resetting 'sample_arg_t' struct! Must set 'min' and 'max' before setting 'size' or 'step'.\n");

        sampleArg -> size = 0;
        sampleArg -> step = 0.;
      }

    sampleArg -> min = min;

    return 0;
}


int sample_arg_set_size(sample_arg_t *sampleArg, size_t size)
{
    /*

        Set the size of the sampling and deduce the step

    */

    sampleArg -> size = size;

    sampleArg -> step = (sampleArg -> size == 1) ? 0. : (sampleArg -> max - sampleArg -> min) / (double) (sampleArg -> size - 1);

    return 0;
}


int sample_arg_set_step(sample_arg_t *sampleArg, double step)
{
    /*

        Set the step of the sampling and deduce the size. Must also alter max!

    */

    sampleArg -> step = fabs(step); // Must be positive

    double size = 1. + (sampleArg -> max - sampleArg -> min) / sampleArg -> step;
    double sizeFloor = floor(size);

    sampleArg -> size = (fabs(size - sizeFloor - 1.) < __ABSTOL__) ? (size_t) sizeFloor + 1
                                                                   : (size_t) sizeFloor;

    double max = sampleArg -> min + (double) (sampleArg -> size - 1) * sampleArg -> step; // Must change max in case (max - min) / step is not an integer

    sampleArg -> max = (fabs(sampleArg -> max - max) < __ABSTOL__) ? sampleArg -> max
                                                                   : max;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Getters   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double sample_arg_get_max(sample_arg_t *sampleArg)
{
    /*

        Get the upper bound of the sampling

    */

    return sampleArg -> max;
}


double sample_arg_get_min(sample_arg_t *sampleArg)
{
    /*

        Get the lower bound of the sampling

    */

    return sampleArg -> min;
}


size_t sample_arg_get_size(sample_arg_t *sampleArg)
{
    /*

        Get the size of the sampling

    */

    return sampleArg -> size;
}


double sample_arg_get_step(sample_arg_t *sampleArg)
{
    /*

        Get the step of the sampling

    */

    return sampleArg -> step;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Compare   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool sample_arg_compare_sub(sample_arg_t *sampleArg1, sample_arg_t *sampleArg2)
{
    /*

        Compare sampleArg1 and sampleArg2 and deduce whether one is a subset of the other

    */

    /* Must have the same step */
    if (sampleArg1 -> step != sampleArg2 -> step)
        return false;

    /* Offset must be integer */
    double offset = fabs(sampleArg1 -> min - sampleArg2 -> min) / sampleArg1 -> step;

    if (fabs(floor(offset) - offset) > __ABSTOL__ && fabs(ceil(offset) - offset) > __ABSTOL__)
        return false;

    return true;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     RAW SAMPLING     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   New / Free / Cp Functions   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_raw_t *sample_raw_new()
{
    /*

        Create a new sample_raw_t struct

    */

    sample_raw_t *sampleRaw = malloc(sizeof(sample_raw_t));

    sampleRaw -> sampleArg = NULL;
    sampleRaw -> array = NULL;

    return sampleRaw;
}


sample_raw_t *sample_raw_free(sample_raw_t *sampleRaw)
{
    /*

        Free sampleRaw

    */

    /* Check for NULL */
    if (sampleRaw == NULL)
        return NULL;

    /* Free contents */
    sampleRaw -> sampleArg = sample_arg_free(sampleRaw -> sampleArg);
    free(sampleRaw -> array);

    /* Free sampleRaw */
    free(sampleRaw);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_raw_t *sample_raw_cp(sample_raw_t *sampleRaw)
{
    /*

        Copy sampleRaw

    */

    /* Check for NULL */
    if (sampleRaw == NULL)
        return NULL;

    sample_raw_t *sampleRawCp = sample_raw_new();

    sampleRawCp -> sampleArg = sample_arg_cp(sampleRaw -> sampleArg);

    if (sampleRaw -> array == NULL)
      {
        sampleRawCp -> array = NULL;
        return sampleRawCp;
      }

    sampleRawCp -> array = malloc(sizeof(double) * sampleRaw -> sampleArg -> size);

    for (size_t i = 0; i < sampleRaw -> sampleArg -> size; i++)
        sampleRawCp -> array[i] = sampleRaw -> array[i];

    return sampleRawCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------------   Getters   --------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_arg_t *sample_raw_get_sample_arg(sample_raw_t *sampleRaw)
{
    /*

        Get the sample_arg_t struct of sampleRaw

    */

    if (sampleRaw == NULL)
        return NULL;

    return sampleRaw -> sampleArg;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double sample_raw_get_value(sample_raw_t *sampleRaw, size_t index)
{
    /*

        Get the index'th element of the sample

    */

    if (sampleRaw == NULL || sampleRaw -> sampleArg -> size < index)
        return NAN;

    return sampleRaw -> array[index];
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------------   Sample Functions   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_raw_t *sample_raw_lin(sample_arg_t *sampleArg)
{
    /*

        Sample size number of points evenly between min and max

    */

    /* Linear sample */
    sample_raw_t *sampleRaw = sample_raw_new();

    sampleRaw -> sampleArg = sample_arg_cp(sampleArg);

    sampleRaw -> array = malloc(sizeof(double) * sampleArg -> size);

    /* Only one element */
    if (sampleArg -> size == 1 || sampleArg -> step == 0.)
      {
        sampleRaw -> array[0] = sampleArg -> min;

        return sampleRaw;
      }

    /* Add elements by filling in elements from both sides */
    sampleRaw -> array[0] = sampleArg -> min;
    sampleRaw -> array[sampleArg -> size - 1] = sampleArg -> max;

    for (size_t i = 0; i < (sampleArg -> size - sampleArg -> size % 2) / 2 - 1; i++)
      {
        sampleRaw -> array[i + 1] = sampleRaw -> array[i] + sampleArg -> step;
        sampleRaw -> array[sampleArg -> size - 2 - i] = sampleRaw -> array[sampleArg -> size - 1 - i] - sampleArg -> step;
      }

    /* If size is odd need to add one more element in the center */
    if (sampleArg -> size % 2 == 1)
        sampleRaw -> array[(sampleArg -> size - sampleArg -> size % 2) / 2] = ( sampleRaw -> array[(sampleArg -> size - sampleArg -> size % 2) / 2 + 1]
                                                                     + sampleRaw -> array[(sampleArg -> size - sampleArg -> size % 2) / 2 - 1] ) / 2.;

    return sampleRaw;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SHAPE SAMPLING     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   New / Free / Cp Functions   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_shape_t *sample_shape_new(size_t dim, size_t size)
{
    /*

        Create a new sample_shape_t struct of given size and dimension (eg 2 for line, 3 for triangle, ...)

    */

    sample_shape_t *sampleShape = malloc(sizeof(sample_shape_t));

    sampleShape -> dim = dim;
    sampleShape -> dimLength = (dim < 3) ? dim - 1 : dim; // For lines only have dim - 1 = 1 independent lengths; higher order shapes have dim independent lengths
    sampleShape -> dimOrientation = (dim < 4) ? dim - 1 : dim; // For lines and triangles have dim - 1 independent orientations; higher order shapes have dim independent orientations

    sampleShape -> size = size;

    sampleShape -> sizeParity = 0;
    sampleShape -> sizeNoParity = 0;
    sampleShape -> sizeFull = 0;

    sampleShape -> arrayShape = malloc(sizeof(shape_t*) * size);

    for (size_t n = 0; n < size; n++)
      {
        sampleShape -> arrayShape[n] = shape_new(dim);
      }

    sampleShape -> sampleRawLength = NULL;
    sampleShape -> sampleRawOrientation = NULL;

    return sampleShape;
}


sample_shape_t *sample_shape_free(sample_shape_t *sampleShape)
{
    /*

        Free sampleShape

    */

    /* Check for NULL */
    if (sampleShape == NULL)
        return NULL;

    for (size_t n = 0; n < sampleShape -> size; n++)
      {
        sampleShape -> arrayShape[n] = shape_free(sampleShape -> arrayShape[n]);
      }

    free(sampleShape -> arrayShape);

    sampleShape -> sampleRawLength = sample_raw_free(sampleShape -> sampleRawLength);
    sampleShape -> sampleRawOrientation = sample_raw_free(sampleShape -> sampleRawOrientation);

    free(sampleShape);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_shape_t *sample_shape_cp(sample_shape_t *sampleShape)
{
    /*

        Copy sampleShape

    */

    /* Check for NULL */
    if (sampleShape == NULL)
        return NULL;

    sample_shape_t *sampleShapeCp = sample_shape_new(sampleShape -> dim, sampleShape -> size);

    sampleShapeCp -> size = sampleShape -> size;

    sampleShapeCp -> sizeParity = sampleShape -> sizeParity;
    sampleShapeCp -> sizeNoParity = sampleShape -> sizeNoParity;
    sampleShapeCp -> sizeFull = sampleShape -> sizeFull;

    sampleShapeCp -> dimLength = sampleShape -> dimLength;
    sampleShapeCp -> dimOrientation = sampleShape -> dimOrientation;

    for (size_t n = 0; n < sampleShapeCp -> size; n++)
      {
        sampleShapeCp -> arrayShape[n] = shape_free(sampleShapeCp -> arrayShape[n]); // Need to free first, due to allocation in new
        sampleShapeCp -> arrayShape[n] = shape_cp(sampleShape -> arrayShape[n]);
      }

    sampleShapeCp -> sampleRawLength = sample_raw_cp(sampleShape -> sampleRawLength);
    sampleShapeCp -> sampleRawOrientation = sample_raw_cp(sampleShape -> sampleRawOrientation);

    return sampleShapeCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool sample_shape_uniq(sample_shape_t *sampleShape, shape_t *shape)
{
    /*

        Test whether a shape is already included in the sampleShape array (i.e. via a permutation)

    */

    for (size_t n = 0; n < sampleShape -> size; n++)
      {
        int comp = shape_comp(sampleShape -> arrayShape[n], shape);

        /* Shapes can be the same or related by parity */
        if (comp != 0)
          {
            /* Keep the shape where the first orientation angle is positive */
            if (sampleShape -> arrayShape[n] -> orientation[0] < 0.)
                shape_parity(sampleShape -> arrayShape[n]);

            return false;
          }
      }

    return true;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Getters   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t sample_shape_get_size(sample_shape_t *sampleShape)
{
    /*

        Get the number of unique shapes in sample

    */

    if (sampleShape == NULL)
        return 0;

    return sampleShape -> size;
}


/*  ------------------------------------------------------------------------------------------------------  */


size_t sample_shape_get_size_parity(sample_shape_t *sampleShape)
{
    /*

        Get the number of parity invariant shapes in sample

    */

    if (sampleShape == NULL)
        return 0;

    return sampleShape -> sizeParity;
}


size_t sample_shape_get_size_nparity(sample_shape_t *sampleShape)
{
    /*

        Get the number of shapes that are not parity invariant in sample

    */

    if (sampleShape == NULL)
        return 0;

    return sampleShape -> sizeNoParity;
}


size_t sample_shape_get_size_full(sample_shape_t *sampleShape)
{
    /*

        Get the full number of shapes in sample

    */

    if (sampleShape == NULL)
        return 0;

    return sampleShape -> sizeFull;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


shape_t *sample_shape_get_shape(sample_shape_t *sampleShape, size_t index)
{
    /*

        Get the index'th shape in sampleShape

    */

    if (sampleShape == NULL || sampleShape -> size < index)
        return NULL;

    return sampleShape -> arrayShape[index];
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_raw_t *sample_shape_get_sample_raw_length(sample_shape_t *sampleShape)
{
    /*

        Get the raw sample of the lengths of the vertexes

    */

    if (sampleShape == NULL)
        return NULL;

    return sampleShape -> sampleRawLength;
}


sample_raw_t *sample_shape_get_sample_raw_orientation(sample_shape_t *sampleShape)
{
    /*

        Get the raw sample of the orientations of the vertexes

    */

    if (sampleShape == NULL)
        return NULL;

    return sampleShape -> sampleRawOrientation;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Sort Sample of Shapes   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _sort_comp(shape_t *shape1, shape_t *shape2);


static int _sort_kern(shape_t **shapes, int low, int high, bool (*comp)(shape_t*, shape_t*));
static int _sort_partition(shape_t **shapes, int low, int high, bool (*comp)(shape_t*, shape_t*));

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int sample_shape_sort(sample_shape_t *sampleShape)
{
    /*

        Sort the shapes in sampleShape such that the first sizeParity shapes are parity invariant and the last
        size - sizeParity are not parity invariant. Within each block the shapes are then further sorted by
        their respective volumes, from small to large

    */

    /* No parity invariant shapes or no not parity invariant shapes */
    if (sampleShape -> sizeParity == 0 || sampleShape -> sizeNoParity == 0)
        goto sortVertices;


    /**  Sort the Shapes by Parity Invariance  **/

    shape_t **shape1 = NULL;
    shape_t **shape2 = NULL;

    bool search1 = true;
    bool search2 = true;

    for (size_t n1 = 0, n2 = sampleShape -> size - 1; n1 < sampleShape -> sizeParity;)
      {
        /* Get next shape if searching */
        if (search1) shape1 = &sampleShape -> arrayShape[n1];
        if (search2) shape2 = &sampleShape -> arrayShape[n2];

        /* Set search flags if shapes are not parity invariant */
        if (!(*shape1) -> parity) search1 = false;
        if ((*shape2) -> parity) search2 = false;

        /* Swap shapes if both search flags are false */
        if (!search1 && !search2)
          {
            shape_swap(shape1, shape2);

            search1 = true;
            search2 = true;
          }

        /* Increment or decrement indexes */
        if (search1) n1++;
        if (search2) n2--;
      }

    sortVertices:


    /**  Sort the Parity Invariant Shapes  **/

      {
        /* Sort the Vertices */
        for (size_t n = 0; n < sampleShape -> sizeParity; n++)
            shape_sort(sampleShape -> arrayShape[n]);

        /* Sort the Shapes */
        _sort_kern(sampleShape -> arrayShape, 0, (int) sampleShape -> sizeParity - 1, _sort_comp);
      }


    /**  Sort the Parity Non-Invariant Shapes  **/

      {
        /* Sort the Vertices */
        for (size_t n = sampleShape -> sizeParity; n < sampleShape -> size; n++)
            shape_sort(sampleShape -> arrayShape[n]);

        /* Sort the Shapes */
        _sort_kern(sampleShape -> arrayShape, (int) sampleShape -> sizeParity, (int) sampleShape -> size - 1, _sort_comp);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static bool _sort_comp(shape_t *shape1, shape_t *shape2)
{
    /*

        Compare shape1 with shape2

    */

    /* Place shape with lowest vertices' lengths to the left */
    for (size_t i = 0; i < shape1 -> dim; i++)
      {
        double length1 = shape_qget_vertex_length(shape1, i);
        double length2 = shape_qget_vertex_length(shape2, i);

        if (length1 == length2)
            continue;

        return (length1 - length2) <= 0.;
      }

    /* Place shape with  */
    for (size_t i = 0; i < shape1 -> dim; i++)
      {
        double orientation1 = shape_qget_vertex_orientation(shape1, i);
        double orientation2 = shape_qget_vertex_orientation(shape2, i);

        if (orientation1 == orientation2)
            continue;

        if (orientation1 == 0.)
            return true;

        if (orientation2 == 0.)
            return false;

        return (orientation1 * orientation2 > 0.) ? orientation1 <= orientation2 : orientation1 > orientation2;
      }

    return true;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _sort_kern(shape_t **shapes, int low, int high, bool (*comp)(shape_t*, shape_t*))
{
    /*

        Kernel of the quicksort function

    */

    /* Sort as long as pivot is larger than low */
    if (low < high)
      {
        /* Place last pivot, and get new pivot */
        int pivot = _sort_partition(shapes, low, high, comp);

        /* Sort to the left of the previous pivot */
        _sort_kern(shapes, low, pivot - 1, comp);

        /* Sort to the right of the last pivot */
        _sort_kern(shapes, pivot + 1, high, comp);
      }

    return 0;
}


static int _sort_partition(shape_t **shapes, int low, int high, bool (*comp)(shape_t*, shape_t*))
{
    /*

        Partition the vertices

    */

    /* Push i'th vertex to the right */
    int i = low - 1;

    for (int j = low; j < high; j++)
      {
        /* Push j'th vertex to the left if comp is true */
        if (comp(shapes[j], shapes[high]))
          {
            shape_swap(&shapes[++i], &shapes[j]);
          }
      }

    /* Swap i'th vertex with pivot's vertex */
    shape_swap(&shapes[i + 1], &shapes[high]);

    return i + 1;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int sample_shape_insort(sample_shape_t *sampleShape, size_t index, bool *success)
{
    /*

        Given an ordered sampleShape, sort a shape at sampleShape -> arrayShape[index] into the remaining sample.

        Set success to true if the shape was sorted or false if it already exists

    */

    /* Success variable */
    bool success_;

    /* Get the index of the value */
    size_t index_ = sample_shape_bsearch(sampleShape, sampleShape -> arrayShape[index], &success_);

    /* Value is already in the array */
    if (success_)
      {
        if (success != NULL) *success = false;

        return 0;
      }


    /* Must insert the shape */

    for (int i = (int) sampleShape -> size; i > (int) index_; i--)
      {
        /* Push the i'th shape forwards */
        shape_swap(&sampleShape -> arrayShape[i], &sampleShape -> arrayShape[i - 1]);
      }

    /* Insert shape */
    shape_swap(&sampleShape -> arrayShape[index_], &sampleShape -> arrayShape[index]);

    /* Increment the size */
    sampleShape -> size += 1;

    /* Set success to true */
    if (success != NULL) *success = true;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Binary Search   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _bsearch_bounds(sample_shape_t *sampleShape, int *low, int *high, shape_t *shape);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t sample_shape_bsearch(sample_shape_t *sampleShape, shape_t *shape, bool *success)
{
    /*

        Perform a binary search for shape in sorted shapeShape of given size.

        This function returns the index of the found shape and sets success to true, or if not found
        returns the index at which the shape would have been expected with success being set to false.

    */

    /* No memory allocated */
    if (sampleShape == NULL || sampleShape -> size == 0)
      {
        if (success != NULL) *success = false;

        return 0;
      }

    /* Index of shapes in sampleShape */
    size_t index;

    /* Bounds of binary search */
    int low = 0;
    int high = (int) sampleShape -> size - 1;

    /* Success of binary search */
    bool success_;

    do
      {
        /* Check the new bounds */
        success_ = _bsearch_bounds(sampleShape, &low, &high, shape);
      }

    while (!success_ && low <= high);

    /* Store outcome of search in success */
    if (success != NULL) *success = success_;

    /* Index of element is low */
    index = (size_t) low;

    return index;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static bool _bsearch_bounds(sample_shape_t *sampleShape, int *low, int *high, shape_t *shape)
{
    /*

        Look up the mid point in sampleShape -> arrayShape between positions low and high and check if shape is equal, higher or lower
        and change low and high accordingly.

    */

    /* Index of the mid point */
    int mid = (*high - *low) / 2 + *low;

    /* Shapes are the same */
    if (shape_ecomp(sampleShape -> arrayShape[mid], shape))
      {
        *high = mid;
        *low = mid;

        return true;
      }

    /* Figure out if shape is to the right or to the left of mid point */
    if (_sort_comp(sampleShape -> arrayShape[mid], shape))
      {
        *low = mid + 1;
      }

    else
      {
        *high = mid - 1;
      }

    return false;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Line Sampling   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_shape_t *sample_shape_line(sample_arg_t *sampleArgLength, sample_arg_t *sampleArgOrientation, const char *fileName)
{
    /*

        Sample lines of different lengths and orientations

    */

    /* Must have orientation in (-1, 1) */
    if (sampleArgOrientation -> min < -1. || sampleArgOrientation -> max > 1. || sampleArgOrientation -> min > sampleArgOrientation -> max)
      {
        printf("Must have -1 < min < max < 1 for sample args!");
        exit(1);

        return NULL;
      }

    /* Must have symmetric orientation args */
    if (fabs(sampleArgOrientation -> min) != sampleArgOrientation -> max && (sampleArgOrientation -> min < 0 && sampleArgOrientation -> max > 0))
      {
        printf("Must have min = -max for sample args if min < 0 and max > 0.\n");
        exit(1);

        return NULL;
      }

    /* Sample length */
    sample_raw_t *sampleRawLength = sample_raw_lin(sampleArgLength);

    /* Sample orientation */
    sample_raw_t *sampleRawOrientation = sample_raw_lin(sampleArgOrientation);

    /* Sample line */
    size_t dim = 2;
    size_t size = sampleRawLength -> sampleArg -> size * sampleRawOrientation -> sampleArg -> size;

    sample_shape_t *sampleShapeLine = sample_shape_new(dim, size + 1);

    sampleShapeLine -> sampleRawLength = sampleRawLength;
    sampleShapeLine -> sampleRawOrientation = sampleRawOrientation;

    /* Set the entries */
    sampleShapeLine -> size = 0;

    bool success;

    size_t indexLength;
    size_t indexOrientation;

    size_t index;

    for (size_t i = 0; i < size; i++)
      {
        /* Decompose the index */
        index = i;

        indexOrientation = index % sampleRawOrientation -> sampleArg -> size;
        index = (index - indexOrientation) / sampleRawOrientation -> sampleArg -> size;

        indexLength = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength) / sampleRawLength -> sampleArg -> size;

        /* Current shape */
        shape_t *shape = sampleShapeLine -> arrayShape[sampleShapeLine -> size + 1];

        /* Variables */
        double l1 = sampleRawLength -> array[indexLength];
        double l2 = sampleRawLength -> array[indexLength];

        double a12 = -1.;

        double o1 = sampleRawOrientation -> array[indexOrientation];
        double o2 = - sampleRawOrientation -> array[indexOrientation];

        /* Set the lengths */
        if (!shape_set_vertex_length(shape, 0, l1)) continue;
        if (!shape_set_vertex_length(shape, 1, l2)) continue;

        /* Set the angle */
        if (!shape_set_vertex_angle(shape, 0, 1, a12)) continue;

        /* Set the orientations */
        if (!shape_set_vertex_orientation(shape, 0, o1)) continue;
        if (!shape_set_vertex_orientation(shape, 1, o2)) continue;

        /* Check the shape */
        if (!shape_check(shape)) continue;

        /* Sort the shape */
        shape_sort(shape);

        /* Insert and sort */
        sample_shape_insort(sampleShapeLine, sampleShapeLine -> size + 1, &success);

        /* Set the parity */
        if (success)
          {
            shape_set_parity(shape);

            sampleShapeLine -> sizeParity += (size_t) (shape -> parity);
            sampleShapeLine -> sizeNoParity += (size_t) (!shape -> parity);
            sampleShapeLine -> sizeFull += 1 + (size_t) (!shape -> parity);
          }
      }

    /* Reallocate memory */
    for (size_t n = sampleShapeLine -> size; n < size + 1; n++)
      {
        sampleShapeLine -> arrayShape[n] = shape_free(sampleShapeLine -> arrayShape[n]);
      }

    sampleShapeLine -> arrayShape = realloc(sampleShapeLine -> arrayShape, sizeof(shape_t*) * sampleShapeLine -> size);

    /* Set the bin-volume */
    avr_shape_sample_vol(sampleShapeLine);

    /* Output the sample */
    if (fileName != NULL)
      {
        sample_shape_output(fileName, sampleShapeLine, NULL);
      }

    return sampleShapeLine;
}


/*  ------------------------------------------------------------------------------------------------------  */


sample_shape_t *sample_shape_line_raw(sample_raw_t *sampleRawLength, sample_raw_t *sampleRawOrientation, const char *fileName)
{
    /*

        Sample lines of different lengths and orientations directly from a sampling of the lengths and orientations

    */

    /* Sample line */
    size_t dim = 2;
    size_t size = sampleRawLength -> sampleArg -> size * sampleRawOrientation -> sampleArg -> size;

    sample_shape_t *sampleShapeLine = sample_shape_new(dim, size + 1);

    sampleShapeLine -> sampleRawLength = sample_raw_cp(sampleRawLength);
    sampleShapeLine -> sampleRawOrientation = sample_raw_cp(sampleRawOrientation);

    /* Set the entries */
    sampleShapeLine -> size = 0;

    bool success;

    size_t indexLength;
    size_t indexOrientation;

    size_t index;

    for (size_t i = 0; i < size; i++)
      {
        /* Decompose the index */
        index = i;

        indexOrientation = index % sampleRawOrientation -> sampleArg -> size;
        index = (index - indexOrientation) / sampleRawOrientation -> sampleArg -> size;

        indexLength = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength) / sampleRawLength -> sampleArg -> size;

        /* Current shape */
        shape_t *shape = sampleShapeLine -> arrayShape[sampleShapeLine -> size + 1];

        /* Variables */
        double l1 = sampleRawLength -> array[indexLength];
        double l2 = sampleRawLength -> array[indexLength];

        double a12 = -1.;

        double o1 = sampleRawOrientation -> array[indexOrientation];
        double o2 = - sampleRawOrientation -> array[indexOrientation];

        /* Set the lengths */
        if (!shape_set_vertex_length(shape, 0, l1)) continue;
        if (!shape_set_vertex_length(shape, 1, l2)) continue;

        /* Set the angle */
        if (!shape_set_vertex_angle(shape, 0, 1, a12)) continue;

        /* Set the orientations */
        if (!shape_set_vertex_orientation(shape, 0, o1)) continue;
        if (!shape_set_vertex_orientation(shape, 1, o2)) continue;

        /* Check the shape */
        if (!shape_check(shape)) continue;

        /* Sort the shape */
        shape_sort(shape);

        /* Insert and sort */
        sample_shape_insort(sampleShapeLine, sampleShapeLine -> size + 1, &success);

        /* Set the parity */
        if (success)
          {
            shape_set_parity(shape);

            sampleShapeLine -> sizeParity += (size_t) (shape -> parity);
            sampleShapeLine -> sizeNoParity += (size_t) (!shape -> parity);
            sampleShapeLine -> sizeFull += 1 + (size_t) (!shape -> parity);
          }
      }

    /* Reallocate memory */
    for (size_t n = sampleShapeLine -> size; n < size + 1; n++)
      {
        sampleShapeLine -> arrayShape[n] = shape_free(sampleShapeLine -> arrayShape[n]);
      }

    sampleShapeLine -> arrayShape = realloc(sampleShapeLine -> arrayShape, sizeof(shape_t*) * sampleShapeLine -> size);

    /* Set the bin-volume */
    avr_shape_sample_vol(sampleShapeLine);

    /* Output the sample */
    if (fileName != NULL)
      {
        sample_shape_output(fileName, sampleShapeLine, NULL);
      }

    return sampleShapeLine;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Line Sampling   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


sample_shape_t *sample_shape_tri(sample_arg_t *sampleArgLength, sample_arg_t *sampleArgOrientation, const char *fileName)
{
    /*

        Sample triangles of different lengths and orientations

    */

    /* Must have orientation in (-1, 1) */
    if (sampleArgOrientation -> min < -1. || sampleArgOrientation -> max > 1. || sampleArgOrientation -> min > sampleArgOrientation -> max)
      {
        printf("Must have -1 < min < max < 1 for sample args!");
        exit(1);

        return NULL;
      }

    /* Must have symmetric orientation args TODO: Keep this? */
    if (fabs(sampleArgOrientation -> min) != sampleArgOrientation -> max && (sampleArgOrientation -> min < 0. && sampleArgOrientation -> min > 0.))
      {
        printf("Must have min = -max for sample args if min < 0 and max > 0.\n");
        exit(1);

        return NULL;
      }

    /* Sample length */
    sample_raw_t *sampleRawLength = sample_raw_lin(sampleArgLength);

    /* Sample orientation */
    sample_raw_t *sampleRawOrientation = sample_raw_lin(sampleArgOrientation);

    /* Sample triangle */
    size_t dim = 3;
    size_t size = sampleRawLength -> sampleArg -> size * sampleRawLength -> sampleArg -> size * sampleRawLength -> sampleArg -> size * sampleRawOrientation -> sampleArg -> size * sampleRawOrientation -> sampleArg -> size;

    sample_shape_t *sampleShapeTri = sample_shape_new(dim, size + 1);

    sampleShapeTri -> sampleRawLength = sampleRawLength;
    sampleShapeTri -> sampleRawOrientation = sampleRawOrientation;

    /* Set the entries */
    sampleShapeTri -> size = 0;

    bool success;

    size_t indexLength1;
    size_t indexLength2;
    size_t indexLength3;

    size_t indexOrientation1;
    size_t indexOrientation2;

    size_t index;

    for (size_t i = 0; i < size; i++)
      {
        /* Decompose the index */
        index = i;

        indexOrientation2 = index % sampleRawOrientation -> sampleArg -> size;
        index = (index - indexOrientation2) / sampleRawOrientation -> sampleArg -> size;

        indexOrientation1 = index % sampleRawOrientation -> sampleArg -> size;
        index = (index - indexOrientation1) / sampleRawOrientation -> sampleArg -> size;

        indexLength3 = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength3) / sampleRawLength -> sampleArg -> size;

        indexLength2 = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength2) / sampleRawLength -> sampleArg -> size;

        indexLength1 = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength1) / sampleRawLength -> sampleArg -> size;

        /* Current shape */
        shape_t *shape = sampleShapeTri -> arrayShape[sampleShapeTri -> size + 1];

        /* Variables */
        double l1 = sampleRawLength -> array[indexLength1];
        double l2 = sampleRawLength -> array[indexLength2];
        double l3 = sampleRawLength -> array[indexLength3];

        double a12 = (l3*l3 - l1*l1 - l2*l2) / (2. * l1 * l2);
        double a13 = (l2*l2 - l1*l1 - l3*l3) / (2. * l1 * l3);
        double a23 = (l1*l1 - l2*l2 - l3*l3) / (2. * l2 * l3);

        double o1 = sampleRawOrientation -> array[indexOrientation1];
        double o2 = sampleRawOrientation -> array[indexOrientation2];
        double o3 = - (o1*l1 + o2*l2) / l3;


        /* Set the lengths */
        if (!shape_set_vertex_length(shape, 0, l1)) continue;
        if (!shape_set_vertex_length(shape, 1, l2)) continue;
        if (!shape_set_vertex_length(shape, 2, l3)) continue;

        /* Set the angles */
        if (!shape_set_vertex_angle(shape, 0, 1, a12)) continue;
        if (!shape_set_vertex_angle(shape, 0, 2, a13)) continue;
        if (!shape_set_vertex_angle(shape, 1, 2, a23)) continue;

        /* Set the orientations */
        if (!shape_set_vertex_orientation(shape, 0, o1)) continue;
        if (!shape_set_vertex_orientation(shape, 1, o2)) continue;
        if (!shape_set_vertex_orientation(shape, 2, o3)) continue;

        /* Check the third orientation angle */
        if (o3 < sampleRawOrientation -> sampleArg -> min || o3 > sampleRawOrientation -> sampleArg -> max) continue;

        /* Check the shape */
        if (!shape_check(shape)) continue;

        /* Sort the shape */
        shape_sort(shape);

        /* Insert and sort */
        sample_shape_insort(sampleShapeTri, sampleShapeTri -> size + 1, &success);

        /* Set the parity */
        if (success)
          {
            shape_set_parity(shape);

            sampleShapeTri -> sizeParity += (size_t) (shape -> parity);
            sampleShapeTri -> sizeNoParity += (size_t) (!shape -> parity);
            sampleShapeTri -> sizeFull += 1 + (size_t) (!shape -> parity);
          }
      }

    /* Reallocate memory */
    for (size_t n = sampleShapeTri -> size; n < size + 1; n++)
      {
        sampleShapeTri -> arrayShape[n] = shape_free(sampleShapeTri -> arrayShape[n]);
      }

    sampleShapeTri -> arrayShape = realloc(sampleShapeTri -> arrayShape, sizeof(shape_t*) * sampleShapeTri -> size);

    /* Set the bin-volume */
    avr_shape_sample_vol(sampleShapeTri);

    /* Output the sample */
    if (fileName != NULL)
      {
        sample_shape_output(fileName, sampleShapeTri, NULL);
      }

    return sampleShapeTri;
}


/*  ------------------------------------------------------------------------------------------------------  */


sample_shape_t *sample_shape_tri_raw(sample_raw_t *sampleRawLength, sample_raw_t *sampleRawOrientation, const char *fileName)
{
    /*

        Sample triangles of different lengths and orientations directly from a sampling of the lengths and orientations

    */

    /* Sample triangle */
    size_t dim = 3;
    size_t size = sampleRawLength -> sampleArg -> size * sampleRawLength -> sampleArg -> size * sampleRawLength -> sampleArg -> size * sampleRawOrientation -> sampleArg -> size * sampleRawOrientation -> sampleArg -> size;

    sample_shape_t *sampleShapeTri = sample_shape_new(dim, size + 1);

    sampleShapeTri -> sampleRawLength = sample_raw_cp(sampleRawLength);
    sampleShapeTri -> sampleRawOrientation = sample_raw_cp(sampleRawOrientation);

    /* Set the entries */
    sampleShapeTri -> size = 0;

    bool success;

    size_t indexLength1;
    size_t indexLength2;
    size_t indexLength3;

    size_t indexOrientation1;
    size_t indexOrientation2;

    size_t index;

    for (size_t i = 0; i < size; i++)
      {
        /* Decompose the index */
        index = i;

        indexOrientation2 = index % sampleRawOrientation -> sampleArg -> size;
        index = (index - indexOrientation2) / sampleRawOrientation -> sampleArg -> size;

        indexOrientation1 = index % sampleRawOrientation -> sampleArg -> size;
        index = (index - indexOrientation1) / sampleRawOrientation -> sampleArg -> size;

        indexLength3 = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength3) / sampleRawLength -> sampleArg -> size;

        indexLength2 = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength2) / sampleRawLength -> sampleArg -> size;

        indexLength1 = index % sampleRawLength -> sampleArg -> size;
        index = (index - indexLength1) / sampleRawLength -> sampleArg -> size;

        /* Current shape */
        shape_t *shape = sampleShapeTri -> arrayShape[sampleShapeTri -> size + 1];

        /* Variables */
        double l1 = sampleRawLength -> array[indexLength1];
        double l2 = sampleRawLength -> array[indexLength2];
        double l3 = sampleRawLength -> array[indexLength3];

        double a12 = (l3*l3 - l1*l1 - l2*l2) / (2. * l1 * l2);
        double a13 = (l2*l2 - l1*l1 - l3*l3) / (2. * l1 * l3);
        double a23 = (l1*l1 - l2*l2 - l3*l3) / (2. * l2 * l3);

        double o1 = sampleRawOrientation -> array[indexOrientation1];
        double o2 = sampleRawOrientation -> array[indexOrientation2];
        double o3 = - (o1*l1 + o2*l2) / l3;


        /* Set the lengths */
        if (!shape_set_vertex_length(shape, 0, l1)) continue;
        if (!shape_set_vertex_length(shape, 1, l2)) continue;
        if (!shape_set_vertex_length(shape, 2, l3)) continue;

        /* Set the angles */
        if (!shape_set_vertex_angle(shape, 0, 1, a12)) continue;
        if (!shape_set_vertex_angle(shape, 0, 2, a13)) continue;
        if (!shape_set_vertex_angle(shape, 1, 2, a23)) continue;

        /* Set the orientations */
        if (!shape_set_vertex_orientation(shape, 0, o1)) continue;
        if (!shape_set_vertex_orientation(shape, 1, o2)) continue;
        if (!shape_set_vertex_orientation(shape, 2, o3)) continue;

        /* Check the third orientation angle */
        if (o3 < sampleRawOrientation -> sampleArg -> min || o3 > sampleRawOrientation -> sampleArg -> max) continue;

        /* Check the shape */
        if (!shape_check(shape)) continue;

        /* Sort the shape */
        shape_sort(shape);

        /* Insert and sort */
        sample_shape_insort(sampleShapeTri, sampleShapeTri -> size + 1, &success);

        /* Set the parity */
        if (success)
          {
            shape_set_parity(shape);

            sampleShapeTri -> sizeParity += (size_t) (shape -> parity);
            sampleShapeTri -> sizeNoParity += (size_t) (!shape -> parity);
            sampleShapeTri -> sizeFull += 1 + (size_t) (!shape -> parity);
          }
      }

    /* Reallocate memory */
    for (size_t n = sampleShapeTri -> size; n < size + 1; n++)
      {
        sampleShapeTri -> arrayShape[n] = shape_free(sampleShapeTri -> arrayShape[n]);
      }

    sampleShapeTri -> arrayShape = realloc(sampleShapeTri -> arrayShape, sizeof(shape_t*) * sampleShapeTri -> size);

    /* Set the bin-volume */
    avr_shape_sample_vol(sampleShapeTri);

    /* Output the sample */
    if (fileName != NULL)
      {
        sample_shape_output(fileName, sampleShapeTri, NULL);
      }

    return sampleShapeTri;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------------   Output and Input   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int sample_shape_output(const char *fileName, sample_shape_t *sampleShape, int *prec)
{
    /*

        Output a sample_shape_t struct to file

    */

    /* Shape dim */
    size_t dim = sampleShape -> dim;

    /* Transform sample_shape_t struct to dat_t struct */
    dat_t *dat = dat_new(2*dim + dim * (dim - 1) / 2, 2, sampleShape -> size);

    for (size_t n = 0; n < sampleShape -> size; n++)
      {
        /* Shape */
        shape_t *shape = sampleShape -> arrayShape[n];

        /* Lengths, Orientations and Angles */
        for (size_t i = 0; i < dim; i++)
          {
            dat_set_xvalue(dat, i, n, shape_get_vertex_length(shape, i));
            dat_set_xvalue(dat, i + dim, n, shape_get_vertex_orientation(shape, i));

            for (size_t j = i + 1; j < dim; j++)
                dat_set_xvalue(dat, i * dim - (i + 1) * (i + 2) / 2 + j + 2*dim, n, shape_get_vertex_angle(shape, i, j));
          }

        /* Bin Volume */
        dat_set_yvalue(dat, 0, n, shape_get_volume(shape));
        dat_set_yvalue(dat, 1, n, shape_get_volume_err(shape));
      }

    /* Output the data */
    char *outFile = misc_scat(3, __FISHERLSSDIR__, _outDir, fileName);
    dat_output(outFile, dat, prec);
    free(outFile);

    /* Free memory */
    dat = dat_free(dat);

    return 0;
}


sample_shape_t *sample_shape_input(const char *fileName, sample_arg_t *sampleArgLength, sample_arg_t *sampleArgOrientation, size_t *size)
{
    /*

        Input a sample_shape_t struct from file. Need sampleArgLength and sampleArgOrientation to obtain the correct step size. (TODO: Maybe output step also to the file...)

        If size is not NULL, draw a random sample from the shapes

    */

    /* Sample length */
    sample_arg_t *sampleArgLength_ = sample_arg_new();
    sample_raw_t *sampleRawLength_ = sample_raw_new();
    sampleRawLength_ -> sampleArg = sampleArgLength_;

    /* Sample orientation */
    sample_arg_t *sampleArgOrientation_ = sample_arg_new();
    sample_raw_t *sampleRawOrientation_ = sample_raw_new();
    sampleRawOrientation_ -> sampleArg = sampleArgOrientation_;

    /* Read the data */
    char *outFile = misc_scat(3, __FISHERLSSDIR__, _outDir, fileName);
    dat_t *dat = dat_input(outFile, NULL, true);
    free(outFile);

    /* Shape dim */
    size_t dim = (size_t) (-3. + sqrt(8. * (double) dat -> xDim + 9)) / 2;

    /* Sample size and indices */
    size_t sampleSize = (size == NULL) ? dat -> size : *size;
    size_t *sampleIndices = misc_random_sample(0, dat -> size - 1, &sampleSize, false);

    /* New sample_shape_t struct */
    sample_shape_t *sampleShape = sample_shape_new(dim, sampleSize);

    for (size_t n = 0; n < sampleSize; n++)
      {
        /* Indices */
        size_t index = sampleIndices[n];

        /* Shape */
        shape_t *shape = sampleShape -> arrayShape[n];

        /* Lengths, Orientations and Angles */
        for (size_t i = 0; i < dim; i++)
          {
            shape_set_vertex_length(shape, i, dat_get_xvalue(dat, i, index));
            shape_set_vertex_orientation(shape, i, dat_get_xvalue(dat, i + dim, index));

            for (size_t j = i + 1; j < dim; j++)
                shape_set_vertex_angle(shape, i, j, dat_get_xvalue(dat, i * dim - (i + 1) * (i + 2) / 2 + j + 2*dim, index));

            misc_insort(&sampleRawLength_ -> array, &sampleArgLength_ -> size, dat_get_xvalue(dat, i, index), __ABSTOL__, NULL);
            misc_insort(&sampleRawOrientation_ -> array, &sampleArgOrientation_ -> size, dat_get_xvalue(dat, i + dim, index), __ABSTOL__, NULL);
          }

        /* Bin Volume */
        shape_set_volume(shape, dat_get_yvalue(dat, 0, index));
        shape_set_volume_err(shape, dat_get_yvalue(dat, 1, index));

        /* Parity */
        shape_set_parity(shape);

        sampleShape -> sizeParity += (size_t) (shape -> parity);
        sampleShape -> sizeNoParity += (size_t) (!shape -> parity);
        sampleShape -> sizeFull += 1 + (size_t) (!shape -> parity);
      }

    /* Sort the shapes (in case shapes in file are not sorted yet) */
    sample_shape_sort(sampleShape);

    /* Set the raw samples */
    sampleArgLength_ -> min = sampleRawLength_ -> array[0];
    sampleArgLength_ -> max = sampleRawLength_ -> array[sampleArgLength_ -> size - 1];
    sampleArgLength_ -> step = sampleArgLength -> step;
    sampleArgLength_ -> size = sampleArgLength -> size; // TODO

    sampleArgOrientation_ -> min = sampleRawOrientation_ -> array[0];
    sampleArgOrientation_ -> max = sampleRawOrientation_ -> array[sampleArgOrientation_ -> size - 1];
    sampleArgOrientation_ -> step = sampleArgOrientation -> step;
    sampleArgOrientation_ -> size = sampleArgOrientation -> size; // TODO

    sampleShape -> sampleRawLength = sampleRawLength_;
    sampleShape -> sampleRawOrientation = sampleRawOrientation_;

    /* Free memory */
    dat = dat_free(dat);
    free(sampleIndices);

    return sampleShape;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


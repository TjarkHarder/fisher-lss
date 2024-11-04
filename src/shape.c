/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SHAPE.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "shape.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SHAPE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   New / Free / Cp Functions   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


shape_t *shape_new(size_t dim)
{
    /*

        Create a new shape_t struct

    */

    shape_t *shape = malloc(sizeof(shape_t));

    shape -> dim = dim;

    shape -> length = calloc(dim, sizeof(double));
    shape -> orientation = calloc(dim, sizeof(double));
    shape -> angle = calloc(dim * (dim - 1) / 2, sizeof(double));

    shape -> volume = 0.;
    shape -> volumeErr = 0.;

    shape -> parity = false;

    return shape;
}


shape_t *shape_free(shape_t *shape)
{
    /*

        Free shape

    */

    /* If shape is NULL return NULL */
    if (shape == NULL)
        return NULL;

    /* Free shape's contents */
    free(shape -> length);
    free(shape -> orientation);
    free(shape -> angle);

    /* Free shape itself */
    free(shape);

    return NULL;

}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


shape_t *shape_cp(shape_t *shape)
{
    /*

        Copy shape

    */

    shape_t *shapeCp = shape_new(shape -> dim);

    for (size_t i = 0; i < shapeCp -> dim; i++)
      {
        shape_set_vertex_length(shapeCp, i, shape_get_vertex_length(shape, i));
        shape_set_vertex_orientation(shapeCp, i, shape_get_vertex_orientation(shape, i));

        for (size_t j = i + 1; j < shapeCp -> dim; j++)
            shape_set_vertex_angle(shapeCp, i, j, shape_get_vertex_angle(shape, i, j));
      }

    shapeCp -> volume = shape -> volume;
    shapeCp -> volumeErr = shape -> volumeErr;

    shapeCp -> parity = shape -> parity;

    return shapeCp;
}


shape_t *shape_cp_parity(shape_t *shape)
{
    /*

        Copy shape under parity

    */

    shape_t *shapeCp = shape_new(shape -> dim);

    for (size_t i = 0; i < shapeCp -> dim; i++)
      {
        shape_set_vertex_length(shapeCp, i, shape_get_vertex_length(shape, i));
        shape_set_vertex_orientation(shapeCp, i, -shape_get_vertex_orientation(shape, i));

        for (size_t j = i + 1; j < shapeCp -> dim; j++)
            shape_set_vertex_angle(shapeCp, i, j, shape_get_vertex_angle(shape, i, j));
      }

    shapeCp -> volume = shape -> volume;
    shapeCp -> volumeErr = shape -> volumeErr;

    shapeCp -> parity = shape -> parity;

    return shapeCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Swap Shapes and Vertices   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int shape_swap(shape_t **shape1, shape_t **shape2)
{
    /*

        Swap the addresses of *shape1 and *shape2

    */

    shape_t *shape = *shape1;

    *shape1 = *shape2;
    *shape2 = shape;

    return 0;
}


int shape_swap_vertex(shape_t *shape, size_t index1, size_t index2)
{
    /*

        Swap the index1'th and index2'th vertices

    */

    double buf;

    /* Swap vertex lengths */
    buf = shape_get_vertex_length(shape, index1);
    shape_set_vertex_length(shape, index1, shape_get_vertex_length(shape, index2));
    shape_set_vertex_length(shape, index2, buf);

    /* Swap vertex orientations */
    buf = shape_get_vertex_orientation(shape, index1);
    shape_set_vertex_orientation(shape, index1, shape_get_vertex_orientation(shape, index2));
    shape_set_vertex_orientation(shape, index2, buf);

    /* Swap vertex angles */
    for (size_t i = 0; i < shape -> dim; i++)
      {
        if (i == index1 || i == index2)
            continue;

        buf = shape_get_vertex_angle(shape, index1, i);
        shape_set_vertex_angle(shape, index1, i, shape_get_vertex_angle(shape, index2, i));
        shape_set_vertex_angle(shape, index2, i, buf);
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Setters   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool shape_set_vertex_length(shape_t *shape, size_t index, double length)
{
    /*

        Set the length of the index'th vertex

    */

    if (index >= shape -> dim)
      {
        printf("Cannot set the %ld'th vertex length if the shape only has %ld vertexes.\n", index, shape -> dim);
        exit(1);

        return false;
      }

    if (!isnormal(length) || length < 0.)
      {
        return false;
      }

    shape -> length[index] = length;

    return true;
}


bool shape_set_vertex_orientation(shape_t *shape, size_t index, double orientation)
{
    /*

        Set the orientation of the index'th vertex

    */

    if (index >= shape -> dim)
      {
        printf("Cannot set the %ld'th vertex orientation if the shape only has %ld vertexes.\n", index, shape -> dim);
        exit(1);

        return false;
      }

    if ((!isnormal(orientation) && orientation != 0.) || fabs(orientation) - 1. > 0.)
      {
        return false;
      }

    shape -> orientation[index] = orientation;

    return true;
}


bool shape_set_vertex_angle(shape_t *shape, size_t index1, size_t index2, double angle)
{
    /*

        Set the angle between the index1'th and index2'th vertexes

    */

    if (index1 >= shape -> dim || index2 >= shape -> dim)
      {
        printf("Cannot set the angle between the %ld'th and %ld'th vertexes if the shape only has %ld vertexes.\n", index1, index2, shape -> dim);
        exit(1);

        return false;
      }

    if ((!isnormal(angle) && angle != 0.) || fabs(angle) - 1. > 0.)
      {
        return false;
      }

    size_t index = shape_get_vertex_angle_index(shape -> dim, index1, index2);

    shape -> angle[index] = angle;

    return true;
}


bool shape_set_volume(shape_t *shape, double volume)
{
    /*

        Set the bin-volume of the shape

    */

    if (!isnormal(volume))
      {
        return false;
      }

    shape -> volume = volume;

    return true;
}


bool shape_set_volume_err(shape_t *shape, double volumeErr)
{
    /*

        Set the error of the bin-volume of the shape

    */

    if (!isnormal(volumeErr))
      {
        return false;
      }

    shape -> volumeErr = volumeErr;

    return true;
}


bool shape_set_parity(shape_t *shape)
{
    /*

        Set whether the shape is parity invariant or not

    */

    /* Compare the shape itself and set the parity flag */
    int comp = shape_comp(shape, shape);
    shape -> parity = (comp == 2); // Invariant under parity only if the shape is equal to itself and its parity transformed version

    return shape -> parity;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------------   Getters   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double shape_get_vertex_length(shape_t *shape, size_t index)
{
    /*

        Get the length of the index'th vertex

    */

    if (index >= shape -> dim)
      {
        printf("Cannot get the %ld'th vertex length if the shape only has %ld vertexes.\n", index, shape -> dim);
        exit(1);

        return NAN;
      }

    return shape -> length[index];
}


double shape_qget_vertex_length(shape_t *shape, size_t index)
{
    /*

        Get the length of the index'th vertex

    */

    return shape -> length[index];
}


/*  ------------------------------------------------------------------------------------------------------  */


double shape_get_vertex_orientation(shape_t *shape, size_t index)
{
    /*

        Get the orientation of the index'th vertex

    */

    if (index >= shape -> dim)
      {
        printf("Cannot get the %ld'th vertex orientation if the shape only has %ld vertexes.\n", index, shape -> dim);
        exit(1);

        return NAN;
      }

    return shape -> orientation[index];
}


double shape_qget_vertex_orientation(shape_t *shape, size_t index)
{
    /*

        Get the orientation of the index'th vertex

    */

    return shape -> orientation[index];
}


/*  ------------------------------------------------------------------------------------------------------  */


double shape_get_vertex_angle(shape_t *shape, size_t index1, size_t index2)
{
    /*

        Get the angle between the index1'th and index2'th vertexes

    */

    if (index1 >= shape -> dim || index2 >= shape -> dim)
      {
        printf("Cannot get the angle between the %ld'th and %ld'th vertexes if the shape only has %ld vertexes.\n", index1, index2, shape -> dim);
        exit(1);

        return NAN;
      }

    size_t index = shape_get_vertex_angle_index(shape -> dim, index1, index2);

    return shape -> angle[index];
}


double shape_qget_vertex_angle(shape_t *shape, size_t index1, size_t index2)
{
    /*

        Get the angle between the index1'th and index2'th vertexes

    */

    size_t index = shape_get_vertex_angle_index(shape -> dim, index1, index2);

    return shape -> angle[index];
}


/*  ------------------------------------------------------------------------------------------------------  */


double shape_get_volume(shape_t *shape)
{
    /*

        Get the bin-volume

    */

    return shape -> volume;
}


double shape_get_volume_err(shape_t *shape)
{
    /*

        Get the bin-volume error

    */

    return shape -> volumeErr;
}


bool shape_get_parity(shape_t *shape)
{
    /*

        Get the parity

    */

    return shape -> parity;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t shape_get_vertex_angle_index(size_t dim, size_t index1, size_t index2)
{
    /*

        Get the index of the angle between the index1'th and index2'th vertexes

    */

    /* Cannot have equal indices */
    if (index1 == index2)
      {
        printf("Cannot get the index of the angle between the same vertexes (angle is always one, thus ignored).\n");
        exit(1);

        return 0;
      }

    /* Too few vertexes */
    if (index1 >= dim || index2 >= dim)
      {
        printf("Cannot get index of the angle between the %ld'th and %ld'th vertexes if there are only %ld vertexes..\n", index1, index2, dim);
        exit(1);

        return 0;
      }

    /* Let i be the smaller index and j the larger */
    size_t i = index1;
    size_t j = index2;

    if (i > j)
      {
        i = index2;
        j = index1;
      }

    return i * (dim - 1) - i * (i + 1) / 2 + j - 1;
}


size_t shape_qget_vertex_angle_index(size_t dim, size_t index1, size_t index2)
{
    /*

        Get the index of the angle between the index1'th and index2'th vertexes

    */

    /* Let i be the smaller index and j the larger */
    size_t i = index1;
    size_t j = index2;

    if (i > j)
      {
        i = index2;
        j = index1;
      }

    return i * (dim - 1) - i * (i + 1) / 2 + j - 1;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Check the Shape   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool shape_check(shape_t *shape)
{
    /*

        Check if the shape's variables obey x1_ + ... = 0_!

    */

    /* Calculate and check the polar angles */
    double *cphi = malloc(sizeof(double) * shape -> dim);

    cphi[0] = 1.; // The first vector can be chosen to lie in x-z plane

    for (size_t i = 1; i < shape -> dim; i++)
      {
        size_t angleIndex = shape_get_vertex_angle_index(shape -> dim, 0, i);
        cphi[i] = (fabs(shape -> angle[angleIndex]) - 1. == 0.)
                    ? shape -> angle[angleIndex]
                    : (shape -> angle[angleIndex] - shape -> orientation[0] * shape -> orientation[i])
                        / sqrt(1. - shape -> orientation[0] * shape -> orientation[0])
                        / sqrt(1. - shape -> orientation[i] * shape -> orientation[i]);

        if ((!isnormal(cphi[i]) && cphi[i] != 0.) || fabs(cphi[i]) - 1. > 0.)
          {
            free(cphi);
            return false;
          }
      }

    /* Verify x1_ + ... = 0 */
    double sum[3] = {0., 0., 0.};

    for (size_t i = 0; i < shape -> dim; i++)
      {
        sum[0] += shape -> length[i] * sqrt(1. - shape -> orientation[i]*shape -> orientation[i]) * cphi[i];
        sum[2] += shape -> length[i] * shape -> orientation[i];
      }

    for (size_t i = 0; i < 3; i++)
      {
        if (fabs(sum[i]) > __ABSTOL__ * (double) shape -> dim)
          {
            free(cphi);
            return false;
          }
      }

    /* Free memory */
    free(cphi);

    return true;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Parity Transform Shape   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int shape_parity(shape_t *shape)
{
    /*

        Parity transform shape

    */

    for (size_t i = 0; i < shape -> dim; i++)
        shape -> orientation[i] *= -1.;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Sort Shape's Vertices   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _sort_comp_length(shape_t *shape, int index1, int index2);
static bool _sort_comp_orientation(shape_t *shape, int index1, int index2);

static int _sort_kern(shape_t *shape, int low, int high, bool (*comp)(shape_t*, int, int));
static int _sort_partition(shape_t *shape, int low, int high, bool (*comp)(shape_t*, int, int));

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int shape_sort(shape_t *shape)
{
    /*

        Sort the vertices in shape from shortest to longest, and the shortest vertices from positive smallest
        to largest orientation

    */

    /**  Sort the vertices' lengths  **/

    _sort_kern(shape, 0, (int) shape -> dim - 1, _sort_comp_length);


    /**  Sort the orientations  **/

    /* Find the smallest value by magnitude */
    double min = shape_qget_vertex_orientation(shape, 0);

    for (size_t i = 1; i < shape -> dim; i++)
      {
        /* Only consider the orientations of the smallest length */
        if (shape_qget_vertex_length(shape, 0) != shape_qget_vertex_length(shape, i))
            break;

        min = (fabs(min) > fabs(shape_qget_vertex_orientation(shape, i))) ? shape_qget_vertex_orientation(shape, i) : min;
      }

    /* If smallest orientation by magnitude is negative, invert the shape */
    if (min < 0)
        shape_parity(shape);

    /* Sort the orientations for the same lengths */
    int low = 0;
    int high = 1;

    for (size_t i = 1; i < shape -> dim; i++)
      {
        /* Length are not the same anymore -> sort the orientations */
        if (shape_qget_vertex_length(shape, (size_t) low) != shape_qget_vertex_length(shape, (size_t) high))
          {
            _sort_kern(shape, low, high - 1, _sort_comp_orientation);

            low = high;
          }

        high++;
      }

    /* Sort one orientations one last time */
    _sort_kern(shape, low, high - 1, _sort_comp_orientation);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static bool _sort_comp_length(shape_t *shape, int index1, int index2)
{
    /*

        Compare the index1'th and index2'th vertices' lengths

    */

    /* Place lowest length to the left */
    return ((shape_qget_vertex_length(shape, (size_t) index1) - shape_qget_vertex_length(shape, (size_t) index2)) <= 0.);
}


static bool _sort_comp_orientation(shape_t *shape, int index1, int index2)
{
    /*

        Compare the index1'th and index2'th vertices' orientations

    */

    /* Orientations */
    double orientation1 = shape_qget_vertex_orientation(shape, (size_t) index1);
    double orientation2 = shape_qget_vertex_orientation(shape, (size_t) index2);

    /* 0 always to the left */
    if (orientation1 == 0.)
        return true;

    /* 0 always to the left */
    if (orientation2 == 0.)
        return false;

    /* Same sign : Place vertices with smaller orientation to the left; Opposite sign : Place vertices with larger orientation to the left (i.e., the positive one) */
    return (orientation1 * orientation2 > 0.) ? orientation1 <= orientation2 : orientation1 > orientation2;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _sort_kern(shape_t *shape, int low, int high, bool (*comp)(shape_t*, int, int))
{
    /*

        Kernel of the quicksort function

    */

    /* Sort as long as pivot is larger than low */
    if (low < high)
      {
        /* Place last pivot, and get new pivot */
        int pivot = _sort_partition(shape, low, high, comp);

        /* Sort to the left of the previous pivot */
        _sort_kern(shape, low, pivot - 1, comp);

        /* Sort to the right of the last pivot */
        _sort_kern(shape, pivot + 1, high, comp);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _sort_partition(shape_t *shape, int low, int high, bool (*comp)(shape_t*, int, int))
{
    /*

        Partition the vertices

    */

    /* Push i'th vertex to the right */
    int i = low - 1;

    for (int j = low; j < high; j++)
      {
        /* Push j'th vertex to the left if comp is true */
        if (comp(shape, j, high))
          {
            shape_swap_vertex(shape, (size_t) ++i, (size_t) j);
          }
      }

    /* Swap i'th vertex with pivot's vertex */
    shape_swap_vertex(shape, (size_t) i + 1, (size_t) high);

    return i + 1;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Compare Shapes   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool shape_ecomp(shape_t *shape1, shape_t *shape2)
{
    /*

        Compare two shapes and true if they are the same and false if not. Assumes that shapes are sorted
        according to shape_sort()

    */

    /* Shapes must have same dimension */
    if (shape1 -> dim != shape2 -> dim)
        return false;

    /* Check if the shapes are the same */
    for (size_t i = 0; i < shape1 -> dim; i++)
      {
        /* Same vertex length and orientation */
        if (fabs(shape1 -> length[i] - shape2 -> length[i]) > __ABSTOL__ || fabs(shape1 -> orientation[i] - shape2 -> orientation[i]) > __ABSTOL__)
            return false;
      }

    return true;
}


int shape_comp(shape_t *shape1, shape_t *shape2)
{
    /*

        Compare two shapes and return 0 if they are not the same, 1 if they are the same but not related by parity,
        -1 if they are related by parity but not the same, and 2 if they are the same and related by parity.

    */

    /* Shapes must have same dimension */
    if (shape1 -> dim != shape2 -> dim)
        return 0;

    /* Copy one shape to permutate the vertexes */
    shape_t *shapeCp = shape_cp(shape2);

    /* Results */
    bool same = false;
    bool parity = false;

    /* Check if the shapes are the same */
    for (size_t i1 = 0; i1 < shape1 -> dim; i1++)
      {
        for (size_t i2 = i1; i2 < shape2 -> dim; i2++)
          {
            /* Same vertex length and orientation */
            if (fabs(shape1 -> length[i1] - shapeCp -> length[i2]) < __ABSTOL__ && fabs(shape1 -> orientation[i1] - shapeCp -> orientation[i2]) < __ABSTOL__)
              {
                /* Rearange shapeCp's arrays as in shape1 */
                misc_swap(&(shapeCp -> length[i1]), &(shapeCp -> length[i2]), "f");
                misc_swap(&(shapeCp -> orientation[i1]), &(shapeCp -> orientation[i2]), "f");

                break;
              }

            /* i1'th vertex of shape1 is not in shape2 -> go to checkParity (cannot return 0 without checking parity configuration) */
            if (i2 == shape2 -> dim - 1)
              {
                goto checkParity;
              }
          }

        /* Can only get to the end if the shapes match */
        if (i1 == shape1 -> dim - 1)
          {
            same = true;
          }
      }

    /* Check if the shapes are related by parity */
    checkParity:

    for (size_t i1 = 0; i1 < shape1 -> dim; i1++)
      {
        for (size_t i2 = i1; i2 < shape2 -> dim; i2++)
          {
            /* Same vertex length and absolute orientation */
            if (fabs(shape1 -> length[i1] - shapeCp -> length[i2]) < __ABSTOL__ && fabs(shape1 -> orientation[i1] + shapeCp -> orientation[i2]) < __ABSTOL__)
              {
                /* Rearange shapeCp's arrays as in shape1 */
                misc_swap(&(shapeCp -> length[i1]), &(shapeCp -> length[i2]), "f");
                misc_swap(&(shapeCp -> orientation[i1]), &(shapeCp -> orientation[i2]), "f");

                break;
              }

            /* i1'th vertex of shape1 is not in parity transformed shape2 */
            if (i2 == shape2 -> dim - 1)
              {
                goto checkFinished;
              }
          }

        /* Can only get to the end if the shapes match */
        if (i1 == shape1 -> dim - 1)
          {
            parity = true;
          }
      }

    /* Finished checking */
    checkFinished:

    /* Free shapeCp */
    shapeCp = shape_free(shapeCp);

    /* Shapes are the same and related by parity */
    if (same && parity)
        return 2;

    /* Shapes are the same but not related by parity */
    if (same)
        return 1;

    /* Shapes are related by parity but not the same */
    if (parity)
        return -1;

    /* Shapes are neither the same nor related by parity */
    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool shape_comp_vertex_bin(shape_t *shape1, size_t index1, shape_t *shape2, size_t index2, double lengthTol, double orientationTol)
{
    /*

        Compare the index1'th vertex's bin of shape1 with the index2'th vertex's bin of shape2

    */

    /* Index out of bounds */
    if (index1 >= shape1 -> dim || index2 >= shape2 -> dim)
      {
        printf("Cannot compare the %ld'th vertex of the first shape with the %ld'th vertex of the second shape if their number of vertexes are %ld and %ld respectively\n", index1, index2, shape1 -> dim, shape2 -> dim);
        exit(1);

        return false;
      }

    /* Must have matching length and orientation */
    if (fabs(shape1 -> length[index1] - shape2 -> length[index2]) >= lengthTol || fabs(shape1 -> orientation[index1] - shape2 -> orientation[index2]) >= orientationTol)
        return false;

    return true;
}


bool **shape_comp_vertexes_bin(shape_t *shape1, shape_t *shape2, double lengthTol, double orientationTol)
{
    /*

        Compare shape1's and shape2's vertexes' bins to see if any match

    */

    bool **comp = malloc(sizeof(bool*) * shape1 -> dim);

    for (size_t i = 0; i < shape1 -> dim; i++)
      {
        comp[i] = malloc(sizeof(bool) * shape2 -> dim);

        for (size_t j = 0; j < shape2 -> dim; j++)
          {
//            comp[i][j] = shape_comp_vertex(shape1, i, shape2, j);
            comp[i][j] = (fabs(shape1 -> length[i] - shape2 -> length[j]) > lengthTol || fabs(shape1 -> orientation[i] - shape2 -> orientation[j]) > orientationTol) ? false : true;
          }
      }

    return comp;
}


/*  ------------------------------------------------------------------------------------------------------  */


bool shape_comp_vertex(shape_t *shape1, size_t index1, shape_t *shape2, size_t index2)
{
    /*

        Compare the index1'th vertex of shape1 with the index2'th vertex of shape2

    */

    return shape_comp_vertex_bin(shape1, index1, shape2, index2, __ABSTOL__, __ABSTOL__);
}


bool **shape_comp_vertexes(shape_t *shape1, shape_t *shape2)
{
    /*

        Compare shape1's and shape2's vertexes to see if any match

    */

    return shape_comp_vertexes_bin(shape1, shape2, __ABSTOL__, __ABSTOL__);
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool shape_comp_vertex_parity_bin(shape_t *shape1, size_t index1, shape_t *shape2, size_t index2, double lengthTol, double orientationTol)
{
    /*

        Compare the index1'th vertex's bin of shape1 with the index2'th parity transformed vertex's bin of shape2

    */

    /* Index out of bounds */
    if (index1 >= shape1 -> dim || index2 >= shape2 -> dim)
      {
        printf("Cannot compare the %ld'th vertex of the first shape with the %ld'th vertex of the second shape if their number of vertexes are %ld and %ld respectively\n", index1, index2, shape1 -> dim, shape2 -> dim);
        exit(1);

        return false;
      }

    /* Must have matching length and orientation */
    if (fabs(shape1 -> length[index1] - shape2 -> length[index2]) >= lengthTol || fabs(shape1 -> orientation[index1] + shape2 -> orientation[index2]) >= orientationTol)
        return false;

    return true;
}


bool **shape_comp_vertexes_parity_bin(shape_t *shape1, shape_t *shape2, double lengthTol, double orientationTol)
{
    /*

        Compare shape1's vertexes' bins to and shape2's parity transformed vertexes' bins to see if any match

    */

    bool **comp = malloc(sizeof(bool*) * shape1 -> dim);

    for (size_t i = 0; i < shape1 -> dim; i++)
      {
        comp[i] = malloc(sizeof(bool) * shape2 -> dim);

        for (size_t j = 0; j < shape2 -> dim; j++)
          {
//            comp[i][j] = shape_comp_vertex(shape1, i, shape2, j);
            comp[i][j] = (fabs(shape1 -> length[i] - shape2 -> length[j]) > lengthTol || fabs(shape1 -> orientation[i] + shape2 -> orientation[j]) > orientationTol) ? false : true;
          }
      }

    return comp;
}


/*  ------------------------------------------------------------------------------------------------------  */


bool shape_comp_vertex_parity(shape_t *shape1, size_t index1, shape_t *shape2, size_t index2)
{
    /*

        Compare the index1'th vertex of shape1 with the index2'th parity transformed vertex of shape2

    */

    return shape_comp_vertex_parity_bin(shape1, index1, shape2, index2, __ABSTOL__, __ABSTOL__);
}


bool **shape_comp_vertexes_parity(shape_t *shape1, shape_t *shape2)
{
    /*

        Compare shape1's vertexes to shape2's parity transformed vertexes to see if any match

    */

    return shape_comp_vertexes_parity_bin(shape1, shape2, __ABSTOL__, __ABSTOL__);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Degeneracy of Edges   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t *shape_degen_vertexes(shape_t *shape)
{
    /*

        Get the degeneracies of the vertexes (in length and orientation).

        Given two degenerate vertexes at indices i and j in the length and orientation arrays, the returned array
        will contain a two at the i'th position and a zero at the j'th position.

    */

    size_t *degen = malloc(sizeof(size_t) * shape -> dim);

    /* Initialise degen to ones */
    for (size_t i = 0; i < shape -> dim; i++)
      {
        degen[i] = 1;
      }

    /* Get the degeneracies */
    for (size_t i = 0; i < shape -> dim; i++)
      {
        /* Skip those that are zero */
        if (degen[i] == 0)
            continue;

        /* Compare the i'th vertex with the other vertexes */
        for (size_t j = i + 1; j < shape -> dim; j++)
          {
            if (shape_comp_vertex(shape, i, shape, j))
              {
                degen[i] += 1;
                degen[j] -= 1;
              }
          }
      }

    return degen;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

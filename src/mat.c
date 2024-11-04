/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     MAT.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "mat.h"


// TODO: Maybe add an additional header to the output function...


/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     MATRIX TYPES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Matrix File Parameters  **/

static mat_file_t *_mat_file_new()
{
    /*

        Create a new mat_file_t struct

    */

    mat_file_t *mfile = malloc(sizeof(mat_file_t));

    mfile -> fileName = NULL;
    mfile -> type = NULL;

    return mfile;
}


static mat_file_t *_mat_file_free(mat_file_t *mfile)
{
    /*

        Free mfile

    */

    if (mfile == NULL)
        return NULL;

    free(mfile -> fileName);
    free(mfile -> type);

    free(mfile);

    return NULL;
}


static mat_file_t *_mat_file_cp(mat_file_t *mfile)
{
    /*

        Copy mfile

    */

    if (mfile == NULL)
        return NULL;

    mat_file_t *mfileCp = _mat_file_new();

    mfileCp -> fileName = misc_scat(1, mfile -> fileName);

    mfileCp -> type = misc_scat(1, mfile -> type);
    mfileCp -> loc = mfile -> loc;
    mfileCp -> binary = mfile -> binary;

    return mfileCp;
}



/**  Matrix Types  **/

#ifndef __MAT_NUMBER_OF_TYPES__
#define __MAT_NUMBER_OF_TYPES__ 6
#endif

/* Number of types */
const size_t _matTypesNum = __MAT_NUMBER_OF_TYPES__;

/* Type ID's (in order of least memory) */
const char *_matTypesID[__MAT_NUMBER_OF_TYPES__] = {"null", "d", "s", "ut", "lt", "f"};

/* (Transposed) Type ID's (in order of least memory) */
const char *_matTpTypesID[__MAT_NUMBER_OF_TYPES__] = {"null", "d", "s", "lt", "ut", "f"}; // ie every type is 'transposed'


/* Null Matrix */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_type_index_null(size_t *dim, size_t *loc, size_t *index);
static int _mat_type_loc_null(size_t *dim, size_t index, size_t *loc);
static bool _mat_type_comp_null(void *matrix, size_t *loc, double tol);
static int _mat_type_rbounds_null(size_t *dim, size_t col, size_t *rbounds);
static int _mat_type_cbounds_null(size_t *dim, size_t row, size_t *cbounds);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_type_t _matTypeNull = {"null", 0, _mat_type_index_null, _mat_type_loc_null, _mat_type_comp_null, _mat_type_rbounds_null, _mat_type_cbounds_null};


/* Full Matrix */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_type_index_f(size_t *dim, size_t *loc, size_t *index);
static int _mat_type_loc_f(size_t *dim, size_t index, size_t *loc);
static bool _mat_type_comp_f(void *matrix, size_t *loc, double tol);
static int _mat_type_rbounds_f(size_t *dim, size_t col, size_t *rbounds);
static int _mat_type_cbounds_f(size_t *dim, size_t row, size_t *cbounds);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_type_t _matTypeF = {"f", 0, _mat_type_index_f, _mat_type_loc_f, _mat_type_comp_f, _mat_type_rbounds_f, _mat_type_cbounds_f};


/* Diagonal Matrix */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_type_index_d(size_t *dim, size_t *loc, size_t *index);
static int _mat_type_loc_d(size_t *dim, size_t index, size_t *loc);
static bool _mat_type_comp_d(void *matrix, size_t *loc, double tol);
static int _mat_type_rbounds_d(size_t *dim, size_t col, size_t *rbounds);
static int _mat_type_cbounds_d(size_t *dim, size_t row, size_t *cbounds);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_type_t _matTypeD = {"d", 0, _mat_type_index_d, _mat_type_loc_d, _mat_type_comp_d, _mat_type_rbounds_d, _mat_type_cbounds_d};


/* Symmetric Matrix */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_type_index_s(size_t *dim, size_t *loc, size_t *index);
static int _mat_type_loc_s(size_t *dim, size_t index, size_t *loc);
static bool _mat_type_comp_s(void *matrix, size_t *loc, double tol);
static int _mat_type_rbounds_s(size_t *dim, size_t col, size_t *rbounds);
static int _mat_type_cbounds_s(size_t *dim, size_t row, size_t *cbounds);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_type_t _matTypeS = {"s", 1, _mat_type_index_s, _mat_type_loc_s, _mat_type_comp_s, _mat_type_rbounds_s, _mat_type_cbounds_s};


/* Upper Triangular Matrix */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_type_index_ut(size_t *dim, size_t *loc, size_t *index);
static int _mat_type_loc_ut(size_t *dim, size_t index, size_t *loc);
static bool _mat_type_comp_ut(void *matrix, size_t *loc, double tol);
static int _mat_type_rbounds_ut(size_t *dim, size_t col, size_t *rbounds);
static int _mat_type_cbounds_ut(size_t *dim, size_t row, size_t *cbounds);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_type_t _matTypeUt = {"ut", 0, _mat_type_index_ut, _mat_type_loc_ut, _mat_type_comp_ut, _mat_type_rbounds_ut, _mat_type_cbounds_ut};


/* Lower Triangular Matrix */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_type_index_lt(size_t *dim, size_t *loc, size_t *index);
static int _mat_type_loc_lt(size_t *dim, size_t index, size_t *loc);
static bool _mat_type_comp_lt(void *matrix, size_t *loc, double tol);
static int _mat_type_rbounds_lt(size_t *dim, size_t col, size_t *rbounds);
static int _mat_type_cbounds_lt(size_t *dim, size_t row, size_t *cbounds);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_type_t _matTypeLt = {"lt", 0, _mat_type_index_lt, _mat_type_loc_lt, _mat_type_comp_lt, _mat_type_rbounds_lt, _mat_type_cbounds_lt};


/* All types in order of least memory (NOTE: ORDER MUST BE THE SAME AS IN _matTypesID) */
const mat_type_t _matTypes[__MAT_NUMBER_OF_TYPES__] = {_matTypeNull, _matTypeD, _matTypeS, _matTypeUt, _matTypeLt, _matTypeF};



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------------   Getters   --------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_type_t *mat_type_get_struct(const char *id)
{
    /*

        Find the type for a matrix with given id

    */

    for (size_t i = 0; i < _matTypesNum; i++)
      {
        if (!strcmp(_matTypesID[i], id))
          {
            return (mat_type_t*) &_matTypes[i];
          }
      }

    printf("Could not find the matrix type for a matrix with ID '%s'.\n", id);
    exit(1);

    return NULL;
}


mat_type_t *mat_type_get_struct_tp(const char *id)
{
    /*

        Find the transposed type for a matrix with given id

    */

    for (size_t i = 0; i < _matTypesNum; i++)
      {
        if (!strcmp(_matTpTypesID[i], id))
          {
            return (mat_type_t*) &_matTypes[i];
          }
      }

    printf("Could not find the transposed matrix type for a matrix with ID '%s'.\n", id);
    exit(1);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Compare Matrices   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _mat_comp_matrix_o(mat_t *mat1, mat_t *mat2, double tol);
static bool _mat_comp_matrix_b(mat_t *mat1, mat_t *mat2, double tol);
static bool _mat_comp_matrix_g(mat_t *mat1, mat_t *mat2, double tol);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


bool mat_comp_matrix(mat_t *mat1, mat_t *mat2, double tol)
{
    /*

        Compare two matrices

    */

    /* Both matrices are NULL */
    if (mat1 == NULL && mat2 == NULL)
        return true;

    /* Only one matrix is NULL */
    if (mat1 == NULL)
        return mat2 -> mtype == &_matTypeNull;

    if (mat2 == NULL)
        return mat1 -> mtype == &_matTypeNull;

    /* Full dimensions */
    size_t fdim1[2];
    mat_get_fdim(mat1, fdim1);

    size_t fdim2[2];
    mat_get_fdim(mat2, fdim2);

    if (fdim1[0] != fdim2[0] || fdim1[1] != fdim2[1])
        return false;

    /* Compare matrix block structures */
    bool comp = mat_comp_bstruct(mat1, mat2);

    /* Ordinary matrices */
    if (comp && mat1 -> mblock == NULL)
      {
        return _mat_comp_matrix_o(mat1, mat2, tol);
      }

    /* Block matrices */
    if (comp && mat1 -> mblock != NULL)
      {
        return _mat_comp_matrix_b(mat1, mat2, tol);
      }

    /* General matrices */
    return _mat_comp_matrix_g(mat1, mat2, tol);
}


/*  ------------------------------------------------------------------------------------------------------  */


static bool _mat_comp_matrix_o(mat_t *mat1, mat_t *mat2, double tol)
{
    /*

        Compare two ordinary matrices

    */

    size_t comp = true;

    #pragma omp parallel for collapse(2) shared(comp)
    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat1 -> dim[1]; j++)
          {
            /* Skip if matrices are the same */
            if (!comp)
                continue;

            size_t loc[2] = {i, j};

            double val1 = mat_get_value(mat1, loc);
            double val2 = mat_get_value(mat2, loc);

            if (fabs(val1 - val2) > tol)
              {
                comp = false;
              }
          }
      }

    return comp;
}


static bool _mat_comp_matrix_b(mat_t *mat1, mat_t *mat2, double tol)
{
    /*

        Compare two block matrices

    */

    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat1 -> dim[1]; j++)
          {
            size_t loc[2] = {i, j};

            mat_t *matBlock1 = mat_fget_block(mat1, loc);
            mat_t *matBlock2 = mat_fget_block(mat2, loc);

            if (!mat_comp_matrix(matBlock1, matBlock2, tol))
              {
                matBlock1 = mat_ffree(matBlock1);
                matBlock2 = mat_ffree(matBlock2);

                return false;
              }
          }
      }

    return true;
}


static bool _mat_comp_matrix_g(mat_t *mat1, mat_t *mat2, double tol)
{
    /*

        Compare two ordinary matrices

    */

    size_t fdim[2] = {0, 0};
    mat_get_fdim(mat1, fdim);

    size_t comp = true;

    #pragma omp parallel for collapse(2) shared(comp)
    for (size_t i = 0; i < fdim[0]; i++)
      {
        for (size_t j = 0; j < fdim[1]; j++)
          {
            /* Skip if matrices are the same */
            if (!comp)
                continue;

            size_t loc[2] = {i, j};

            double val1 = mat_get_value(mat1, loc);
            double val2 = mat_get_value(mat2, loc);

            if (fabs(val1 - val2) > tol)
              {
                comp = false;
              }
          }
      }

    return comp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool mat_comp_bstruct(mat_t *mat1, mat_t *mat2)
{
    /*

        Compare the matrix block structs of mat1 and mat2

    */

    /* Check Dimensions */
    if (mat1 -> dim[0] != mat2 -> dim[0] || mat1 -> dim[1] != mat2 -> dim[1])
        return false;

    /* Check Block Status */
    if (mat1 -> mblock == NULL && mat2 -> mblock == NULL)
        return true;

    if ((mat1 -> mblock == NULL && mat2 -> mblock != NULL) || (mat1 -> mblock != NULL && mat2 -> mblock == NULL))
        return false;

    /* Check Blocks */
    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat1 -> dim[1]; j++)
          {
            /* Blocks at location */
            size_t loc[2] = {i, j};
            mat_t *mat1Block = mat_fget_block(mat1, loc);
            mat_t *mat2Block = mat_fget_block(mat2, loc);

            /* Atleast one of the blocks is NULL -> inconclusive, must check other blocks */
            if (mat1Block == NULL || mat2Block == NULL)
              {
                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);

                continue;
              }

            /* Compare blocks */
            if (!mat_comp_bstruct(mat1Block, mat2Block))
              {
                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);

                return false;
              }

            /* Free memory */
            mat1Block = mat_ffree(mat1Block);
            mat2Block = mat_ffree(mat2Block);
          }
      }

    return true;
}


/*  ------------------------------------------------------------------------------------------------------  */


bool mat_comp_bstruct_cr(mat_t *mat1, mat_t *mat2)
{
    /*

        Compare the matrix block structs of the columns of mat1 and the rows of mat2

    */

    /* Check Dimensions */
    if (mat1 -> dim[1] != mat2 -> dim[0])
        return false;

    /* Check Block Status */
    if (mat1 -> mblock == NULL && mat2 -> mblock == NULL)
        return true;

    if ((mat1 -> mblock == NULL && mat2 -> mblock != NULL) || (mat1 -> mblock != NULL && mat2 -> mblock == NULL))
        return false;

    /* Check Blocks */
    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat2 -> dim[1]; j++)
          {
            for (size_t k = 0; k < mat1 -> dim[1]; k++)
              {
                /* Blocks at location */
                size_t loc1[2] = {i, k};
                mat_t *mat1Block = mat_fget_block(mat1, loc1);

                size_t loc2[2] = {k, j};
                mat_t *mat2Block = mat_fget_block(mat2, loc2);

                /* Atleast one of the blocks is NULL -> inconclusive, must check other blocks */
                if (mat1Block == NULL || mat2Block == NULL)
                  {
                    mat1Block = mat_ffree(mat1Block);
                    mat2Block = mat_ffree(mat2Block);

                    continue;
                  }

                /* Compare blocks */
                if (!mat_comp_bstruct_cr(mat1Block, mat2Block))
                  {
                    mat1Block = mat_ffree(mat1Block);
                    mat2Block = mat_ffree(mat2Block);

                    return false;
                  }

                /* Free memory */
                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);
              }
          }
      }

    return true;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Matrix Types   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_type_index_null(size_t *dim, size_t *loc, size_t *index)
{
    /*

        Get the index of the (i,j) position (loc) for a 1-d array representation of a null matrix of dimensions dim.

    */

    (void) index;

    /* Make sure location is within dimension bounds */
    if (loc[0] >= dim[0] || loc[1] >= dim[1])
      {
        printf("Location (%ld,%ld) is out of bounds if the dimensions are (%ld,%ld).\n", loc[0], loc[1], dim[0], dim[1]);
        exit(1);

        return 1;
      }

    return -1;
}


static int _mat_type_loc_null(size_t *dim, size_t index, size_t *loc)
{
    /*

        Get the (i, j) position of the index position in a 1-d array representation of a null matrix of dimensions dim.

    */

    (void) loc;

    /* Make sure location is within dimension bounds */
    if (index >= dim[0] * dim[1])
      {
        printf("Index position %ld is out of bounds if the size of the matrix array is %ld.\n", index, dim[0] * dim[1]);
        exit(1);

        return 1;
      }

    return -1;
}


static bool _mat_type_comp_null(void *matrix, size_t *loc, double tol)
{
    /*

        Comparison function to determine if elements at loc in matrix satisfy requirements to be empty.

        If loc is NULL, only check if matrix' type satisfies requirements to be empty.

    */

    mat_t *mat = (mat_t*) matrix;

    /* Matrix already has type 'null' */
    if (mat -> mtype == &_matTypeNull)
        return true;

    /* loc is NULL, ie only check for types */
    if (loc == NULL)
        return false;

    /* Element does not exist */
    size_t index = 0;

    if (mat_get_index(mat, loc, &index) == -1)
        return true;

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Element must be close to zero */
        if (fabs(mat_get_value(mat, loc)) > tol)
            return false;
      }

    /* Block Matrix */
    else
      {
        mat_t *matBlock = mat_fget_block(mat, loc);

        /* Every element must be close to zero */
        for (size_t i = 0; i < matBlock -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matBlock, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                size_t locBlock[2] = {i, j};

                if (!_mat_type_comp_null(matBlock, locBlock, tol))
                  {
                    /* Free memory */
                    matBlock = mat_ffree(matBlock);

                    return false;
                  }
              }
          }

        /* Free memory */
        matBlock = mat_ffree(matBlock);
      }

    return true;
}


static int _mat_type_rbounds_null(size_t *dim, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds for a column of a null matrix

    */

    (void) dim;
    (void) col;

    rbounds[0] = 0;
    rbounds[1] = 0;

    return -1;
}


static int _mat_type_cbounds_null(size_t *dim, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds for a row of a null matrix

    */

    (void) dim;
    (void) row;

    cbounds[0] = 0;
    cbounds[1] = 0;

    return -1;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_type_index_f(size_t *dim, size_t *loc, size_t *index)
{
    /*

        Get the index of the (i,j) position (loc) for a 1-d array representation of a full matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (loc[0] >= dim[0] || loc[1] >= dim[1])
      {
        printf("Location (%ld,%ld) is out of bounds if the dimensions are (%ld,%ld).\n", loc[0], loc[1], dim[0], dim[1]);
        exit(1);

        return 1;
      }

    *index = dim[1] * loc[0] + loc[1];

    return 0;
}


static int _mat_type_loc_f(size_t *dim, size_t index, size_t *loc)
{
    /*

        Get the (i, j) position of the index position in a 1-d array representation of a full matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (index >= dim[0] * dim[1])
      {
        printf("Index position %ld is out of bounds if the size of the matrix array is %ld.\n", index, dim[0] * dim[1]);
        exit(1);

        return 1;
      }

    loc[1] = index % dim[1];
    loc[0] = (index - loc[1]) / dim[1];

    return 0;
}


static bool _mat_type_comp_f(void *matrix, size_t *loc, double tol)
{
    /*

        Comparison function to determine if elements at loc in matrix satisfy requirements to be full.

        Always true.

   */

    (void) matrix;
    (void) loc;
    (void) tol;

    return true;
}


static int _mat_type_rbounds_f(size_t *dim, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds for a column of a full matrix

    */

    (void) col;

    rbounds[0] = 0;
    rbounds[1] = dim[0];

    return 0;
}


static int _mat_type_cbounds_f(size_t *dim, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds for a row of a full matrix

    */

    (void) row;

    cbounds[0] = 0;
    cbounds[1] = dim[1];

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_type_index_d(size_t *dim, size_t *loc, size_t *index)
{
    /*

        Get the index of the (i,j) position (loc) for a 1-d array representation of a diagonal matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (loc[0] >= dim[0] || loc[1] >= dim[1])
      {
        printf("Location (%ld,%ld) is out of bounds if the dimensions are (%ld,%ld).\n", loc[0], loc[1], dim[0], dim[1]);
        exit(1);

        return 1;
      }

    if (loc[0] == loc[1])
      {
        *index = *loc;

        return 0;
      }

    return -1;
}


static int _mat_type_loc_d(size_t *dim, size_t index, size_t *loc)
{
    /*

        Get the (i, j) position of the index position in a 1-d array representation of a diagonal matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (index >= dim[0])
      {
        printf("Index position %ld is out of bounds if the size of the matrix array is %ld.\n", index, dim[0]);
        exit(1);

        return 1;
      }

    loc[0] = index;
    loc[1] = index;

    return 0;
}


static bool _mat_type_comp_d(void *matrix, size_t *loc, double tol)
{
    /*

        Comparison function to determine if elements at loc in matrix satisfy requirements to be diagonal.

        If loc is NULL, only check if the matrix' type satisfies requiremenets to be diagonal

    */

    mat_t *mat = (mat_t*) matrix;

    /* Non-square matrices can never be diagonal */
    if (mat -> dim[0] != mat -> dim[1])
        return false;

    /* Matrix already has type 'd' or 'null' */
    if (mat -> mtype == &_matTypeNull || mat -> mtype == &_matTypeD)
        return true;

    /* loc is NULL, ie only check for types */
    if (loc == NULL)
        return false;

    /* Diagonal terms must not be checked */
    if (loc[0] == loc[1])
        return true;

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Off-diagonal terms must be close to zero */
        if (fabs(mat_get_value(mat, loc)) > tol)
            return false;
      }

    /* Block Matrix */
    else
      {
        /* Off-diagonal matrices must be empty */
        return _mat_type_comp_null(mat, loc, tol);
      }

    return true;
}


static int _mat_type_rbounds_d(size_t *dim, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds for a column of a diagonal matrix

    */

    (void) dim;

    rbounds[0] = col;
    rbounds[1] = col + 1;

    return 0;
}


static int _mat_type_cbounds_d(size_t *dim, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds for a row of a diagonal matrix

    */

    (void) dim;

    cbounds[0] = row;
    cbounds[1] = row + 1;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_type_index_s(size_t *dim, size_t *loc, size_t *index)
{
    /*

        Get the index of the (i,j) position (loc) for a 1-d array representation of a symmetric matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (loc[0] >= dim[0] || loc[1] >= dim[1])
      {
        printf("Location (%ld,%ld) is out of bounds if the dimensions are (%ld,%ld).\n", loc[0], loc[1], dim[0], dim[1]);
        exit(1);

        return 1;
      }

    if (loc[0] <= loc[1])
        *index = dim[1] * loc[0] - loc[0] * (loc[0] + 1) / 2 + loc[1];

    else
        *index = dim[0] * loc[1] - loc[1] * (loc[1] + 1) / 2 + loc[0];

    return 0;
}


static int _mat_type_loc_s(size_t *dim, size_t index, size_t *loc)
{
    /*

        Get the (i, j) position of the index position in a 1-d array representation of a symmetric matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (index >= dim[0] * (dim[0] + 1) / 2)
      {
        printf("Index position %ld is out of bounds if the size of the matrix array is %ld.\n", index, dim[0] * (dim[0] + 1) / 2);
        exit(1);

        return 1;
      }

    loc[0] = 0;
    loc[1] = index;

    if (index < dim[1])
      {
        return 0;
      }

    do
      {
        loc[1] -= dim[1] - loc[0];
        loc[0]++;
      }

    while (loc[1] >= dim[1] - loc[0]);

    loc[1] += loc[0];

    return 0;
}


static bool _mat_type_comp_s(void *matrix, size_t *loc, double tol)
{
    /*

        Comparison function to determine if elements at loc in matrix satisfy requirements to be symmetric.

        If loc is NULL, only check if matrix' type satisfies requirements to be symmetric.

    */

    mat_t *mat = (mat_t*) matrix;

    /* Non-square matrices can never be symmetric */
    if (mat -> dim[0] != mat -> dim[1])
        return false;

    /* Matrix already has type 'd', 's' or 'null' */
    if (mat -> mtype == &_matTypeNull || mat -> mtype == &_matTypeD || mat -> mtype == &_matTypeS)
        return true;

    /* loc is NULL, ie only check for types */
    if (loc == NULL)
        return false;

    /* Diagonal terms are always symmetric */
    if (loc[0] == loc[1])
        return true;

    /* Tranposed location */
    size_t locTp[2] = {loc[1], loc[0]};

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Off-diagonal terms must be equal */
        if (fabs(mat_get_value(mat, loc) - mat_get_value(mat, locTp)) > tol)
            return false;
      }

    /* Block Matrix */
    else
      {
        mat_t *matBlock = mat_mget_block(mat, loc);
        mat_t *matBlockTp = mat_ftp(mat_fget_block(mat, locTp), NULL);

        if (!mat_comp_matrix(matBlock, matBlockTp, tol))
          {
            matBlockTp = mat_ffree(matBlockTp);
            return false;
          }

        matBlockTp = mat_ffree(matBlockTp);
      }

    return true;
}


static int _mat_type_rbounds_s(size_t *dim, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds for a column of a symmetric matrix

    */

    (void) col;

    rbounds[0] = 0;
    rbounds[1] = dim[0];

    return 0;
}


static int _mat_type_cbounds_s(size_t *dim, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds for a row of a symmetric matrix

    */

    (void) row;

    cbounds[0] = 0;
    cbounds[1] = dim[1];

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_type_index_ut(size_t *dim, size_t *loc, size_t *index)
{
    /*

        Get the index of the (i,j) position (loc) for a 1-d array representation of an upper triangular matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (loc[0] >= dim[0] || loc[1] >= dim[1])
      {
        printf("Location (%ld,%ld) is out of bounds if the dimensions are (%ld,%ld).\n", loc[0], loc[1], dim[0], dim[1]);
        exit(1);

        return 1;
      }

    if (loc[0] <= loc[1])
      {
        *index = dim[1] * loc[0] - loc[0] * (loc[0] + 1) / 2 + loc[1];

        return 0;
      }

    return -1;
}


static int _mat_type_loc_ut(size_t *dim, size_t index, size_t *loc)
{
    /*

        Get the (i, j) position of the index position in a 1-d array representation of an upper triangular matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (index >= dim[0] * (dim[0] + 1) / 2)
      {
        printf("Index position %ld is out of bounds if the size of the matrix array is %ld.\n", index, dim[0] * (dim[0] + 1) / 2);
        exit(1);

        return 1;
      }

    loc[0] = 0;
    loc[1] = index;

    if (index < dim[1])
      {
        return 0;
      }

    do
      {
        loc[1] -= dim[1] - loc[0];
        loc[0]++;
      }

    while (loc[1] >= dim[1] - loc[0]);

    loc[1] += loc[0];

    return 0;
}


static bool _mat_type_comp_ut(void *matrix, size_t *loc, double tol)
{
    /*

        Comparison function to determine if elements at loc in matrix satisfy requirements to be upper triangular.

        If loc is NULL, only check if matrix' type satisfies requirements to be upper triangular.

    */

    mat_t *mat = (mat_t*) matrix;

    /* Non-square matrices can never be upper triangular */
    if (mat -> dim[0] != mat -> dim[1])
        return false;

    /* Matrix already has type 'null', 'd' or 'ut' */
    if (mat -> mtype == &_matTypeNull || mat -> mtype == &_matTypeD || mat -> mtype == &_matTypeUt)
        return true;

    /* loc is NULL, ie only check for types */
    if (loc == NULL)
        return false;

    /* Upper-triangular terms can be skipped */
    if (loc[0] <= loc[1])
        return true;

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Lower triangular terms must be zero */
        if (fabs(mat_get_value(mat, loc)) > tol)
            return false;
      }

    /* Block Matrix */
    else
      {
        /* Lower triangular matrices must be empty */
        return _mat_type_comp_null(mat, loc, tol);
      }

    return true;
}


static int _mat_type_rbounds_ut(size_t *dim, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds for a column of an upper triangular matrix

    */

    (void) dim;

    rbounds[0] = 0;
    rbounds[1] = col + 1;

    return 0;
}


static int _mat_type_cbounds_ut(size_t *dim, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds for a row of an upper triangular matrix

    */

    cbounds[0] = row;
    cbounds[1] = dim[1];

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_type_index_lt(size_t *dim, size_t *loc, size_t *index)
{
    /*

        Get the index of the (i,j) position (loc) for a 1-d array representation of a symmetric matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (loc[0] >= dim[0] || loc[1] >= dim[1])
      {
        printf("Location (%ld,%ld) is out of bounds if the dimensions are (%ld,%ld).\n", loc[0], loc[1], dim[0], dim[1]);
        exit(1);

        return 1;
      }

    if (loc[0] >= loc[1])
      {
        *index = dim[0] * loc[1] - loc[1] * (loc[1] + 1) / 2 + loc[0];

        return 0;
      }

    return -1;
}


static int _mat_type_loc_lt(size_t *dim, size_t index, size_t *loc)
{
    /*

        Get the (i, j) position of the index position in a 1-d array representation of a lower triangular matrix of dimensions dim.

    */

    /* Make sure location is within dimension bounds */
    if (index >= dim[0] * (dim[0] + 1) / 2)
      {
        printf("Index position %ld is out of bounds if the size of the matrix array is %ld.\n", index, dim[0] * (dim[0] + 1) / 2);
        exit(1);

        return 1;
      }

    loc[0] = index;
    loc[1] = 0;

    if (index < dim[0])
      {
        return 0;
      }

    do
      {
        loc[0] -= dim[0] - loc[1];
        loc[1]++;
      }

    while (loc[0] >= dim[0] - loc[1]);

    loc[0] += loc[1];

    return 0;
}


static bool _mat_type_comp_lt(void *matrix, size_t *loc, double tol)
{
    /*

        Comparison function to determine if elements at loc in matrix satisfy requirements to be lower triangular.

        If loc is NULL, only check if matrix' type satisfies requirements to be lower triangular.

    */

    mat_t *mat = (mat_t*) matrix;

    /* Non-square matrices can never be lower triangular */
    if (mat -> dim[0] != mat -> dim[1])
        return false;

    /* Matrix already has type 'null', 'd' or 'lt' */
    if (mat -> mtype == &_matTypeNull || mat -> mtype == &_matTypeD || mat -> mtype == &_matTypeLt)
        return true;

    /* loc is NULL, ie only check for types */
    if (loc == NULL)
        return false;

    /* Lower-triangular terms can be skipped */
    if (loc[0] >= loc[1])
        return true;

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Upper triangular terms must be zero */
        if (fabs(mat_get_value(mat, loc)) > tol)
            return false;
      }

    /* Block Matrix */
    else
      {
        /* Upper triangular matrices must be empty */
        return _mat_type_comp_null(mat, loc, tol);
      }

    return true;
}


static int _mat_type_rbounds_lt(size_t *dim, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds for a column of a lower triangular matrix

    */

    rbounds[0] = col;
    rbounds[1] = dim[0];

    return 0;
}


static int _mat_type_cbounds_lt(size_t *dim, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds for a row of a lower triangular matrix

    */

    (void) dim;

    cbounds[0] = 0;
    cbounds[1] = row + 1;

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     MATRIX LOCATION TRANSFORMATIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/**  Location Tranformations  **/

#ifndef __MAT_NUMBER_OF_LOCATION_TRAFOS__
#define __MAT_NUMBER_OF_LOCATION_TRAFOS__ 1
#endif

/* Number of trafos */
const size_t _matNumLocTrfs = __MAT_NUMBER_OF_LOCATION_TRAFOS__;

/* Trafo ID's */
const char *_matLocTrfsID[__MAT_NUMBER_OF_LOCATION_TRAFOS__] = {"tp"};


/* Transpose Transformation */

/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_ltrf_loc_tp(void *matrix, size_t *loc);
static int _mat_ltrf_dim_tp(void *matrix, size_t *dim);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */

const mat_ltrf_t _matLocTrfTp = {"tp", _mat_ltrf_loc_tp, _mat_ltrf_dim_tp};


/* All transformations */
const mat_ltrf_t _matLocTrfs[__MAT_NUMBER_OF_LOCATION_TRAFOS__] = {_matLocTrfTp};


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_ltrf_t *mat_ltrf_get_struct(const char *id)
{
    /*

        Find the location transformation with given id

    */

    for (size_t i = 0; i < _matNumLocTrfs; i++)
      {
        if (!strcmp(_matLocTrfsID[i], id))
          {
            return (mat_ltrf_t*) &_matLocTrfs[i];
          }
      }

    /* Did not find the correct location trafo */
    printf("Could not find the matrix location trafo with trafo ID '%s'.\n", id);
    exit(1);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------   Matrix Location Transformations   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_ltrf_loc_tp(void *matrix, size_t *loc)
{
    /*

        Transpose location

    */

    (void) matrix;

    misc_swap(&loc[0], &loc[1], "ld");

    return 0;
}


static int _mat_ltrf_dim_tp(void *matrix, size_t *dim)
{
    /*

        Transpose dimensions

    */

    (void) matrix;

    misc_swap(&dim[0], &dim[1], "ld");

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   New / Free / Cp   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_ltrfs_t *mat_ltrfs_new(void)
{
    /*

        Create a new mat_ltrfs_t struct

    */

    mat_ltrfs_t *mltrfs = malloc(sizeof(mat_ltrfs_t));

    mltrfs -> size = 0;
    mltrfs -> mltrf = NULL;

    return mltrfs;
}


mat_ltrfs_t *mat_ltrfs_free(mat_ltrfs_t *mltrfs)
{
    /*

        Create a new mat_ltrfs_t struct

    */

    /* Check for NULL */
    if (mltrfs == NULL)
        return NULL;

    /* Free contents */
    free(mltrfs -> mltrf);

    /* Free mltrfs itself */
    free(mltrfs);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


mat_ltrfs_t *mat_ltrfs_cp(mat_ltrfs_t *mltrfs)
{
    /*

        Create a new mat_ltrfs_t struct

    */

    mat_ltrfs_t *mltrfsCp = mat_ltrfs_new();

    mltrfsCp -> size = mltrfs -> size;
    mltrfsCp -> mltrf = malloc(sizeof(mat_ltrf_t*) * mltrfs -> size);

    for (size_t i = 0; i < mltrfs -> size; i++)
        mltrfsCp -> mltrf[i] = mltrfs -> mltrf[i];

    return mltrfsCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Add a Transformation   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_ltrfs_add(mat_ltrfs_t *mltrfs, const char *id)
{
    /*

        Add a location trafo to mltrfs

        TODO: Add possibility of trafos cancelling previous ones

    */

    mat_ltrf_t *mltrf = mat_ltrf_get_struct(id);

    mltrfs -> size += 1;
    mltrfs -> mltrf = realloc(mltrfs -> mltrf, sizeof(mat_ltrf_t*) * mltrfs -> size);

    mltrfs -> mltrf[mltrfs -> size - 1] = mltrf;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Block Matrix Insertion   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static void *_mat_blockrep_id(void *matrix, size_t *loc)
{
    /*

        Matrix' representation in parent's matrix is itself

    */

    (void) loc;

    return matrix;
}


/*  ------------------------------------------------------------------------------------------------------  */


static void *_mat_blockrep_s(void *matrix, size_t *loc)
{
    /*

        Matrix' representation in parent's matrix is itself in the upper triangular half, and its transpose in the true lower triangular half.

    */

    if (loc[0] > loc[1])
      {
        matrix =  mat_ftp((mat_t*) matrix, NULL);
      }

    return matrix;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     MATRIX STRUCT     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   New / Free Struct   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_fnew_null(size_t *dim, bool block);
static mat_t *_mat_new_null(size_t *dim, bool block);

static mat_t *_mat_fnew_f(size_t *dim, bool block);
static mat_t *_mat_new_f(size_t *dim, bool block);

static mat_t *_mat_fnew_d(size_t *dim, bool block);
static mat_t *_mat_new_d(size_t *dim, bool block);

static mat_t *_mat_fnew_s(size_t *dim, bool block);
static mat_t *_mat_new_s(size_t *dim, bool block);

static mat_t *_mat_fnew_ut(size_t *dim, bool block);
static mat_t *_mat_new_ut(size_t *dim, bool block);

static mat_t *_mat_fnew_lt(size_t *dim, bool block);
static mat_t *_mat_new_lt(size_t *dim, bool block);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_new(const char *type, size_t *dim, bool block)
{
    /*

        Create a new mat_t struct

    */


    mat_t *mat;

    /* Null (empty) Matrix */
    if (!strcmp(type, "null"))
      {
        mat = _mat_new_null(dim, block);

        return mat;
      }

    /* Full Matrix */
    if (!strcmp(type, "full") || !strcmp(type, "f"))
      {
        mat = _mat_new_f(dim, block);

        return mat;
      }

    /* Diagonal Matrix */
    if (!strcmp(type, "diagonal") || !strcmp(type, "d"))
      {
        mat = _mat_new_d(dim, block);

        return mat;
      }

    /* Symmetric Matrix */
    if (!strcmp(type, "symmetric") || !strcmp(type, "s"))
      {
        mat = _mat_new_s(dim, block);

        return mat;
      }

    /* Upper Triangular Matrix */
    if (!strcmp(type, "upper-triangular") || !strcmp(type, "ut"))
      {
        mat = _mat_new_ut(dim, block);

        return mat;
      }

    /* Lower Triangular Matrix */
    if (!strcmp(type, "lower-triangular") || !strcmp(type, "lt"))
      {
        mat = _mat_new_lt(dim, block);

        return mat;
      }


    /* Wrong input */

    printf("The matrix type '%s' does not match any of the expected types:\n", type);

    printf("  -  'null'\n");
    printf("  -  'full' or 'f'\n");
    printf("  -  'diagonal' or 'd'\n");
    printf("  -  'symmetric' or 's'\n");
    printf("  -  'upper-triangular' or 'ut'\n");
    printf("  -  'lower-triangular' or 'lt'\n");

    exit(1);

    return NULL;
}


mat_t *mat_fnew(const char *type, size_t *dim, bool block)
{
    /*

        Create a new mat_t struct

        Note: This function does not allocate memory to mat -> matrix

    */


    mat_t *mat;

    /* Null (empty) Matrix */
    if (!strcmp(type, "null"))
      {
        mat = _mat_fnew_null(dim, block);

        return mat;
      }

    /* Full Matrix */
    if (!strcmp(type, "full") || !strcmp(type, "f"))
      {
        mat = _mat_fnew_f(dim, block);

        return mat;
      }

    /* Diagonal Matrix */
    if (!strcmp(type, "diagonal") || !strcmp(type, "d"))
      {
        mat = _mat_fnew_d(dim, block);

        return mat;
      }

    /* Symmetric Matrix */
    if (!strcmp(type, "symmetric") || !strcmp(type, "s"))
      {
        mat = _mat_fnew_s(dim, block);

        return mat;
      }

    /* Upper Triangular Matrix */
    if (!strcmp(type, "upper-triangular") || !strcmp(type, "ut"))
      {
        mat = _mat_fnew_ut(dim, block);

        return mat;
      }

    /* Lower Triangular Matrix */
    if (!strcmp(type, "lower-triangular") || !strcmp(type, "lt"))
      {
        mat = _mat_fnew_lt(dim, block);

        return mat;
      }


    /* Wrong input */

    printf("The matrix type '%s' does not match any of the expected types:\n", type);

    printf("  -  'null'\n");
    printf("  -  'full' or 'f'\n");
    printf("  -  'diagonal' or 'd'\n");
    printf("  -  'symmetric' or 's'\n");
    printf("  -  'upper-triangular' or 'ut'\n");
    printf("  -  'lower-triangular' or 'lt'\n");

    exit(1);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_new_null(size_t *dim, bool block)
{
    /*

        Create a new null matrix

    */

    /* Matrix struct */
    mat_t *mat = _mat_fnew_null(dim, block);

    return mat;
}


static mat_t *_mat_fnew_null(size_t *dim, bool block)
{
    /*

        Create a new null matrix.

        Note: Does not allocate to mat -> matrix

    */

    /* Matrix struct */
    mat_t *mat = malloc(sizeof(mat_t));

    /* Matrix file */
    mat -> mfile = NULL;

    /* Matrix label */
    mat -> label = NULL;

    /* Matrix type */
    mat -> mtype = (mat_type_t*) &_matTypeNull;

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_new();

    /* Set dimensions */
    mat -> dim = malloc(sizeof(size_t) * 3);

    mat -> dim[0] = dim[0];
    mat -> dim[1] = dim[1];
    mat -> dim[2] = 0;

    /* Block Matrix */
    if (block)
      {
        printf("Cannot create a null matrix with block structure.\n");
        exit(1);

        return NULL;
      }

    /* Ordinary Matrix */
    else
      {
        /* No block variables */
        mat -> mblock = NULL;
      }

    /* Initialise matrix to NULL */
    mat -> matrix = NULL;

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_new_f(size_t *dim, bool block)
{
    /*

        Create a new full matrix

    */

    /* Matrix struct */
    mat_t *mat = _mat_fnew_f(dim, block);

    /* Block Matrix */
    if (block)
      {
        /* Matrix */
        mat -> matrix = malloc(sizeof(mat_t*) * mat -> dim[2]);

        for (size_t n = 0; n < mat -> dim[2]; n++)
          {
            /* Set block to NULL */
            ((mat_t**) mat -> matrix)[n] = NULL;
          }
      }

    /* Ordinary Matrix */
    else
      {
        /* Initialise matrix elements to zero */
        mat -> matrix = calloc(mat -> dim[2], sizeof(double));
      }

    return mat;
}


static mat_t *_mat_fnew_f(size_t *dim, bool block)
{
    /*

        Create a new full matrix

        Note: This function does not allocate memory to mat -> matrix

    */

    /* Matrix struct */
    mat_t *mat = malloc(sizeof(mat_t));

    /* Matrix file */
    mat -> mfile = NULL;

    /* Matrix label */
    mat -> label = NULL;

    /* Matrix type */
    mat -> mtype = (mat_type_t*) &_matTypeF;

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_new();

    /* Set dimensions */
    mat -> dim = malloc(sizeof(size_t) * 3);

    mat -> dim[0] = dim[0];
    mat -> dim[1] = dim[1];
    mat -> dim[2] = dim[0] * dim[1];

    /* Block Matrix */
    if (block)
      {
        /* Block variables */
        mat_block_t *block = malloc(sizeof(mat_block_t));

        /* Block dimensions */
        block -> bdim = malloc(sizeof(size_t*) * 2);
        block -> bdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> bfdim = malloc(sizeof(size_t*) * 2);
        block -> bfdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bfdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> fdim = calloc(2, sizeof(size_t));

        /* Block representations */
        block -> rep = _mat_blockrep_id;

        /* Block struct */
        mat -> mblock = block;
      }

    /* Ordinary Matrix */
    else
      {
        /* No block variables */
        mat -> mblock = NULL;
      }

    /* Initialise matrix to NULL */
    mat -> matrix = NULL;

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_new_d(size_t *dim, bool block)
{
    /*

        Create a new diagonal matrix

    */

    /* Matrix struct */
    mat_t *mat = _mat_fnew_d(dim, block);

    /* Block Matrix */
    if (block)
      {
        /* Matrix */
        mat -> matrix = malloc(sizeof(mat_t*) * mat -> dim[2]);

        for (size_t n = 0; n < mat -> dim[2]; n++)
          {
            /* Set block to NULL */
            ((mat_t**) mat -> matrix)[n] = NULL;
          }
      }

    /* Ordinary Matrix */
    else
      {
        /* Initialise matrix to zero */
        mat -> matrix = calloc(mat -> dim[2], sizeof(double));
      }

    return mat;
}


static mat_t *_mat_fnew_d(size_t *dim, bool block)
{
    /*

        Create a new diagonal matrix

        Note: This function does not allocate memory to mat -> matrix

    */

    /* Must have dim1 = dim2 */
    if (dim[0] != dim[1])
      {
        printf("Cannot create a diagonal matrix that is not square.\n");

        exit(1);

        return NULL;
      }

    /* Matrix struct */
    mat_t *mat = malloc(sizeof(mat_t));

    /* Matrix file */
    mat -> mfile = NULL;

    /* Matrix label */
    mat -> label = NULL;

    /* Matrix type */
    mat -> mtype = (mat_type_t*) &_matTypeD;

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_new();

    /* Set dimensions */
    mat -> dim = malloc(sizeof(size_t) * 3);

    mat -> dim[0] = *dim;
    mat -> dim[1] = *dim;
    mat -> dim[2] = *dim;

    /* Block Matrix */
    if (block)
      {
        /* Block variables */
        mat_block_t *block = malloc(sizeof(mat_block_t));

        /* Block dimensions */
        block -> bdim = malloc(sizeof(size_t*) * 2);
        block -> bdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> bfdim = malloc(sizeof(size_t*) * 2);
        block -> bfdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bfdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> fdim = calloc(2, sizeof(size_t));

        /* Block representation */
        block -> rep = _mat_blockrep_id;

        /* Block struct */
        mat -> mblock = block;
      }

    /* Ordinary Matrix */
    else
      {
        /* No block variables */
        mat -> mblock = NULL;
      }

    /* Initialise matrix to NULL */
    mat -> matrix = NULL;

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_new_s(size_t *dim, bool block)
{
    /*

        Create a new symmetric matrix

    */

    /* Matrix struct */
    mat_t *mat = _mat_fnew_s(dim, block);

    /* Block Matrix */
    if (block)
      {
        /* Matrix */
        mat -> matrix = malloc(sizeof(mat_t*) * mat -> dim[2]);

        for (size_t n = 0; n < mat -> dim[2]; n++)
          {
            /* Set block to NULL */
            ((mat_t**) mat -> matrix)[n] = NULL;
          }
      }

    /* Ordinary Matrix */
    else
      {
        /* Initialise matrix to zero */
        mat -> matrix = calloc(mat -> dim[2], sizeof(double));
      }

    return mat;
}


static mat_t *_mat_fnew_s(size_t *dim, bool block)
{
    /*

        Create a new symmetric matrix

        Note: This function does not allocate memory to mat -> matrix

    */

    /* Must have dim1 = dim2 */
    if (dim[0] != dim[1])
      {
        printf("Cannot create a symmetric matrix that is not square.\n");

        exit(1);

        return NULL;
      }

    /* Matrix struct */
    mat_t *mat = malloc(sizeof(mat_t));

    /* Matrix file */
    mat -> mfile = NULL;

    /* Matrix label */
    mat -> label = NULL;

    /* Matrix type */
    mat -> mtype = (mat_type_t*) &_matTypeS;

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_new();

    /* Set dimensions */
    mat -> dim = malloc(sizeof(size_t) * 3);

    mat -> dim[0] = *dim;
    mat -> dim[1] = *dim;
    mat -> dim[2] = *dim * (*dim + 1) / 2;

    /* Block Matrix */
    if (block)
      {
        /* Block variables */
        mat_block_t *block = malloc(sizeof(mat_block_t));

        /* Block dimensions */
        block -> bdim = malloc(sizeof(size_t*) * 2);
        block -> bdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> bfdim = malloc(sizeof(size_t*) * 2);
        block -> bfdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bfdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> fdim = calloc(2, sizeof(size_t));

        /* Block representation */
        block -> rep = _mat_blockrep_s;

        /* Block struct */
        mat -> mblock = block;
      }

    /* Ordinary Matrix */
    else
      {
        /* No block variables */
        mat -> mblock = NULL;
      }

    /* Initialise matrix to NULL */
    mat -> matrix = NULL;

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_new_ut(size_t *dim, bool block)
{
    /*

        Create a new upper triangular matrix

    */

    /* Matrix struct */
    mat_t *mat = _mat_fnew_ut(dim, block);

    /* Block Matrix */
    if (block)
      {
        /* Matrix */
        mat -> matrix = malloc(sizeof(mat_t*) * mat -> dim[2]);

        for (size_t n = 0; n < mat -> dim[2]; n++)
          {
            /* Set block to NULL */
            ((mat_t**) mat -> matrix)[n] = NULL;
          }
      }

    else
      {
        /* Initialise matrix to zero */
        mat -> matrix = calloc(mat -> dim[2], sizeof(double));
      }

    return mat;
}


static mat_t *_mat_fnew_ut(size_t *dim, bool block)
{
    /*

        Create a new upper triangular matrix

        Note: This function does not allocate memory to mat -> matrix

    */

    /* Must have dim1 = dim2 */
    if (dim[0] != dim[1])
      {
        printf("Cannot create an upper triangular matrix that is not square.\n");

        exit(1);

        return NULL;
      }

    /* Matrix struct */
    mat_t *mat = malloc(sizeof(mat_t));

    /* Matrix file */
    mat -> mfile = NULL;

    /* Matrix label */
    mat -> label = NULL;

    /* Matrix type */
    mat -> mtype = (mat_type_t*) &_matTypeUt;

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_new();

    /* Set dimensions */
    mat -> dim = malloc(sizeof(size_t) * 3);

    mat -> dim[0] = *dim;
    mat -> dim[1] = *dim;
    mat -> dim[2] = *dim * (*dim + 1) / 2;

    /* Block Matrix */
    if (block)
      {
        /* Block variables */
        mat_block_t *block = malloc(sizeof(mat_block_t));

        /* Block dimensions */
        block -> bdim = malloc(sizeof(size_t*) * 2);
        block -> bdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> bfdim = malloc(sizeof(size_t*) * 2);
        block -> bfdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bfdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> fdim = calloc(2, sizeof(size_t));

        /* Block representation */
        block -> rep = _mat_blockrep_id;

        /* Block struct */
        mat -> mblock = block;
      }

    /* Ordinary Matrix */
    else
      {
        /* No block varriables */
        mat -> mblock = NULL;
      }

    /* Initialise matrix to NULL */
    mat -> matrix = NULL;

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_new_lt(size_t *dim, bool block)
{
    /*

        Create a new lower triangular matrix

    */

    /* Matrix struct */
    mat_t *mat = _mat_fnew_lt(dim, block);

    /* Block Matrix */
    if (block)
      {
        /* Matrix */
        mat -> matrix = malloc(sizeof(mat_t*) * mat -> dim[2]);

        for (size_t n = 0; n < mat -> dim[2]; n++)
          {
            /* Set block to NULL */
            ((mat_t**) mat -> matrix)[n] = NULL;
          }
      }

    else
      {
        /* Initialise matrix to zero */
        mat -> matrix = calloc(mat -> dim[2], sizeof(double));
      }

    return mat;
}


static mat_t *_mat_fnew_lt(size_t *dim, bool block)
{
    /*

        Create a new lower triangular matrix

    */

    /* Must have dim1 = dim2 */
    if (dim[0] != dim[1])
      {
        printf("Cannot create a lower triangular matrix that is not square.\n");

        exit(1);

        return NULL;
      }

    /* Matrix struct */
    mat_t *mat = malloc(sizeof(mat_t));

    /* Matrix file */
    mat -> mfile = NULL;

    /* Matrix label */
    mat -> label = NULL;

    /* Matrix type */
    mat -> mtype = (mat_type_t*) &_matTypeLt;

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_new();

    /* Set dimensions */
    mat -> dim = malloc(sizeof(size_t) * 3);

    mat -> dim[0] = *dim;
    mat -> dim[1] = *dim;
    mat -> dim[2] = *dim * (*dim + 1) / 2;

    /* Block Matrix */
    if (block)
      {
        /* Block variables */
        mat_block_t *block = malloc(sizeof(mat_block_t));

        /* Block dimensions */
        block -> bdim = malloc(sizeof(size_t*) * 2);
        block -> bdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> bfdim = malloc(sizeof(size_t*) * 2);
        block -> bfdim[0] = calloc(mat -> dim[0], sizeof(size_t));
        block -> bfdim[1] = calloc(mat -> dim[1], sizeof(size_t));

        block -> fdim = calloc(2, sizeof(size_t));

        /* Block representation */
        block -> rep = _mat_blockrep_id;

        /* Block struct */
        mat -> mblock = block;
      }

    /* Ordinary Matrix */
    else
      {
        /* No block variables */
        mat -> mblock = NULL;
      }

    /* Initialise matrix to NULL */
    mat -> matrix = NULL;

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_new_id(mat_t *mat)
{
    /*

        Create a new identity matrix while copying mat's block structure

    */

    if (mat -> dim[0] != mat -> dim[1])
      {
        printf("Cannot create a non-square identity matrix.\n");
        exit(1);

        return NULL;
      }

    mat_t *matId = mat_new(_matTypeD.id, mat -> dim, mat -> mblock != NULL);

    /* Ordinary Matrix */
    if (matId -> mblock == NULL)
      {
        for (size_t i = 0; i < *(matId -> dim); i++)
            ((double*) matId -> matrix)[i] = 1.;
      }

    /* Block Matrix */
    else
      {
        for (size_t i = 0; i < matId -> dim[2]; i++)
          {
            size_t loc[2] = {i, i};

            mat_t *matBlock = mat_fget_block(mat, loc);
            mat_t *matBlockId = mat_new_id(matBlock);

            mat_fset_block(matId, loc, matBlockId);

            matBlock = mat_ffree(matBlock);
          }
      }

    return matId;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_free(mat_t *mat)
{
    /*

        Free a mat_t struct

    */

    /* If mat is NULL simply return NULL */
    if (mat == NULL)
        return NULL;

    /* File */
    mat -> mfile = _mat_file_free(mat -> mfile);

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_free(mat -> mltrfs);

    /* Label */
    free(mat -> label);

    /* Block matrix */
    if (mat -> mblock != NULL)
      {
        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            ((mat_t**) mat -> matrix)[i] = mat_free(((mat_t**) mat -> matrix)[i]);
          }

        for (size_t i = 0; i < 2; i++)
          {
            free(mat -> mblock -> bdim[i]);
            free(mat -> mblock -> bfdim[i]);
          }

        free(mat -> mblock -> bdim);
        free(mat -> mblock -> bfdim);
        free(mat -> mblock -> fdim);
      }

    free(mat -> matrix);
    free(mat -> dim);
    free(mat -> mblock);

    /* Free mat itself */
    free(mat);

    return NULL;
}


mat_t *mat_ffree(mat_t *mat)
{
    /*

        Free a mat_t struct.

        Note: Does not free mat -> matrix

    */

    /* If mat is NULL simply return NULL */
    if (mat == NULL)
        return NULL;

    /* File */
    mat -> mfile = _mat_file_free(mat -> mfile);

    /* Location transformations */
    mat -> mltrfs = mat_ltrfs_free(mat -> mltrfs);

    /* Label */
    free(mat -> label);

    /* Block matrix */
    if (mat -> mblock != NULL)
      {
        for (size_t i = 0; i < 2; i++)
          {
            free(mat -> mblock -> bdim[i]);
            free(mat -> mblock -> bfdim[i]);
          }

        free(mat -> mblock -> bdim);
        free(mat -> mblock -> bfdim);
        free(mat -> mblock -> fdim);
      }

    free(mat -> dim);
    free(mat -> mblock);

    /* Free mat itself */
    free(mat);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Submatrix   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_sub(mat_t *mat, size_t *indices1, size_t size1, size_t *indices2, size_t size2, mat_reduce_t *reduce)
{
    /*

        Get the submatrix of mat at given indices

    */

    /* Check sizes */
    if (size1 > mat -> dim[0] || size2 > mat -> dim[1])
      {
        printf("Cannot get a submatrix of dimensions (%ld,%ld) from a matrix of dimensions (%ld,%ld).\n", size1, size2, mat -> dim[0], mat -> dim[1]);
        exit(1);

        return NULL;
      }

    /* Check indices */
    for (size_t i = 0; i < size1; i++)
      {
        if (indices1[i] >= mat -> dim[0])
          {
            printf("Cannot get a submatrix of dimensions (%ld,%ld) from a matrix of dimensions (%ld,%ld) at row %ld.\n", size1, size2, mat -> dim[0], mat -> dim[1], indices1[i]);
            exit(1);

            return NULL;
          }
      }

    for (size_t i = 0; i < size2; i++)
      {
        if (indices2[i] >= mat -> dim[1])
          {
            printf("Cannot get a submatrix of dimensions (%ld,%ld) from a matrix of dimensions (%ld,%ld) at column %ld.\n", size1, size2, mat -> dim[0], mat -> dim[1], indices2[i]);
            exit(1);

            return NULL;
          }
      }

    /* Submatrix dimensions */
    size_t dimSub[2] = {size1, size2};

    /* Submatrix */
    mat_t *matSub = mat_new(_matTypeF.id, dimSub, mat -> mblock != NULL);

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {

        for (size_t i = 0; i < size1; i++)
          {
            for (size_t j = 0; j < size2; j++)
              {
                size_t locSub[2] = {i, j};
                size_t loc[2] = {indices1[i], indices2[j]};

                mat_set_value(matSub, locSub, mat_get_value(mat, loc));
              }
          }

      }

    /* Block Matrix */
    else
      {

        for (size_t i = 0; i < size1; i++)
          {
            for (size_t j = 0; j < size2; j++)
              {
                size_t locSub[2] = {i, j};
                size_t loc[2] = {indices1[i], indices2[j]};

                mat_t *matBlock = mat_fget_block(mat, loc);

                mat_set_block(matSub, locSub, matBlock);

                matBlock = mat_ffree(matBlock);
              }
          }

      }

    /* Attempt to reduce the matrix */
    matSub = mat_reduce(matSub, reduce);

    return matSub;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Copy Matrix   -------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_cp_o(mat_t *matCp, mat_t *mat);
static int _mat_cp_b(mat_t *matCp, mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_cp(mat_t *mat)
{
    /*

        Deep copy a matrix (mat)

    */

    /* New matrix */
    mat_t *matCp = mat_new((char*) mat -> mtype -> id, mat -> dim, mat -> mblock != NULL);

    /* Copy the label */
    mat_set_label(matCp, mat -> label);

    /* Copy the file parameters */
    matCp -> mfile = _mat_file_cp(mat -> mfile);

    /* Copy the location trafos */
    for (size_t i = 0; i < mat -> mltrfs -> size; i++)
        mat_ltrfs_add(matCp -> mltrfs, mat -> mltrfs -> mltrf[i] -> id);

    /* Oridnary Matrix */
    if (matCp -> mblock == NULL)
      {
        _mat_cp_o(matCp, mat);
      }

    /* Block Matrix */
    else
      {
        _mat_cp_b(matCp, mat);
      }

    return matCp;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_cp_o(mat_t *matCp, mat_t *mat)
{
    /*

        Deep copy an ordinary matrix

    */

    #pragma omp parallel for schedule(dynamic) shared(mat, matCp)
    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds[2] = {0, 0};
        mat_get_cbounds(mat, i, cBounds);

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            /* Index of location in memory */
            size_t loc[2] = {i, j};
            size_t index = 0;
            int success = mat_mget_index(mat, loc, &index);

            /* Skip if an element does not exist */
            if (success == -1)
                continue;

            /* Copy the element */
            ((double*) matCp -> matrix)[index] = ((double*) mat -> matrix)[index];
          }
      }

    return 0;
}


static int _mat_cp_b(mat_t *matCp, mat_t *mat)
{
    /*

        Deep copy a block matrix

    */

    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Index of location in memory */
            size_t loc[2] = {i, j};
            size_t index = 0;
            int success = mat_mget_index(mat, loc, &index);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* Block to be copied */
            mat_t *matBlock = mat_get_block(mat, loc);

            /* Can simply insert NULL */
            if (matBlock == NULL)
              {
                ((mat_t**) matCp -> matrix)[index] = NULL;

                continue;
              }

            /* Copy and insert the block */
            mat_fset_block(matCp, loc, matBlock);
          }
      }

    free(visitedIndices);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_fcp(mat_t *mat)
{
    /*

        Copy mat_t struct

        Note: Deep copies everything, except for mat -> matrix

    */

    /* Check for NULL */
    if (mat == NULL)
        return NULL;

    /* New matrix */
    mat_t *matCp = mat_fnew((char*) mat -> mtype -> id, mat -> dim, mat -> mblock != NULL);

    /* (Deep) Copy the block dimensions */
    if (mat -> mblock != NULL)
      {
        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            matCp -> mblock -> bdim[0][i] = mat -> mblock -> bdim[0][i];
            matCp -> mblock -> bfdim[0][i] = mat -> mblock -> bfdim[0][i];
          }

        for (size_t i = 0; i < mat -> dim[1]; i++)
          {
            matCp -> mblock -> bdim[1][i] = mat -> mblock -> bdim[1][i];
            matCp -> mblock -> bfdim[1][i] = mat -> mblock -> bfdim[1][i];
          }

        matCp -> mblock -> fdim[0] = mat -> mblock -> fdim[0];
        matCp -> mblock -> fdim[1] = mat -> mblock -> fdim[1];
      }

    /* (Deep) Copy the label */
    mat_set_label(matCp, mat -> label);

    /* (Deep) Copy the file parameters */
    matCp -> mfile = _mat_file_cp(mat -> mfile);

    /* (Deep) Copy the location trafos */
    for (size_t i = 0; i < mat -> mltrfs -> size; i++)
        mat_ltrfs_add(matCp -> mltrfs, mat -> mltrfs -> mltrf[i] -> id);

    /* (Shallow) Copy Matrix */
    matCp -> matrix = mat -> matrix;

    return matCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_cp_struct(mat_t *mat, const char *id, int idDepth)
{
    /*

        Copy the structure of mat.

        If id is NULL and idDepth is 0 copy the matrix ids of mat, and if idDepth is not 0 set the ids to "f".

        If id is not NULL, and idDepth is positive, set the matrix ids to id up until a depth of idDepth - 1,
        with the deeper blocks' ids copied from mat. If id is not NULL, and idDepth is negative, set the matrix
        ids to id up until a depth of - (1 + idDepth), with the deeper blocks' ids set to "f". If id is not NULL,
        and idDepth is zero, set all blocks' ids to id.

        If id is 'null', the entire matrix is set to a 'null' type matrix.

    */

    /* Check for mat NULL */
    if (mat == NULL)
        return NULL;

    /* Id */
    const char *idMat = (id != NULL) ? id
                                     : (idDepth == 0) ? mat -> mtype -> id
                                                      : _matTypeF.id;

    /* Full dimensions */
    size_t fdim[2] = {0, 0};
    mat_get_fdim(mat, fdim);

    /* Block matrix flag */
    bool block = (strcmp(idMat, _matTypeNull.id)) ? mat -> mblock != NULL : false; // 'null' matrix must be ordinary

    /* Dimension */
    size_t dim[2];

    dim[0] = (block) ? mat -> dim[0] : fdim[0];
    dim[1] = (block) ? mat -> dim[1] : fdim[1];

    /* New matrix */
    mat_t *matCp = mat_new(idMat, dim, block);


    /* Ordinary Matrix */

    if (matCp -> mblock == NULL)
        return matCp;


    /* Block matrix */

    /* Depth and locations */
    size_t depth = 0;
    size_t **loc = malloc(sizeof(size_t*) * (depth + 1));

    for (size_t i = 0; i < depth + 1; i++)
      {
        loc[i] = malloc(sizeof(size_t) * 2);
      }

    /* Keep track of the visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < matCp -> dim[0]; i++)
      {
        for (size_t j = 0; j < matCp -> dim[1]; j++)
          {
            /* Location in matCp's matrix and index in matCp's memory */
            loc[depth][0] = i;
            loc[depth][1] = j;

            size_t indexMem = 0;
            int success = mat_mget_index(matCp, loc[depth], &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            mat_t *matBlock = mat_cp_struct_block(mat, depth, loc,
                                                    (id == NULL || idDepth == 1 || idDepth == -1) ? NULL : id,
                                                    (idDepth > 0 && id != NULL) ? idDepth - 1 : (idDepth < -1 && id != NULL) ? idDepth + 1 : idDepth); // Need to exlcude -1 since we cannot have idDepth = 0 in the next recursion to have type $f$ (could also pass "f" directly)

            mat_fset_block(matCp, loc[depth], matBlock);
          }
      }

    /* Free memory */
    for (size_t i = 0; i < depth + 1; i++)
        free(loc[i]);

    free(loc);
    free(visitedIndices);

    return matCp;
}


mat_t *mat_cp_structs(mats_t *mats, const char *id, int idDepth)
{
    /*

        Copy the similarities of the structures of mats.

        If id is NULL and idDepth is 0 copy the matrix ids of mats, and if idDepth is not 0 set the ids to "f".

        If id is not NULL, and idDepth is positive, set the matrix ids to id up until a depth of idDepth - 1,
        with the deeper blocks' ids copied from mat. If id is not NULL, and idDepth is negative, set the matrix
        ids to id up until a depth of - (1 + idDepth), with the deeper blocks' ids set to "f". If id is not NULL,
        and idDepth is zero, set all blocks' ids to id.

    */

    /* Check for mat NULL */
    for (size_t i = 0; i < mats -> size; i++)
      {
        if (mats_get_mat(mats, i) == NULL)
          {
            /* Exclude the i'th matrix */
            mats_t *matsCp = mats_new(0);

            for (size_t j = 0; j < mats -> size; j++)
              {
                if (i == j)
                    continue;

                mats_set_mat(matsCp, matsCp -> size, mats_get_mat(mats, j));
              }

            /* Get the matrix */
            mat_t *matCp = mat_cp_structs(matsCp, id, idDepth);

            /* Free memory */
            matsCp = mats_free(matsCp);

            return matCp;
          }
      }

    /* Check dimensions */
    mat_t *mat = mats_get_mat(mats, 0); // First matrix (arbitrary choice)

    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    for (size_t i = 1; i < mats -> size; i++)
      {
        size_t fdim_[2];
        mat_get_fdim(mats_get_mat(mats, i), fdim_);

        if (fdim[0] != fdim_[0] || fdim[1] != fdim_[1])
          {
            printf("Cannot copy the matrix structures of a collection of matrices when the total dimensions do not match.");
            exit(1);

            return NULL;
          }
      }

    /* Matrix id */
    const char *idMat = mat -> mtype -> id;

    for (size_t i = 1; i < mats -> size; i++)
      {
        /* Final matrix has full type if the matrices have different types */
        if (strcmp(idMat, mats_get_mat(mats, i) -> mtype -> id))
          {
            idMat = _matTypeF.id;

            break;
          }
      }

    idMat = (id != NULL) ? id
                         : (idDepth == 0) ? idMat
                                          : _matTypeF.id;


    /* Block matrix flag */
    bool block = mat -> mblock != NULL;

    for (size_t i = 1; i < mats -> size; i++)
      {
        /* Final matrix is ordinary if the matrices have different block structures */
        if (!mat_comp_bstruct(mat, mats_get_mat(mats, i)))
          {
            block = false;

            break;
          }
      }

    block = (id == NULL || strcmp(id, _matTypeNull.id)) ? block : false; // 'null' matrix must be ordinary

    /* Matrix dimensions */
    size_t dim[2];

    dim[0] = (block) ? mat -> dim[0] : fdim[0];
    dim[1] = (block) ? mat -> dim[1] : fdim[1];

    /* New matrix */
    mat_t *matCp = mat_new(idMat, dim, block);


    /* Ordinary Matrix */

    if (matCp -> mblock == NULL)
        return matCp;


    /* Block matrix */

    /* Depth and locations */
    size_t depth = 0;
    size_t **loc = malloc(sizeof(size_t*) * (depth + 1));

    for (size_t i = 0; i < depth + 1; i++)
      {
        loc[i] = malloc(sizeof(size_t) * 2);
      }

    /* Keep track of the visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < matCp -> dim[0]; i++)
      {
        for (size_t j = 0; j < matCp -> dim[1]; j++)
          {
            /* Location in matCp's matrix and index in matCp's memory */
            loc[depth][0] = i;
            loc[depth][1] = j;

            size_t indexMem = 0;
            int success = mat_mget_index(matCp, loc[depth], &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            mat_t *matBlock = mat_cp_structs_block(mats, depth, loc,
                                                    (id == NULL || idDepth - 1 == 0 || idDepth + 1 == 0) ? NULL : id,
                                                    (idDepth > 0 && id != NULL) ? idDepth - 1 : (idDepth < -1 && id != NULL) ? idDepth + 1 : idDepth);

            mat_fset_block(matCp, loc[depth], matBlock);
          }
      }

    /* Free memory */
    for (size_t i = 0; i < depth + 1; i++)
        free(loc[i]);

    free(loc);
    free(visitedIndices);

    return matCp;
}


mat_t *mat_cp_struct_rc(mat_t *mat1, mat_t *mat2)
{
    /*

        Copy the structure of mat

    */

    /* Check for mat NULL */
    if (mat1 == NULL || mat2 == NULL)
        return NULL;

    /* Check dimensions */
    size_t fdim1[2];
    mat_get_fdim(mat1, fdim1);

    size_t fdim2[2];
    mat_get_fdim(mat2, fdim2);

    if (fdim1[1] != fdim2[0])
      {
        printf("Cannot copy the row/column matrix structures of two matrices with total dimensions (%ld, %ld) and (%ld, %ld).", fdim1[0], fdim1[1], fdim2[0], fdim2[1]);
        exit(1);

        return NULL;
      }

    /* Matrix id */
    const char *idMat = _matTypeF.id;

    /* Block matrix flag */
    bool block = mat1 -> mblock != NULL && mat_comp_bstruct_cr(mat1, mat2);

    /* Dimensions */
    size_t dim[2] = {mat1 -> dim[0], mat2 -> dim[1]};

    /* New matrix */
    mat_t *matCp = mat_new(idMat, dim, block);


    /* Ordinary Matrix */

    if (matCp -> mblock == NULL)
        return matCp;


    /* Block matrix */

    /* Depth and locations */
    size_t depth = 0;
    size_t **loc = malloc(sizeof(size_t*) * (depth + 1));

    for (size_t i = 0; i < depth + 1; i++)
      {
        loc[i] = malloc(sizeof(size_t) * 2);
      }

    /* Keep track of the visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < matCp -> dim[0]; i++)
      {
        for (size_t j = 0; j < matCp -> dim[1]; j++)
          {
            /* Location in matCp's matrix and index in matCp's memory */
            loc[depth][0] = i;
            loc[depth][1] = j;

            size_t indexMem = 0;
            int success = mat_mget_index(matCp, loc[depth], &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            mat_fset_block(matCp, loc[depth], mat_cp_struct_block_rc(mat1, mat2, depth, loc));
          }
      }

    /* Free memory */
    for (size_t i = 0; i < depth + 1; i++)
        free(loc[i]);

    free(loc);
    free(visitedIndices);

    return matCp;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_cp_struct_block_dim(mat_t *mat, size_t depth, size_t depthTot, size_t **loc, size_t *locBlock, const char *id, size_t index, size_t *dim);
static int _mat_cp_struct_block_id(mat_t *mat, size_t depth, size_t depthTot, size_t **loc, const char **id);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_cp_struct_block(mat_t *mat, size_t depth, size_t **loc, const char *id, int idDepth)
{
    /*

        Copy the structure of mat's block at given location

        If id is NULL and idDepth is 0 copy the matrix ids of mat, and if idDepth is not 0 set the ids to "f".

        If id is not NULL, and idDepth is positive, set the matrix ids to id up until a depth of idDepth - 1,
        with the deeper blocks' ids copied from mat. If id is not NULL, and idDepth is negative, set the matrix
        ids to id up until a depth of - (1 + idDepth), with the deeper blocks' ids set to "f". If id is not NULL,
        and idDepth is zero, set all blocks' ids to id.

    */

    /* Check for mat NULL */
    if (mat == NULL)
        return NULL;

    /* Check for ordinary mat */
    if (mat -> mblock == NULL)
      {
        printf("Cannot copy mat's block structure if mat is an ordinary matrix.\n");
        exit(1);

        return NULL;
      }

    /* Check for depth */
    size_t depthTot = mat_get_depth(mat);

    if (depthTot < depth)
      {
        printf("Cannot copy a block in mat at depth %ld and location (%ld, %ld) if the max depth is %ld.\n", depth, loc[0][0], loc[0][1], depthTot);
        exit(1);

        return NULL;
      }

    /* Type of block */
    const char *idBlock = _matTypeF.id;
    _mat_cp_struct_block_id(mat, depth, depthTot, loc, &idBlock);

    idBlock = (id != NULL) ? id
                           : (idDepth == 0) ? idBlock
                                            : _matTypeF.id;

    /* Dimensions of block */
    size_t dimBlock[2] = {0, 0};

    /* Row dimension */
    for (size_t i = 0; i < mat -> dim[1]; i++)
      {
        size_t locBlock[2] = {loc[0][0], i};
        int success = _mat_cp_struct_block_dim(mat, depth, depthTot, loc, locBlock, idBlock, 0, &dimBlock[0]);

        /* Must consider the next block */
        if (success == -1)
            continue;

        /* Found the dimension */
        else
            break;
      }

    /* Column dimension */
    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        size_t locBlock[2] = {i, loc[0][1]};
        int success = _mat_cp_struct_block_dim(mat, depth, depthTot, loc, locBlock, idBlock, 1, &dimBlock[1]);

        /* Must consider the next block */
        if (success == -1)
            continue;

        /* Found the dimension */
        else
            break;
      }

    /* New matrix */
    mat_t *matBlockCp = mat_new(idBlock, dimBlock, (strcmp(idBlock, _matTypeNull.id)) ? depthTot > depth + 1 : false);


    /* Ordinary Matrix */

    if (matBlockCp -> mblock == NULL)
        return matBlockCp;


    /* Block Matrix */

    /* Go one block deeper */
    size_t **locCp = malloc(sizeof(size_t*) * (depth + 2));

    for (size_t i = 0; i < depth + 1; i++)
      {
        locCp[i] = malloc(sizeof(size_t) * 2);

        locCp[i][0] = loc[i][0];
        locCp[i][1] = loc[i][1];
      }

    locCp[depth + 1] = malloc(sizeof(size_t) * 2);

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < dimBlock[0]; i++)
      {
        for (size_t j = 0; j < dimBlock[1]; j++)
          {
            /* Location of next block */
            locCp[depth + 1][0] = i;
            locCp[depth + 1][1] = j;

            size_t indexMem = 0;
            int success = mat_mget_index(matBlockCp, locCp[depth + 1], &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            mat_t *matBlock = mat_cp_struct_block(mat, depth + 1, locCp,
                                                    (id == NULL || idDepth == 1 || idDepth == -1) ? NULL : id,
                                                    (idDepth > 0 && id != NULL) ? idDepth - 1 : (idDepth < -1 && id != NULL) ? idDepth + 1 : idDepth);

            mat_fset_block(matBlockCp, locCp[depth + 1], matBlock);
          }
      }

    /* Free memory */
    free(visitedIndices);

    for (size_t i = 0; i < depth + 2; i++)
        free(locCp[i]);

    free(locCp);

    return matBlockCp;
}


mat_t *mat_cp_structs_block(mats_t *mats, size_t depth, size_t **loc, const char *id, int idDepth)
{
    /*

        Copy the structure of mat's block at given location

    */

    /* Check for mat NULL */
    for (size_t i = 0; i < mats -> size; i++)
      {
        if (mats_get_mat(mats, i) == NULL)
          {
            /* Exclude the i'th matrix */
            mats_t *matsCp = mats_new(0);

            for (size_t j = 0; j < mats -> size; j++)
              {
                if (i == j)
                    continue;

                mats_set_mat(matsCp, matsCp -> size, mats_get_mat(mats, j));
              }

            /* Get the matrix */
            mat_t *matCp = mat_cp_structs_block(matsCp, depth, loc, id, idDepth);

            /* Free memory */
            matsCp = mats_free(matsCp);

            return matCp;
          }
      }

    /* First matrix (arbitrary choice) */
    mat_t *mat = mats_get_mat(mats, 0);

    /* Check for the same block structure */
    for (size_t i = 1; i < mats -> size; i++)
      {
        if (!mat_comp_bstruct(mat, mats_get_mat(mats, i)))
          {
            printf("Cannot copy the similar block at depth %ld and location (%ld, %ld) of a collection of matrices if they do not have the same block structure.\n", depth, loc[0][0], loc[0][1]);
            exit(1);

            return NULL;
          }
      }

    /* Check for ordinary matrix */
    if (mat -> mblock == NULL)
      {
        printf("Cannot copy the similar block at depth %ld and location (%ld, %ld) of a collection of matrices if they are ordinary matrices.\n", depth, loc[0][0], loc[0][1]);
        exit(1);

        return NULL;
      }

    /* Check for depth */
    size_t depthTot = mat_get_depth(mat);

    if (depthTot < depth)
      {
        printf("Cannot copy the similar block at depth %ld and location (%ld, %ld) of a collection of matrices if the total depth is just %ld.\n", depth, loc[0][0], loc[0][1], depthTot);
        exit(1);

        return NULL;
      }

    /* Type of block */
    const char *idBlock = _matTypeF.id;
    _mat_cp_struct_block_id(mat, depth, depthTot, loc, &idBlock);

    for (size_t i = 1; i < mats -> size; i++)
      {
        const char *idBlock_ = _matTypeF.id;
        _mat_cp_struct_block_id(mats_get_mat(mats, i), depth, depthTot, loc, &idBlock_);

        idBlock = (!strcmp(idBlock, idBlock_)) ? idBlock : _matTypeF.id;
      }

    idBlock = (id != NULL) ? id
                           : (idDepth == 0) ? idBlock
                                            : _matTypeF.id;

    /* Dimensions of block */
    size_t dimBlock[2] = {0, 0};

    /* Row dimension */
    for (size_t i = 0; i < mat -> dim[1]; i++)
      {
        size_t locBlock[2] = {loc[0][0], i};
        int success = _mat_cp_struct_block_dim(mat, depth, depthTot, loc, locBlock, idBlock, 0, &dimBlock[0]);

        /* Must consider the next block */
        if (success == -1)
            continue;

        /* Found the dimension */
        else
            break;
      }

    /* Column dimension */
    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        size_t locBlock[2] = {i, loc[0][1]};
        int success = _mat_cp_struct_block_dim(mat, depth, depthTot, loc, locBlock, idBlock, 1, &dimBlock[1]);

        /* Must consider the next block */
        if (success == -1)
            continue;

        /* Found the dimension */
        else
            break;
      }

    /* New matrix */
    mat_t *matBlockCp = mat_new(idBlock, dimBlock, (strcmp(idBlock, _matTypeNull.id)) ? depthTot > depth + 1 : false);


    /* Ordinary Matrix */

    if (matBlockCp -> mblock == NULL)
        return matBlockCp;


    /* Block Matrix */

    /* Go one block deeper */
    size_t **locCp = malloc(sizeof(size_t*) * (depth + 2));

    for (size_t i = 0; i < depth + 1; i++)
      {
        locCp[i] = malloc(sizeof(size_t) * 2);

        locCp[i][0] = loc[i][0];
        locCp[i][1] = loc[i][1];
      }

    locCp[depth + 1] = malloc(sizeof(size_t) * 2);

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < dimBlock[0]; i++)
      {
        for (size_t j = 0; j < dimBlock[1]; j++)
          {
            /* Location of next block */
            locCp[depth + 1][0] = i;
            locCp[depth + 1][1] = j;

            size_t indexMem = 0;
            int success = mat_mget_index(matBlockCp, locCp[depth + 1], &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            mat_t *matBlock = mat_cp_structs_block(mats, depth + 1, locCp,
                                                    (id == NULL || idDepth == 1 || idDepth == -1) ? NULL : id,
                                                    (idDepth > 0 && id != NULL) ? idDepth - 1 : (idDepth < -1 && id != NULL) ? idDepth + 1 : idDepth);

            mat_fset_block(matBlockCp, locCp[depth + 1], matBlock);
          }
      }

    /* Free memory */
    free(visitedIndices);

    for (size_t i = 0; i < depth + 2; i++)
        free(locCp[i]);

    free(locCp);

    return matBlockCp;
}


mat_t *mat_cp_struct_block_rc(mat_t *mat1, mat_t *mat2, size_t depth, size_t **loc)
{
    /*

        Copy the structure of mat's block at given location by looking at its row and column

    */

    /* Check for mat NULL */
    if (mat1 == NULL || mat2 == NULL)
        return NULL;

    /* Check for ordinary mat */
    if (mat1 -> mblock == NULL || mat2 -> mblock == NULL)
      {
        printf("Cannot copy mat's block structure if mat is an ordinary matrix.\n");
        exit(1);

        return NULL;
      }

    /* Check for depth */
    size_t depthTot = mat_get_depth(mat1);

    if (depthTot < depth)
      {
        printf("Cannot copy a block in mat at depth %ld and location (%ld, %ld) if the max depth is %ld.\n", depth, loc[0][0], loc[0][1], depthTot);
        exit(1);

        return NULL;
      }

    /* Type of block */
    const char *idBlock = _matTypeF.id;

    /* Dimensions of block */
    size_t dimBlock[2] = {0, 0};

    /* Row dimension */
    for (size_t i = 0; i < mat2 -> dim[1]; i++)
      {
        size_t locBlock[2] = {loc[0][0], i};
        int success = _mat_cp_struct_block_dim(mat1, depth, depthTot, loc, locBlock, NULL, 0, &dimBlock[0]);

        /* Must consider the next block */
        if (success == -1)
            continue;

        /* Found the dimension */
        else
            break;
      }

    /* Column dimension */
    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        size_t locBlock[2] = {i, loc[0][1]};
        int success = _mat_cp_struct_block_dim(mat2, depth, depthTot, loc, locBlock, NULL, 1, &dimBlock[1]);

        /* Must consider the next block */
        if (success == -1)
            continue;

        /* Found the dimension */
        else
            break;
      }

    /* New matrix */
    mat_t *matBlockCp = mat_new(idBlock, dimBlock, depthTot > depth + 1);


    /* Ordinary Matrix */

    if (matBlockCp -> mblock == NULL)
        return matBlockCp;


    /* Block Matrix */

    /* Go one block deeper */
    size_t **locCp = malloc(sizeof(size_t*) * (depth + 2));

    for (size_t i = 0; i < depth + 1; i++)
      {
        locCp[i] = malloc(sizeof(size_t) * 2);

        locCp[i][0] = loc[i][0];
        locCp[i][1] = loc[i][1];
      }

    locCp[depth + 1] = malloc(sizeof(size_t) * 2);

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < dimBlock[0]; i++)
      {
        for (size_t j = 0; j < dimBlock[1]; j++)
          {
            /* Location of next block */
            locCp[depth + 1][0] = i;
            locCp[depth + 1][1] = j;

            size_t indexMem = 0;
            int success = mat_mget_index(matBlockCp, locCp[depth + 1], &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            mat_fset_block(matBlockCp, locCp[depth + 1], mat_cp_struct_block_rc(mat1, mat2, depth + 1, locCp));
          }
      }

    /* Free memory */
    free(visitedIndices);

    for (size_t i = 0; i < depth + 2; i++)
        free(locCp[i]);

    free(locCp);

    return matBlockCp;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_cp_struct_block_dim(mat_t *mat, size_t depth, size_t depthTot, size_t **loc, size_t *locBlock, const char *id, size_t index, size_t *dim)
{
    /*

        Get the dimensions

    */

    mat_t *matBlock = mat_fget_block(mat, locBlock);

    /* Block does not exist */
    if (matBlock == NULL || (depth != depthTot && matBlock -> mtype == &_matTypeNull))
      {
        matBlock = mat_ffree(matBlock);

        return -1;
      }

    /* Can get the dimension */
    if (depth == 0)
      {
        /* id either NULL (default) or not 'null' */
        if (id == NULL || strcmp(id, _matTypeNull.id))
            *dim = matBlock -> dim[index];

        /* If id is 'null', we need to get the full dimension */
        else
          {
            size_t fdim[2] = {0, 0};
            mat_get_fdim(matBlock, fdim);
            *dim = fdim[index];
          }
      }

    /* Must go to deeper blocks */
    else
      {
        loc++;

        for (size_t i = 0; i < (index == 0) ? matBlock -> dim[1] : matBlock -> dim[0]; i++)
          {
            /* Location of the next block */
            locBlock[0] = (index == 0) ? loc[0][0] : i;
            locBlock[1] = (index == 0) ? i : loc[0][1];

            /* Get the dimension */
            int success = _mat_cp_struct_block_dim(matBlock, depth - 1, depthTot, loc, locBlock, id, index, dim);

            /* Must consider the next block */
            if (success == -1)
                continue;

            /* Found the dimension */
            else
                break;
          }

        --loc;
      }

    matBlock = mat_ffree(matBlock);

    return 0;
}


static int _mat_cp_struct_block_id(mat_t *mat, size_t depth, size_t depthTot, size_t **loc, const char **id)
{
    /*

        Get the dimensions

    */

    mat_t *matBlock = mat_fget_block(mat, loc[0]);

    /* Block does not exist */
    if (matBlock == NULL || (depth != depthTot && matBlock -> mtype == &_matTypeNull))
      {
        matBlock = mat_ffree(matBlock);

        return -1;
      }

    int success = 0;

    /* Can get the type id */
    if (depth == 0)
      {
        *id = matBlock -> mtype -> id;
      }

    /* Must go to deeper blocks */
    else
      {
        loc++;

        /* Get the dimension */
        success = _mat_cp_struct_block_id(matBlock, depth - 1, depthTot, loc, id);

        --loc;
      }

    matBlock = mat_ffree(matBlock);

    return success;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_cp_struct_sim(mat_t *mat1, mat_t *mat2, const char *id, int idDepth)
{
    /*

        Copy the similarities of the structures of mat1 and mat2

    */

    /* Check for mat1 and mat2 NULL */
    if (mat1 == NULL)
        return mat_cp_struct(mat2, id, idDepth);

    if (mat2 == NULL)
        return mat_cp_struct(mat1, id, idDepth);

    /* Check dimensions */
    size_t fdim1[2], fdim2[2];
    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);

    if (fdim1[0] != fdim2[0] || fdim1[1] != fdim2[1])
      {
        printf("Cannot copy the structs of mat1 and mat2 with total dimensions (%ld,%ld) and (%ld,%ld).", fdim1[0], fdim1[1], fdim2[0], fdim2[1]);
        exit(1);

        return NULL;
      }

    /* Matrix id */
    const char *idMat = (id == NULL || idDepth == -1) ? (!strcmp(mat1 -> mtype -> id, mat2 -> mtype -> id)) ? mat1 -> mtype -> id : _matTypeF.id
                                                      : id;

    /* Compare block structs */
    bool comp = mat_comp_bstruct(mat1, mat2);
    bool block = comp && mat1 -> mblock != NULL;

    /* Matrix dimensions */
    size_t dim[2];
    dim[0] = (block) ? mat1 -> dim[0] : fdim1[0];
    dim[1] = (block) ? mat1 -> dim[1] : fdim1[1];

    /* New matrix */
    mat_t *matCp = mat_new((char*) idMat, dim, block);


    /* Ordinary Matrix */

    if (!block)
        return matCp;


    /* Block matrix */

    /* Block dimensions */
    for (size_t i = 0; i < matCp -> dim[0]; i++)
      {
        matCp -> mblock -> bdim[0][i] = mat1 -> mblock -> bdim[0][i];
        matCp -> mblock -> bfdim[0][i] = mat1 -> mblock -> bfdim[0][i];
      }

    for (size_t i = 0; i < matCp -> dim[1]; i++)
      {
        matCp -> mblock -> bdim[1][i] = mat1 -> mblock -> bdim[1][i];
        matCp -> mblock -> bfdim[1][i] = mat1 -> mblock -> bfdim[1][i];
      }

    matCp -> mblock -> fdim[0] = mat1 -> mblock -> fdim[0];
    matCp -> mblock -> fdim[1] = mat1 -> mblock -> fdim[1];

    /* Blocks */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < matCp -> dim[0]; i++)
      {
        for (size_t j = 0; j < matCp -> dim[1]; j++)
          {
            size_t loc[2] = {i, j};
            size_t indexMem = 0;

            int success = mat_mget_index(matCp, loc, &indexMem);

            /* Skip if a block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* Block in mat1's and mat2's matrices */
            mat_t *mat1Block = mat_fget_block(mat1, loc);
            mat_t *mat2Block = mat_fget_block(mat2, loc);

            /* Block does not exist */
            if ((mat1Block == NULL && mat2Block != NULL) || (mat2Block == NULL && mat1Block != NULL))
              {
                mat_t *matBlock = (mat1Block == NULL) ? mat2Block : mat1Block;

                ((mat_t**) matCp -> matrix)[indexMem] = matCp -> mblock -> rep(mat_cp_struct(matBlock, id, (idDepth <= -1) ? idDepth : idDepth - 1), loc);

                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);

                continue;
              }

            if (mat1Block == NULL && mat2Block == NULL)
              {
                size_t dimBlock[2] = {mat1 -> mblock -> bfdim[0][loc[0]], mat1 -> mblock -> bfdim[1][loc[1]]};
                mat1Block = mat_fnew(_matTypeNull.id, dimBlock, false);
                mat2Block = mat_fnew(_matTypeNull.id, dimBlock, false);
              }

            /* Copy the struct of mat1's and mat2's blocks */
            mat_t *matBlock = mat_cp_struct_sim(mat1Block, mat2Block,
                                                    (id == NULL || idDepth - 1 >= 0 || idDepth + 1 <= 0) ? NULL : id,
                                                    (idDepth > 0 && id != NULL) ? idDepth - 1 : (idDepth < 0 && id != NULL) ? idDepth + 1 : idDepth);

            ((mat_t**) matCp -> matrix)[indexMem] = matCp -> mblock -> rep(matBlock, loc);

            /* Free memory */
            mat1Block = mat_ffree(mat1Block);
            mat2Block = mat_ffree(mat2Block);
          }
      }

    /* Free memory */
    free(visitedIndices);

    return matCp;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------------   Setters   --------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_set_label(mat_t *mat, char *label)
{
    /*

        Set the label of a matrix (mat)

    */

    /* Unset the label */
    if (label == NULL)
      {
        free(mat -> label);
        mat -> label = NULL;

        return 0;
      }

    /* Reallocate memory */
    mat -> label = realloc(mat -> label, strlen(label) + 1);
    snprintf(mat -> label, strlen(label) + 1, "%s", label);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_set_value(mat_t *mat, size_t *loc, double value)
{
    /*

        Insert an element (value) into a matrix (mat) at location loc.

    */

    /* No matrix provided */
    if (mat == NULL)
        return 0;

    /* Block matrix */
    if (mat -> mblock != NULL)
      {
        /* Location of block in mat */
        size_t locBlock[2] = {0, 0};

        /* Location withing block */
        size_t locInBlock[2] = {loc[0], loc[1]};

        /* Get the row block location */
        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            if (mat -> mblock -> bfdim[0][i] <= locInBlock[0])
              {
                locBlock[0] += 1;
                locInBlock[0] -= mat -> mblock -> bfdim[0][i];

                continue;
              }

            break;
          }

        /* Get the column block location */
        for (size_t i = 0; i < mat -> dim[1]; i++)
          {
            if (mat -> mblock -> bfdim[1][i] <= locInBlock[1])
              {
                locBlock[1] += 1;
                locInBlock[1] -= mat -> mblock -> bfdim[1][i];

                continue;
              }

            break;
          }

        /* Set the value in the block */
        mat_t *matBlock = mat_fget_block(mat, locBlock);
        mat_set_value(matBlock, locInBlock, value);
        matBlock = mat_ffree(matBlock);

        return 0;
      }


    /* Index of location in memory */
    size_t index = 0;
    int success = mat_mget_index(mat, loc, &index);

    /* Element does not exist */
    if (success == -1)
        return 0;

    ((double*) mat -> matrix)[index] = value;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_set_matrix_o(mat_t *mat1, mat_t *mat2);
static int _mat_set_matrix_b(mat_t *mat1, mat_t *mat2);
static int _mat_set_matrix_g(mat_t *mat1, mat_t *mat2);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_set_matrix(mat_t *mat1, mat_t *mat2)
{
    /*

        Insert mat1 into mat2.

        Note: Only inserts values that are expected in mat1

    */

    /* Check for NULL */
    if (mat1 == NULL)
        return 0;

    if (mat2 == NULL)
      {
        printf("Cannot insert a matrix into a NULL pointer.\n");
        exit(1);

        return 1;
      }

    /* Check for full dimensions */
    size_t fdim1[2] = {0, 0};
    size_t fdim2[2] = {0, 0};

    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);

    if (fdim1[0] != fdim2[0] || fdim1[1] != fdim2[1])
      {
        printf("Cannot insert a matrix of dimensions (%ld, %ld) into a matrix of dimensions (%ld, %ld).\n", fdim2[0], fdim2[1], fdim1[0], fdim1[1]);
        exit(1);

        return 1;
      }

    bool comp = mat_comp_bstruct(mat1, mat2);

    /* Ordinary Matrices */
    if (comp && mat1 -> mblock == NULL)
      {
        _mat_set_matrix_o(mat1, mat2);
      }

    /* Block Matrices of the same Structure */
    else if (comp && mat1 -> mblock != NULL)
      {
        _mat_set_matrix_b(mat1, mat2);
      }

    /* Block Matrices of different structure or mixed Ordinary and Block matrix */
    else
      {
        _mat_set_matrix_g(mat1, mat2);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_set_matrix_o(mat_t *mat1, mat_t *mat2)
{
    /*

        Insert ordinary matrix mat1 into ordinary matrix mat2.

        Note: Only inserts values that are expected in mat1

    */

    /* Insert mat2 into mat1 */
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds1[2] = {0, 0};
        mat_get_cbounds(mat1, i, cBounds1);

        size_t cBounds2[2] = {0, 0};
        mat_get_cbounds(mat2, i, cBounds2);

        size_t cBounds[2] = {0, 0};
        cBounds[0] = (cBounds1[0] >= cBounds2[0]) ? cBounds1[0] : cBounds2[0];
        cBounds[1] = (cBounds1[1] <= cBounds2[1]) ? cBounds1[1] : cBounds2[1];

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            size_t loc[2] = {i, j};
            mat_set_value(mat1, loc, mat_get_value(mat2, loc));
          }
      }

    return 0;
}


static int _mat_set_matrix_b(mat_t *mat1, mat_t *mat2)
{
    /*

        Insert block matrix mat1 into block matrix mat2.

        Note: Only inserts values that are expected in mat1

    */

    /* Insert mat2 into mat1 */
    for (size_t i = 0; i < mat1 -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds1[2] = {0, 0};
        mat_get_cbounds(mat1, i, cBounds1);

        size_t cBounds2[2] = {0, 0};
        mat_get_cbounds(mat2, i, cBounds2);

        size_t cBounds[2] = {0, 0};
        cBounds[0] = (cBounds1[0] >= cBounds2[0]) ? cBounds1[0] : cBounds2[0];
        cBounds[1] = (cBounds1[1] <= cBounds2[1]) ? cBounds1[1] : cBounds2[1];

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            size_t loc[2] = {i, j};

            mat_t *matBlock1 = mat_fget_block(mat1, loc);
            mat_t *matBlock2 = mat_fget_block(mat2, loc);

            mat_set_matrix(matBlock1, matBlock2);

            matBlock1 = mat_ffree(matBlock1);
            matBlock2 = mat_ffree(matBlock2);
          }
      }

    return 0;
}


static int _mat_set_matrix_g(mat_t *mat1, mat_t *mat2)
{
    /*

        Insert general matrix mat2 into general matrix mat1.

        Note: Only inserts values that are expected in mat1

    */

    size_t fdim1[2] = {0, 0};
    size_t fdim2[2] = {0, 0};

    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);

    /* Insert mat2 into mat1 */
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < fdim1[0]; i++)
      {
        for (size_t j = 0; j < fdim1[1]; j++)
          {
            size_t loc[2] = {i, j};
            mat_set_value(mat1, loc, mat_get_value(mat2, loc));
          }
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _mat_set_block_test(mat_t *mat, mat_t *matBlock, size_t *loc);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_set_block(mat_t *mat, size_t *loc, mat_t *matBlock)
{
    /*

        Insert a block into a mat at location loc.

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        printf("Cannot set a block in an ordinary matrix. Use 'mat_set_value' to set a value instead.\n");
        exit(1);

        return 1;
      }

    /* Index of location in memory */
    size_t indexMem;
    int success = mat_mget_index(mat, loc, &indexMem);

    /* Block does not exist */
    if (success == -1)
        return 0;

    /* Location in memory */
    size_t locMem[2];
    mat_mget_loc(mat, loc, locMem);

    /* Set and test the dimensions of the blocks; exit if an error was found */
    bool test = _mat_set_block_test(mat, matBlock, locMem);
    (void) test;

    /* Must remove existing matrices */
    if (((mat_t**) mat -> matrix)[indexMem] != NULL || matBlock == NULL)
      {
        ((mat_t**) mat -> matrix)[indexMem] = mat_free(((mat_t**) mat -> matrix)[indexMem]);
      }

    /* Insert the matrix */
    ((mat_t**) mat -> matrix)[indexMem] = mat -> mblock -> rep(mat_cp(matBlock), locMem);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int mat_fset_block(mat_t *mat, size_t *loc, mat_t *matBlock)
{
    /*

        Insert a block into a mat at location loc.

        This function inserts block directly and does not create a copy.

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        printf("Cannot set a block in an ordinary matrix. Use 'mat_set_value' to set a value instead.\n");
        exit(1);

        return 1;
      }

    /* Index of location in memory */
    size_t indexMem;
    int success = mat_mget_index(mat, loc, &indexMem);

    /* Block does not exist */
    if (success == -1)
        return 0;

    /* Location in memory */
    size_t locMem[2];
    mat_mget_loc(mat, loc, locMem);

    /* Set and test the dimensions of the blocks; exit if an error was found */
    bool test = _mat_set_block_test(mat, matBlock, locMem);
    (void) test;

    /* Insert the matrix */
    ((mat_t**) mat -> matrix)[indexMem] = mat -> mblock -> rep(matBlock, locMem);

    return 0;
}


static bool _mat_set_block_test(mat_t *mat, mat_t *matBlock, size_t *loc)
{
    /*

        Test if the block matBlock can be inserted into the matrix mat at location loc

        TODO: Getting issues with fdim if lower level block matrices are not yet filled with ordinary ones...

    */

    if (mat -> mblock == NULL)
      {
        printf("Can only insert a block into a matrix with block structure.\n");
        exit(1);

        return false;
      }

    if (matBlock == NULL)
        return true;


    /* Make sure the row dimension is the same for all matrices in a row */

    if (matBlock -> mtype == &_matTypeNull && mat -> mfile != NULL)
      {
        if (mat -> mblock -> bfdim[0][loc[0]] == 0)
          {
            mat -> mblock -> bfdim[0][loc[0]] = matBlock -> dim[0];
            mat -> mblock -> fdim[0] += mat -> mblock -> bfdim[0][loc[0]];
          }

        else
          {
            if (mat -> mblock -> bfdim[0][loc[0]] != matBlock -> dim[0])
              {
                printf("Cannot insert a matrix into a block matrix at location (%ld,%ld) with total row dimension %ld if %ld is expected.\n",
                       loc[0], loc[1], (matBlock -> mblock == NULL) ? matBlock -> dim[0] : matBlock -> mblock -> fdim[0], mat -> mblock -> bfdim[0][loc[0]]);

                exit(1);

                return false;
              }
          }
      }

    else
      {
        if (mat -> mblock -> bdim[0][loc[0]] == 0)
          {
            mat -> mblock -> bdim[0][loc[0]] = matBlock -> dim[0];
          }

        else
          {
            if (mat -> mblock -> bdim[0][loc[0]] != matBlock -> dim[0])
              {
                printf("Cannot insert a matrix into a block matrix at location (%ld,%ld) with row dimension %ld if %ld is expected.\n",
                       loc[0], loc[1], matBlock -> dim[0], mat -> mblock -> bdim[0][loc[0]]);

                exit(1);

                return false;
              }
          }

        if (mat -> mblock -> bfdim[0][loc[0]] == 0)
          {
            mat -> mblock -> bfdim[0][loc[0]] = (matBlock -> mblock == NULL) ? matBlock -> dim[0] : matBlock -> mblock -> fdim[0];
            mat -> mblock -> fdim[0] += mat -> mblock -> bfdim[0][loc[0]];
          }

        else
          {
            if (mat -> mblock -> bfdim[0][loc[0]] != ((matBlock -> mblock == NULL) ? matBlock -> dim[0] : matBlock -> mblock -> fdim[0]))
              {
                printf("Cannot insert a matrix into a block matrix at location (%ld,%ld) with total row dimension %ld if %ld is expected.\n",
                       loc[0], loc[1], (matBlock -> mblock == NULL) ? matBlock -> dim[0] : matBlock -> mblock -> fdim[0], mat -> mblock -> bfdim[0][loc[0]]);

                exit(1);

                return false;
              }
          }
      }


    /* Make sure the column dimension is the same for all matrices in a column */

    if (matBlock -> mtype == &_matTypeNull && mat -> mfile != NULL)
      {
        if (mat -> mblock -> bfdim[1][loc[1]] == 0)
          {
            mat -> mblock -> bfdim[1][loc[1]] = matBlock -> dim[1];
            mat -> mblock -> fdim[1] += mat -> mblock -> bfdim[1][loc[1]];
          }

        else
          {
            if (mat -> mblock -> bfdim[1][loc[1]] != matBlock -> dim[1])
              {
                printf("Cannot insert a matrix into a block matrix at location (%ld,%ld) with total row dimension %ld if %ld is expected.\n",
                       loc[0], loc[1], (matBlock -> mblock == NULL) ? matBlock -> dim[1] : matBlock -> mblock -> fdim[1], mat -> mblock -> bfdim[1][loc[1]]);

                exit(1);

                return false;
              }
          }
      }

    else
      {
        if (mat -> mblock -> bdim[1][loc[1]] == 0)
          {
            mat -> mblock -> bdim[1][loc[1]] = matBlock -> dim[1];
          }

        else
          {
            if (mat -> mblock -> bdim[1][loc[1]] != matBlock -> dim[1])
              {
                printf("Cannot insert a matrix into a block matrix at location (%ld,%ld) with row dimension %ld if %ld is expected.\n",
                       loc[0], loc[1], matBlock -> dim[1], mat -> mblock -> bdim[1][loc[1]]);

                exit(1);

                return false;
              }
          }

        if (mat -> mblock -> bfdim[1][loc[1]] == 0)
          {
            mat -> mblock -> bfdim[1][loc[1]] = (matBlock -> mblock == NULL) ? matBlock -> dim[1] : matBlock -> mblock -> fdim[1];
            mat -> mblock -> fdim[1] += mat -> mblock -> bfdim[1][loc[1]];
          }

        else
          {
            if (mat -> mblock -> bfdim[1][loc[1]] != ((matBlock -> mblock == NULL) ? matBlock -> dim[1] : matBlock -> mblock -> fdim[1]))
              {
                printf("Cannot insert a matrix into a block matrix at location (%ld,%ld) with total row dimension %ld if %ld is expected.\n",
                       loc[0], loc[1], (matBlock -> mblock == NULL) ? matBlock -> dim[1] : matBlock -> mblock -> fdim[1], mat -> mblock -> bfdim[1][loc[1]]);

                exit(1);

                return false;
              }
          }
      }

    return true;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_set_array(mat_t *mat, double *array)
{
    /*

        Set the matrix array in mat.

        Note: Assumes matrix to be a 1-d array of size mat -> dim[0] * mat -> dim[1] (ie NOT mat -> dim[2] for e.g. diagonal matrix)

    */

    /* Block Matrix */
    if (mat -> mblock != NULL)
      {
        printf("Cannot set the matrix (double array) in a block matrix. Use 'mat_set_block' to set a block instead.\n");
        exit(1);

        return 1;
      }

    size_t index = 0;

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            size_t loc[2] = {i, j};

            mat_set_value(mat, loc, array[index++]);
          }
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------------   Getters   --------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *mat_get_label(mat_t *mat)
{
    /*

        Get the matrix label

    */

    return mat -> label;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_dim(mat_t *mat, size_t *dimMat)
{
    /*

        Get the dimensions in mat's matrix

    */

    dimMat[0] = mat -> dim[0];
    dimMat[1] = mat -> dim[1];

    return 0;
}


int mat_mget_dim(mat_t *mat, size_t *dimMem)
{
    /*

        Get the dimensions in mat's memory

    */

    dimMem[0] = mat -> dim[0];
    dimMem[1] = mat -> dim[1];

    for (size_t i = 0; i < mat -> mltrfs -> size; i++)
      {
        mat -> mltrfs -> mltrf[i] -> dim(mat, dimMem);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_bdim(mat_t *mat, size_t depth, size_t **loc, size_t *bdimMat)
{
    /*

        Get the dimensions of a block in mat's matrix at given depth and location

    */

    if (mat -> mblock == NULL)
      {
        printf("Cannot get the block dimensions of an ordinary matrix.\n");
        exit(1);

        return 1;
      }

    mat_t *matBlock = mat_fget_block(mat, loc[0]);

    /* Block does not exist */
    if (matBlock == NULL)
        return -1;

    /* Get the dimensions */
    if (depth == 0)
      {
        bdimMat[0] = matBlock -> dim[0];
        bdimMat[1] = matBlock -> dim[1];
      }

    /* Need to go deeper */
    else
      {
        mat_get_bdim(matBlock, depth - 1, ++loc, bdimMat);
        loc--; // Reset the pointer location in memory
      }

    matBlock = mat_ffree(matBlock);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_bfdim(mat_t *mat, size_t *loc, size_t *bfdimMat)
{
    /*

        Get the full dimensions of a block in mat's matrix

    */

    if (mat -> mblock == NULL)
      {
        printf("Cannot get the full dimensions of all blocks in an ordinary matrix.\n");
        exit(1);

        return 1;
      }

    bfdimMat[0] = mat -> mblock -> bfdim[0][loc[0]];
    bfdimMat[1] = mat -> mblock -> bfdim[1][loc[1]];

    return 0;
}


int mat_mget_bfdim(mat_t *mat, size_t *loc, size_t *bfdimMem)
{
    /*

        Get the full dimensions of a block in mat's memory

    */

    /* Dimension in matrix */
    mat_get_bfdim(mat, loc, bfdimMem);

    /* Dimension in memory */
    for (size_t i = 0; i < mat -> mltrfs -> size; i++)
      {
        mat -> mltrfs -> mltrf[i] -> dim(mat, bfdimMem);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_fdim(mat_t *mat, size_t *fdimMat)
{
    /*

        Get the full dimensions in mat's matrix

    */

    fdimMat[0] = (mat -> mblock == NULL) ? mat -> dim[0] : mat -> mblock -> fdim[0];
    fdimMat[1] = (mat -> mblock == NULL) ? mat -> dim[1] : mat -> mblock -> fdim[1];

    return 0;
}


int mat_mget_fdim(mat_t *mat, size_t *fdimMem)
{
    /*

        Get the full dimensions in mat's memory

    */

    /* Dimension in matrix */
    mat_get_fdim(mat, fdimMem);

    /* Dimension in memory */
    for (size_t i = 0; i < mat -> mltrfs -> size; i++)
      {
        mat -> mltrfs -> mltrf[i] -> dim(mat, fdimMem);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_rbounds(mat_t *mat, size_t col, size_t *rbounds)
{
    /*

        Get the row bounds of mat for a given column

    */

    return mat -> mtype -> rbounds(mat -> dim, col, rbounds);
}


int mat_get_cbounds(mat_t *mat, size_t row, size_t *cbounds)
{
    /*

        Get the column bounds of mat for a given row

    */

    return mat -> mtype -> cbounds(mat -> dim, row, cbounds);
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_loc(mat_t *mat, size_t indexMat, size_t *locMat)
{
    /*

        Get the location in mat's matrix, given the index in mat's matrix

    */

    mat -> mtype -> loc(mat -> dim, indexMat, locMat);

    return 0;
}


int mat_mget_loc(mat_t *mat, size_t *locMat, size_t *locMem)
{
    /*

        Get the location in mat's memory, given the location in mat's matrix

    */

    locMem[0] = locMat[0];
    locMem[1] = locMat[1];

    for (size_t i = 0; i < mat -> mltrfs -> size; i++)
      {
        mat -> mltrfs -> mltrf[i] -> loc(mat, locMem);
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int mat_get_index(mat_t *mat, size_t *locMat, size_t *indexMat)
{
    /*

        Get the index in mat's matrix, given the location in mat's matrix

    */

    return mat -> mtype -> index(mat -> dim, locMat, indexMat);
}


int mat_mget_index(mat_t *mat, size_t *locMat, size_t *indexMem)
{
    /*

        Get the index in mat's memory, given the location in mat's matrix

    */

    /* Location in memory */
    size_t locMem[2];
    mat_mget_loc(mat, locMat, locMem);

    /* Dimension in memory */
    size_t dimMem[2];
    mat_mget_dim(mat, dimMem);

    /* Index */
    return mat -> mtype -> index(dimMem, locMem, indexMem);
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double mat_get_value(mat_t *mat, size_t *loc)
{
    /*

        Get the value of mat at location loc.

    */

    return mat_get_ivalue(mat, loc, NULL);
}


double mat_get_ivalue(mat_t *mat, size_t *loc, size_t *index)
{
    /*

        Get the value and its index of mat at location loc.

    */

    /* No matrix provided */
    if (mat == NULL)
        return 0.;

    /* Block matrix */
    if (mat -> mblock != NULL)
      {
        /* Location of block in mat */
        size_t locBlock[2] = {0, 0};

        /* Location withing block */
        size_t locInBlock[2] = {loc[0], loc[1]};

        /* Get the row block location */
        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            if (mat -> mblock -> bfdim[0][i] <= locInBlock[0])
              {
                locBlock[0] += 1;
                locInBlock[0] -= mat -> mblock -> bfdim[0][i];

                continue;
              }

            break;
          }

        /* Get the column block location */
        for (size_t i = 0; i < mat -> dim[1]; i++)
          {
            if (mat -> mblock -> bfdim[1][i] <= locInBlock[1])
              {
                locBlock[1] += 1;
                locInBlock[1] -= mat -> mblock -> bfdim[1][i];

                continue;
              }

            break;
          }

        /* Get the value in the block */
        mat_t *matBlock = mat_fget_block(mat, locBlock);
        double val = mat_get_value(matBlock, locInBlock);
        matBlock = mat_ffree(matBlock);

        return val;
      }

    /* Index of location in memory */
    size_t indexMem;
    int success = mat_mget_index(mat, loc, &indexMem);

    /* Element does not exist */
    if (success == -1)
        return 0.;

    /* Store index index */
    if (index != NULL)
        *index = indexMem;

    return ((double*) mat -> matrix)[indexMem];
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_get_block(mat_t *mat, size_t *loc)
{
    /*

        Get a block in mat's matrix.

    */

    return mat_get_iblock(mat, loc, NULL);
}


mat_t *mat_get_iblock(mat_t *mat, size_t *loc, size_t *index)
{
    /*

        Get a block and its index in mat's matrix.

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        printf("Cannot get a block of an ordinary matrix. Use 'mat_get_value' to get a value instead.\n");
        exit(1);

        return NULL;
      }

    /* Index of location in memory */
    size_t indexMem;
    int success = mat_mget_index(mat, loc, &indexMem);

    /* Block does not exist */
    if (success == -1)
        return NULL;

    /* Location in memory */
    size_t locMem[2];
    mat_mget_loc(mat, loc, locMem);

    /* Get the block */
    mat_t *matBlock = mat -> mblock -> rep(mat_cp(((mat_t**) mat -> matrix)[indexMem]), locMem);

    /* Store the index */
    if (index != NULL)
        *index = indexMem;

    return matBlock;
}


/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_fget_block(mat_t *mat, size_t *loc)
{
    /*

        Get the block of mat at location loc.

        Note: This function returns a shallow copy of the block's representation in the parent's matrix. It does not return the raw block in memory.

    */

    return mat_fget_iblock(mat, loc, NULL);
}


mat_t *mat_fget_iblock(mat_t *mat, size_t *loc, size_t *index)
{
    /*

        Get the block and its index of mat at location loc.

        Note: This function returns a shallow copy of the block's representation in the parent's matrix. It does not return the raw block in memory.

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        printf("Cannot get a block of an ordinary matrix. Use 'mat_get_value' to get a value instead.\n");
        exit(1);

        return NULL;
      }

    /* Index of location in memory */
    size_t indexMem;
    int success = mat_mget_index(mat, loc, &indexMem);

    /* Block does not exist */
    if (success == -1)
      {
        size_t bfDim[2];
        mat_get_bfdim(mat, loc, bfDim);

        mat_t *matNull = mat_fnew(_matTypeNull.id, bfDim, false);

        return matNull;
      }

    /* Location in memory */
    size_t locMem[2];
    mat_mget_loc(mat, loc, locMem);

    /* Get the block */
    mat_t *matBlock = mat -> mblock -> rep(mat_fcp(((mat_t**) mat -> matrix)[indexMem]), locMem);

    /* Store the index */
    if (index != NULL)
        *index = indexMem;

    return matBlock;
}


/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_mget_block(mat_t *mat, size_t *loc)
{
    /*

        Get a block in mat's memory

    */

    return mat_mget_iblock(mat, loc, NULL);
}


mat_t *mat_mget_iblock(mat_t *mat, size_t *loc, size_t *index)
{
    /*

        Get a block and its index in mat's memory

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        printf("Cannot get a block of an ordinary matrix. Use 'mat_get_value' to get a value instead.\n");
        exit(1);

        return NULL;
      }

    /* Index of location in memory */
    size_t indexMem;
    int success = mat_mget_index(mat, loc, &indexMem);

    /* Block does not exist */
    if (success == -1)
        return NULL;

    /* Get the block */
    mat_t *matBlock = ((mat_t**) mat -> matrix)[indexMem];

    /* Store the index */
    if (index != NULL)
        *index = indexMem;

    return matBlock;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t mat_get_depth(mat_t *mat)
{
    /*

        Get the block structure's depth

    */

    size_t depth = 0;

    /* Ordinary matrix or NULL pointer */
    if (mat == NULL || mat -> mblock == NULL)
        return depth;

    for (size_t i = 0; i < mat -> dim[2]; i++)
      {
        /* Depth of block */
        size_t depthBlock = mat_get_depth(((mat_t**) mat -> matrix)[i]);

        /* Adjust depth */
        depth = (depth < depthBlock + 1) ? depthBlock + 1 : depth;
      }

    return depth;
}


//mat_t mat_fget_depth_block(mat_t *mat, size_t depth, size_t **loc)
//{
//    /*
//
//        Get the mat's block at location loc
//
//    */
//
//    if (mat == NULL)
//        return NULL;
//
//    if (mat -> mblock == NULL)
//        return NULL;
//
//    mat_t *matBlock = mat_fget_block(mat, loc[0]);
//
//    for (size_t i = 0; i < depth; i++)
//      {
//        mat_t *matBlock_ = mat_fget_block(matBlock, loc[i + 1]);
//        mat_swap(matBlock, matBlock_);
//        matBlock = mat_ffree(matBlock_);
//      }
//
//    return matBlock;
//}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Swap Matrices   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_swap(mat_t *mat1, mat_t *mat2)
{
    /*

        Swap two matrices

    */

    mat_t matTemp = *mat1;

    *mat1 = *mat2;
    *mat2 = matTemp;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   p-Norm of a Matrix   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static double _mat_pnorm_o(mat_t *mat, size_t p);
static double _mat_pnorm_b(mat_t *mat, size_t p);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


double mat_pnorm(mat_t *mat, size_t p)
{
    /*

        Calculate the p-norm of a matrix (mat):

            |A| = |a_11|^p + ... + |a_nm|^p

    */

    /* Matrix is NULL */
    if (mat == NULL)
        return 0.;

    double norm = 0.;

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        norm = _mat_pnorm_o(mat, p);
      }

    /* Block Matrix */
    else
      {
        norm = _mat_pnorm_b(mat, p);
      }

    return norm;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static double _mat_pnorm_o(mat_t *mat, size_t p)
{
    /*

        Calculate the p-norm of an ordinary matrix (mat):

            |A| = (|a_11|^p + ... + |a_nm|^p)^(1/p)

    */

    double norm = 0;

    #pragma omp parallel for schedule(dynamic) reduction(+:norm)
    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds[2] = {0, 0};
        mat_get_cbounds(mat, i, cBounds);

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            /* Value at location */
            size_t loc[2] = {i, j};
            double val = mat_get_value(mat, loc);

            /* Skip if an element is zero */
            if (val == 0.)
                continue;

            /* Calculate a_ij^p */
            double prod = 1.;

            for (size_t k = 0; k < p; k++)
                prod *= val;

            /* Add the absolute value of a_ij^p */
            norm += fabs(prod);
          }
      }

    return pow(norm, 1. / ((double) p));
}


static double _mat_pnorm_b(mat_t *mat, size_t p)
{
    /*

        Calculate the p-norm of a block matrix (mat):

            |A| = (|a_11|^p + ... + |a_nm|^p)^(1/p)

    */

    double norm = 0;

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds[2] = {0, 0};
        mat_get_cbounds(mat, i, cBounds);

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            /* Block at location */
            size_t loc[2] = {i, j};
            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Skip if a block is NULL (ie does not exist) */
            if (matBlock == NULL)
                continue;

            /* Add the norm of a block */
            double normBlock = mat_pnorm(matBlock, p);
            norm += normBlock;

            /* Free memory */
            matBlock = mat_ffree(matBlock);
          }
      }

    return pow(norm, 1. / ((double) p));
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Reduce Matrix   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_reduce_o(mat_t *mat, mat_reduce_t *reduce);
static mat_t *_mat_reduce_b(mat_t *mat, mat_reduce_t *reduce);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_reduce(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Return the matrix that takes up the least amount of memory

    */

    /* No reduction requested */
    if (reduce == NULL || mat == NULL)
      {
        return mat;
      }

    /* Ordinary matrix */
    if (mat -> mblock == NULL)
      {
        return _mat_reduce_o(mat, reduce);
      }

    /* Block matrix */
    if (mat -> mblock != NULL)
      {
        return _mat_reduce_b(mat, reduce);
      }

    return NULL;
}


static mat_t *_mat_reduce_o(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Reduce an ordinary matrix

    */

    /* Absolute error allowed for two elements to differ */
    double err = reduce -> err * mat_pnorm(mat, 1) / ((double) (mat -> dim[0] * mat -> dim[1]));

    /* Check which types can support mat */
    for (size_t n = 0; n < _matTypesNum; n++)
      {
        /* Do not need to check for the same type */
        if (&_matTypes[n] == mat -> mtype)
            break;

        /* Do not need to check for full type */
        if (&_matTypes[n] == &_matTypeF)
            continue;

        /* Next type flag */
        bool nextType = false;

        /* Check for the type */
        #pragma omp parallel for
        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            /* Skip to the next type */
            if (nextType)
                continue;

            for (size_t j = 0; j < mat -> dim[1]; j++)
              {
                /* Skip to the next type */
                if (nextType)
                    continue;

                size_t loc[2] = {i, j};

                /* Matrix does not have correct type */
                if (!_matTypes[n].comp(mat, loc, err))
                  {
                    nextType = true;
                  }
              }
          }

        /* Found the correct type */
        if (!nextType)
          {
            mat = mat_trafo_ip(mat, (char*) _matTypes[n].id, mat -> dim, false, NULL);

            break;
          }
      }

    return mat;
}


static mat_t *_mat_reduce_b(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Reduce a block matrix

    */

    /* First, reduce blocks of mat */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Index of location in memory */
            size_t loc[2] = {i, j};
            size_t index = 0;
            int success = mat_mget_index(mat, loc, &index);

            /* Block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* Block at location in memory */
            mat_t *matBlock = mat_mget_block(mat, loc);

            /* Skip if block is NULL (ie does not exist) */
            if (matBlock == NULL)
                continue;

            /* Reduce the block */
            matBlock = mat_reduce(matBlock, reduce);
          }
      }

    /* Free memory */
    free(visitedIndices);

    /* Absolute error allowed for two elements to differ */
    double err = reduce -> err * mat_pnorm(mat, 1) / ((double) (mat -> mblock -> fdim[0] * mat -> mblock -> fdim[1]));

    /* Check which types can support mat */
    for (size_t n = 0; n < _matTypesNum; n++)
      {
        /* Do not need to check for the same type */
        if (&_matTypes[n] == mat -> mtype)
            break;

        /* Do not need to check for full type */
        if (&_matTypes[n] == &_matTypeF)
            continue;

        /* Next type flag */
        bool nextType = false;

        /* Check for the type */
        #pragma omp parallel for
        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            /* Skip to the next type */
            if (nextType)
                continue;

            for (size_t j = 0; j < mat -> dim[1]; j++)
              {
                /* Skip to the next type */
                if (nextType)
                    continue;

                /* Location in mat */
                size_t loc[2] = {i, j};

                /* Matrix does not have correct type */
                if (!_matTypes[n].comp(mat, loc, err))
                  {
                    nextType = true;
                  }
              }
          }

        /* Found the correct type */
        if (!nextType)
          {
            /* Type is 'null' -> cannot have block structure */
            if (!strcmp(_matTypeNull.id, _matTypes[n].id))
                mat = mat_trafo_ip(mat, (char*) _matTypes[n].id, mat -> mblock -> fdim, false, NULL);

            /* Other type can have block structure */
            else
                mat = mat_trafo_ip(mat, (char*) _matTypes[n].id, mat -> dim, true, NULL);

            break;
          }
      }

    return mat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Count Matrices in Struct   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_count_rec(mat_t *mat, size_t *count);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t mat_count(mat_t *mat)
{
    /*

        Count the number of matrices in mat (only count mat itself if mat is not a collection)

    */

    size_t count = 0;

    _mat_count_rec(mat, &count);

    return count;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_count_rec(mat_t *mat, size_t *count)
{
    /*

        Count the number of matrices in mat (also count mat itself)

    */

    /* Ordinary matrix */
    if (mat -> mblock == NULL)
      {
        *count += 1;
      }

    else
      {
        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            _mat_count_rec(((mat_t**) mat -> matrix)[i], count);
          }

        *count += 1;
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Diagonal Matrix   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_diag(mat_t *mat)
{
    /*

        Copy only the diagonal part of mat

    */

    /* Copy the struct */
    mat_t *matDiag = mat_cp_struct(mat, _matTypeD.id, 0);

    /* Set the diagonal */
    mat_set_matrix(matDiag, mat);

    return matDiag;
}


mat_t *mat_diag_ip(mat_t *mat)
{
    /*

        Copy only the diagonal part of mat

    */

    /* Get the diagonal matrix */
    mat_t *matDiag = mat_cp_struct(mat, _matTypeD.id, 0);

    /* Swap the matrices */
    mat_swap(mat, matDiag);

    /* Free the previous matrix */
    matDiag = mat_free(matDiag);

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_diag_func(mat_t *mat, double (*func)(double))
{
    /*

        Transform the elements of the diagonal matrix with func

    */

    /* Must have diagonal matrix */
    if (mat -> mtype != &_matTypeD)
      {
        printf("Cannot tranform the elements of a non-diagonal matrix with a scalar function.\n");
        exit(1);

        return NULL;
      }

    /* Full dimensions */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    /* Transform the elements */
    for (size_t i = 0; i < fdim[0]; i++)
      {
        size_t loc[2] = {i, i};

        mat_set_value(mat, loc, func(mat_get_value(mat, loc)));
      }

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


double mat_diag_func_sqrt(double val)
{
    /*

        Square root function for mat_diag_func

    */

    return sqrt(val);
}


double mat_diag_func_inv(double val)
{
    /*

        Inverse function for mat_diag_func

    */

    return 1. / val;
}


double mat_diag_func_invsqrt(double val)
{
    /*

        Inverse square root function for mat_diag_func (used in mat_corr)

    */

    return 1. / sqrt(val);
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Transform a Matrix Struct   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_trafo_ip(mat_t *mat, const char *type, size_t *dim, bool block, mat_reduce_t *reduce)
{
    /*

        Transform a matrix (mat) into a given type with dimensions dim and perhaps in block structure.
        If type is NULL, mat's types will be copied.

        This function overwrites mat with the transformed matrix

    */

    /* Get the trafo */
    mat_t *matTrf = mat_trafo(mat, (type != NULL) ? type : mat -> mtype -> id, dim, block, reduce);

    /* Swap matrices */
    mat_swap(mat, matTrf);

    /* Free memory */
    matTrf = mat_free(matTrf);

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_trafo_rec(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce);

static int _mat_trafo_oo(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce);
static int _mat_trafo_ob(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce);
static int _mat_trafo_bo(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce);
static int _mat_trafo_bb(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_trafo(mat_t *mat, const char *type, size_t *dim, bool block, mat_reduce_t *reduce)
{
    /*

        Transform a matrix (mat) into a given type with dimensions dim and perhaps in block structure.
        If type is NULL, mat's types will be copied.

    */

    /* New matrix */
    mat_t *matTrf = mat_new((type != NULL) ? type : mat -> mtype -> id, dim, block);

    /* Insert mat at location (0,0) into matTrf at location (0,0) */
    size_t loc[2] = {0, 0};
    size_t locTrf[2] = {0, 0};

    /* Transform the matrix */
    _mat_trafo_rec(mat, loc, matTrf, locTrf, dim, reduce);

    return matTrf;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_trafo_rec(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce)
{
    /*

        Insert a dim[0] x dim[1] part of a matrix (mat) starting at location loc into another matrix (matTrf) starting at location locTrf.

        This function is always called when a recursion is required and redirects to the oo, ob, bo, bb functions.

    */


    /* Make sure the matrices are legitimate */

    if (matTrf -> mblock == NULL)
      {
        size_t fdim[2];
        mat_get_fdim(mat, fdim);

        /* Total dimensions must match given dimensions */
        if (fdim[0] != dim[0] || fdim[1] != dim[1])
          {
            printf("Cannot transform a matrix with total dimensions (%ld,%ld) into an ordinary matrix with dimensions (%ld,%ld).\n", fdim[0], fdim[1], dim[0], dim[1]);
            exit(1);

            return 1;
          }
      }

    else
      {
        /*  Ordinary matrices must have dimensions divisible by dim */
        if (mat -> mblock == NULL && (mat -> dim[0] % dim[0] != 0 || mat -> dim[1] % dim[1] != 0))
          {
            printf("Cannot trasform an ordinary matrix with dimensions (%ld,%ld) into a block matrix with dimensions (%ld,%ld).\n", mat -> dim[0], mat -> dim[1], dim[0], dim[1]);
            printf("The dimensions of the ordinary matrix must be divisible by the dimensions of the block matrix.\n");

            exit(1);

            return 1;
          }

        /* Block matrices must have same dimensions as dim */
        if (mat -> mblock != NULL && (mat -> dim[0] != dim[0] || mat -> dim[1] != dim[1]))
          {
            printf("Cannot trasform a block matrix with dimensions (%ld,%ld) into a block matrix with dimensions (%ld,%ld).\n", mat -> dim[0], mat -> dim[1], dim[0], dim[1]);
            printf("The dimensions must be the same only the type may differ.\n");

            exit(1);

            return 1;
          }
      }


    /* Trafo functions */

    if (mat -> mblock == NULL && matTrf -> mblock == NULL)
      {
        _mat_trafo_oo(mat, loc, matTrf, locTrf, dim, reduce);
      }

    else if (mat -> mblock == NULL && matTrf -> mblock != NULL)
      {
        _mat_trafo_ob(mat, loc, matTrf, locTrf, dim, reduce);
      }

    else if (mat -> mblock != NULL && matTrf -> mblock == NULL)
      {
        _mat_trafo_bo(mat, loc, matTrf, locTrf, dim, reduce);
      }

    else
      {
        _mat_trafo_bb(mat, loc, matTrf, locTrf, dim, reduce);
      }


    /* Set the label */

    mat_set_label(matTrf, mat -> label);


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_trafo_oo(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce)
{
    /*

        Insert a dim[0] x dim[1] part of a matrix (mat) starting at location loc into another matrix (matTrf) starting at location locTrf.

        Note that both matrices must be ordinary (o-o) and that the locations are within bounds in regard to the dimensions of single matrices
        and the requested dim[0] x dim[1] submatrix.

    */

    (void) reduce;

    /* Initial locations */
    size_t loc1 = loc[0];
    size_t loc2 = loc[1];

    size_t locTrf1 = locTrf[0];
    size_t locTrf2 = locTrf[1];

    /* Insert part of mat into matTrf */
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < dim[0]; i++)
      {
        for (size_t j = 0; j < dim[1]; j++)
          {
            /* Index of location in mat's memory */
            size_t loc_[2] =  {loc1 + i, loc2 + j};
            size_t index = 0;
            int success = mat_mget_index(mat, loc_, &index);

            /* Index of location in matTrf's memory */
            size_t locTrf_[2] = {locTrf1 + i, locTrf2 + j};
            size_t indexTrf = 0;
            int successTrf = mat_mget_index(matTrf, locTrf_, &indexTrf);

            /* Element does not exist */
            if (success == -1 || successTrf == -1)
                continue;

            /* Insert the element */
            ((double*) matTrf -> matrix)[indexTrf] = ((double*) mat -> matrix)[index];
          }
      }

    /* Reset locations */
//    loc[0] = loc1;
//    loc[1] = loc2;
//
//    locTrf[0] = locTrf1;
//    locTrf[1] = locTrf2;

    return 0;
}


static int _mat_trafo_ob(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce)
{
    /*

        Insert a dim[0] x dim[1] part of a matrix (mat) starting at location loc into a block of a block matrix (matTrf) at location locTrf.

        Note that mat must be ordinary and matTrf must be block (o-b), that the locations are within bounds in regard to the dimensions of
        single matrices and the requested dim[0] x dim[1] submatrix as well as that dim is compatible with the block dimensions of matTrf.

    */

    (void) dim;

    /* Dimension of block */
    size_t dimTrfBlock[2] = {mat -> dim[0] / matTrf -> dim[0], mat -> dim[1] / matTrf -> dim[1]};


    /* Insert mat into blocks of matTrf */

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    /* Initial locations */
    size_t loc1 = loc[0];
    size_t loc2 = loc[1];

    size_t locTrf1 = locTrf[0];
    size_t locTrf2 = locTrf[1];

    for (size_t i = 0; i < matTrf -> dim[0]; i++)
      {
        for (size_t j = 0; j < matTrf -> dim[1]; j++)
          {
            /* Location of next block in mat */
            if (j > 0)
                loc[1] += dimTrfBlock[1];

            /* Index of location in matTrf's memory */
            locTrf[0] = i;
            locTrf[1] = j;
            size_t indexTrf;
            int successTrf = mat_mget_index(matTrf, locTrf, &indexTrf);

            /* Element does not exist */
            if (successTrf == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexTrf, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* New block (must start with full; can reduce later) */
            mat_t *matTrfBlock = mat_new(_matTypeF.id, dimTrfBlock, false);

            /* Insert part of mat into matTrfBlock */
            size_t locTrfBlock[2] = {0, 0};
            _mat_trafo_oo(mat, loc, matTrfBlock, locTrfBlock, dimTrfBlock, reduce);

            /* Try to reduce the block */
            matTrfBlock = mat_reduce(matTrfBlock, reduce);

            /* Insert the reduced block */
            mat_fset_block(matTrf, locTrf, matTrfBlock);
          }

        /* Location of next block in mat */
        loc[0] += dimTrfBlock[0];
        loc[1] = loc2;
      }

    /* Reset locations */
    loc[0] = loc1;
    loc[1] = loc2;

    locTrf[0] = locTrf1;
    locTrf[1] = locTrf2;

    /* Free memory */
    free(visitedIndices);

    return 0;
}


static int _mat_trafo_bo(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce)
{
    /*

        Insert a block of a block matrix (mat) location loc into an ordinary matrix (matTrf) starting at location locTrf.

        Note that mat must be block and matTrf must be ordinary (b-o), that the locations are within bounds in regard to the dimensions of
        single matrices.

    */

    (void) dim;

    /* Insert blocks of mat into matTrf */

    /* Initial locations */
    size_t loc1 = loc[0];
    size_t loc2 = loc[1];

    size_t locTrf1 = locTrf[0];
    size_t locTrf2 = locTrf[1];

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Location of next block in mat */
            if (j > 0)
                locTrf[1] += mat -> mblock -> bfdim[1][j - 1];

            /* Index of location in mat's memory */
            loc[0] = i;
            loc[1] = j;
            size_t index = 0;
            int success = mat_mget_index(mat, loc, &index);

            /* Element does not exist */
            if (success == -1)
                continue;

            /* Block to be inserted into part of matTrf */
            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Insert part of mat into matBlock (cannot simply call oo since matBlock might itself be another block matrix) */
            size_t locBlock[2] = {0, 0};
            _mat_trafo_rec(matBlock, locBlock, matTrf, locTrf, (matBlock -> mblock == NULL) ? matBlock -> dim : matBlock -> mblock -> fdim, reduce);

            /* Free memory */
            matBlock = mat_ffree(matBlock);
          }

        /* Location of next block in mat */
        locTrf[0] += mat -> mblock -> bfdim[0][i];
        locTrf[1] = locTrf2;
      }

    /* Reset locations */
    loc[0] = loc1;
    loc[1] = loc2;

    locTrf[0] = locTrf1;
    locTrf[1] = locTrf2;

    return 0;
}


static int _mat_trafo_bb(mat_t *mat, size_t *loc, mat_t *matTrf, size_t *locTrf, size_t *dim, mat_reduce_t *reduce)
{
    /*

        Insert a block of a block matrix (mat) at location loc into a block of a block matrix (matTrf) at location locTrf.

        Note that both mat and matTrf must be block (b-b), that the locations are within bounds in regard to the dimensions of
        single matrices.

    */

    (void) dim;

    /* Insert blocks of mat into blocks of matTrf */

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    /* Initial locations */
    size_t loc1 = loc[0];
    size_t loc2 = loc[1];

    size_t locTrf1 = locTrf[0];
    size_t locTrf2 = locTrf[1];

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Index of location in mat's memory */
            loc[0] = i;
            loc[1] = j;
            size_t index = 0;
            int success = mat_mget_index(mat, loc, &index);

            /* Index of location in matTrf's memory */
            locTrf[0] = i;
            locTrf[1] = j;
            size_t indexTrf = 0;
            int successTrf = mat_mget_index(matTrf, locTrf, &indexTrf);

            /* Block in matTrf does not exist */
            if (successTrf == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexTrf, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* Insert an empty matrix if a block in mat does not exist but should exist in matTrf */
            if (success == -1)
              {
                /* Empty block */
                size_t dim[2] = {mat -> mblock -> bfdim[0][i], mat -> mblock -> bfdim[0][j]};
                mat_t *matBlock = mat_new(_matTypeNull.id, dim, false);

                /* Insert the empty block */
                mat_fset_block(matTrf, locTrf, matBlock);

                continue;
              }

            /* Block to be inserted into block of matTrf */
            mat_t *matBlock = mat_get_block(mat, loc);

            /* Try to reduce the block */
            matBlock = mat_reduce(matBlock, reduce);

            /* Insert the reduced block */
            mat_fset_block(matTrf, locTrf, matBlock);
          }
      }

    /* Reset locations */
    loc[0] = loc1;
    loc[1] = loc2;

    locTrf[0] = locTrf1;
    locTrf[1] = locTrf2;

    /* Free memory */
    free(visitedIndices);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Print Matrix   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mat_print(mat_t *mat, int precision)
{
    /*

        Print a matrix

    */

    /* If mat is NULL print 'NULL' */
    if (mat == NULL)
      {
        printf("'NULL'\n\n");

        return 0;
      }

    /* Get the full matrix */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    /* Print the matrix */
    for (size_t i = 0; i < fdim[0]; i++)
      {
        for (size_t j = 0; j < fdim[1]; j++)
          {
            size_t loc[2] = {i, j};

            printf("%.*e ", precision, mat_get_value(mat, loc));
          }

        printf("\n");
      }

    printf("\n");

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LINEAR ALGEBRA     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Transpose Matrix   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_tp_o(mat_t *mat, mat_t *matTp);
static int _mat_tp_b(mat_t *mat, mat_t *matTp);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_tp(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Transpose a matrix; transposes matrix by creating new matrix

    */

    /* Check for NULL */
    if (mat == NULL)
        return NULL;

    /* Dimensions of transposed matrix */
    size_t dim[2] = {mat -> dim[1], mat -> dim[0]};

    /* Transposed matrix */
    mat_t *matTp = mat_new((char*) (mat_type_get_struct_tp(mat -> mtype -> id) -> id), dim, mat -> mblock != NULL);

    /* Transpose an Ordinary Matrix */
    if (matTp -> mblock == NULL)
      {
        _mat_tp_o(mat, matTp);
      }

    /* Transpose a Block Matrix */
    else
      {
        _mat_tp_b(mat, matTp);
      }

    /* Try to reduce the matrix */
    matTp = mat_reduce(matTp, reduce);

    return matTp;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_tp_o(mat_t *mat, mat_t *matTp)
{
    /*

        Transpose an ordinary matrix

    */

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Index of location in mat's memory */
            size_t loc[2] = {i, j};
            size_t index = 0;
            int success = mat_mget_index(mat, loc, &index);

            /* Index of location in matTp's memory */
            size_t locTp[2] = {j, i};
            size_t indexTp = 0;
            int successTp = mat_mget_index(matTp, locTp, &indexTp);

            /* Element does not exist */
            if (success == -1 || successTp == -1)
                continue;

            /* Copy the element */
            ((double*) matTp -> matrix)[indexTp] = ((double*) mat -> matrix)[index];
          }
      }

    return 0;
}


static int _mat_tp_b(mat_t *mat, mat_t *matTp)
{
    /*

        Transpose a block matrix

    */

    /* Transpose blocks */

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Index of location in mat's memory */
            size_t loc[2] = {i, j};
            size_t indexMem = 0;
            int success = mat_mget_index(mat, loc, &indexMem);

            /* Block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* Location of block in matTp */
            size_t locTp[2] = {j, i};

            /* Block to be transposed */
            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Transpose the block */
            mat_t *matTpBlock = mat_tp(matBlock, NULL);

            /* Insert the block */
            mat_fset_block(matTp, locTp, matTpBlock);

            /* Free memory */
            matBlock = mat_ffree(matBlock);
          }
      }

    /* Free memory */
    free(visitedIndices);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_ftp_b(mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_ftp(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        (Fast) Transpose a matrix; Transposes matrix via index / location trafo

    */

    /* Check for NULL */
    if (mat == NULL)
        return NULL;

    /* Dimensions of transposed matrix */
    misc_swap(&mat -> dim[0], &mat -> dim[1], "ld");

    /* 'Transpose' the matrix by adding a transpose location trafo */
    mat_ltrfs_add(mat -> mltrfs, "tp");

    /* 'Transpose' blocks in block matrix */
    if (mat -> mblock != NULL)
      {
        _mat_ftp_b(mat);
      }

    /* Try to reduce mat */
    mat = mat_reduce(mat, reduce);

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_ftp_b(mat_t *mat)
{
    /*

        Transpose a block matrix

    */

    /* Dimensions */
    misc_swap(&mat -> mblock -> bdim[0], &mat -> mblock -> bdim[1], "*ld");
    misc_swap(&mat -> mblock -> bfdim[0], &mat -> mblock -> bfdim[1], "*ld");
    misc_swap(&mat -> mblock -> fdim[0], &mat -> mblock -> fdim[1], "ld");

    /* 'Transpose blocks' */

    /* Keep track of visited indices */
    size_t visitedIndicesSize = 0;
    double *visitedIndices = NULL;
    bool successInsert;

    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        for (size_t j = 0; j < mat -> dim[1]; j++)
          {
            /* Index of location in mat's memory */
            size_t loc[2] = {i, j};
            size_t indexMem = 0;
            int success = mat_mget_index(mat, loc, &indexMem);

            /* Block does not exist */
            if (success == -1)
                continue;

            /* Search for and insert the index (as double) in visitedIndices */
            misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

            /* Skip if an index has already been visited */
            if (!successInsert)
                continue;

            /* Block */
            mat_t *matBlock = ((mat_t**) mat -> matrix)[indexMem];

            /* Transpose the block */
            matBlock = mat_ftp(matBlock, NULL);
          }
      }

    /* Free memory */
    free(visitedIndices);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Add / Subtract Matrices   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mat_add_o(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd);
static int _mat_add_b(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd);
static int _mat_add_g(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_add(char job, mat_t *mat1, mat_t *mat2, mat_reduce_t *reduce)
{
    /*

        Add / subtract two matrices (mat1, mat2) depending on whether job is 'a' or 's'.

    */

    /* Must have job == 'a' or 's' */
    if (job != 'a' && job != 's')
      {
        printf("'job' must be either 'a' or 's' but have '%c' instead.\n", job);
        exit(1);

        return NULL;
      }

    /* Can only add / subtract matrices that have the same dimensions */
    size_t fdim1[2], fdim2[2];
    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);

    if (fdim1[0] != fdim2[0] || fdim1[1] != fdim2[1])
      {
        printf("Cannot add / subtract two matrices with total dimensions (%ld,%ld) and (%ld,%ld).\n", fdim1[0], fdim1[1], fdim2[0], fdim2[1]);
        exit(1);

        return NULL;
      }

    /* Result */
    mat_t *matAdd = NULL;

    /* Group mat1 and mat2 together */
    mats_t *mats = mats_new(0);

    mats_set_mat(mats, 0, mat1);
    mats_set_mat(mats, 1, mat2);

    /* Ordinary Matrices */
    if (mat1 -> mblock == NULL && mat2 -> mblock == NULL)
      {
        matAdd = mat_cp_structs(mats, NULL, 0); // 0 : Copy the similar ids of mats

        _mat_add_o(job, mat1, mat2, matAdd);
      }

    /* Block Matrices with the same Block Structure */
    else if ((mat1 -> mblock != NULL && mat2 -> mblock != NULL) && mat_comp_bstruct(mat1, mat2))
      {
        matAdd = mat_cp_structs(mats, NULL, 0); // 0 : Copy the similar ids of mats

        _mat_add_b(job, mat1, mat2, matAdd);
      }

    /* Any other Matrices */
    else
      {
        matAdd = mat_cp_structs(mats, _matTypeF.id, 0); // 0 : Set all ids to "f"

        _mat_add_g(job, mat1, mat2, matAdd);
      }

    /* Try to reduce the Matrix */
    matAdd = mat_reduce(matAdd, reduce);

    /* Free memory */
    mats = mats_free(mats);

    return matAdd;
}


/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_add_ip(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd, mat_reduce_t *reduce)
{
    /*

        Add / subtract two matrices (mat1, mat2) depending on whether job is 'a' or 's' in place (matAdd).

    */

    /* Must have job == 'a' or 's' */
    if (job != 'a' && job != 's')
      {
        printf("'job' must be either 'a' or 's' but have '%c' instead.\n", job);
        exit(1);

        return NULL;
      }

    /* Can only add / subtract matrices that have the same dimensions */
    size_t fdim1[2], fdim2[2], fdimAdd[2];
    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);
    mat_get_fdim(matAdd, fdimAdd);

    if (fdim1[0] != fdim2[0] || fdim1[1] != fdim2[1])
      {
        printf("Cannot add / subtract two matrices with total dimensions (%ld,%ld) and (%ld,%ld).\n", fdim1[0], fdim1[1], fdim2[0], fdim2[1]);
        exit(1);

        return NULL;
      }

    if (fdim1[0] != fdimAdd[0] || fdim1[1] != fdimAdd[1])
      {
        printf("Cannot add / subtract two matrices with total dimensions (%ld,%ld) and store the result in a matrix with total dimensions (%ld,%ld).\n", fdim1[0], fdim1[1], fdimAdd[0], fdimAdd[1]);
        exit(1);

        return NULL;
      }

    /* Add the Matrices */
    _mat_add_g(job, mat1, mat2, matAdd);

    /* Try to reduce the matrix (reduces everything -> no need to reduce in previous functions) */
    matAdd = mat_reduce(matAdd, reduce);

    return matAdd;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mat_add_o(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd)
{
    /*

        Add / subtract two ordinary matrices

    */

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < matAdd -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds[2] = {0, 0};
        mat_get_cbounds(matAdd, i, cBounds);

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            /* Value at location in mat1's nd mat2's matrices */
            size_t loc[2] = {i, j};
            double val1 = mat_get_value(mat1, loc);
            double val2 = mat_get_value(mat2, loc);

            /* Add values */
            double valAdd = (job == 'a') ? val1 + val2 : val1 - val2;

            /* Insert the new value into valAdd */
            mat_set_value(matAdd, loc, valAdd);
          }
      }

    return 0;
}


static int _mat_add_b(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd)
{
    /*

        Add two block matrices

    */

    for (size_t i = 0; i < matAdd -> dim[0]; i++)
      {
        /* Column bounds */
        size_t cBounds[2] = {0, 0};
        mat_get_cbounds(matAdd, i, cBounds);

        for (size_t j = cBounds[0]; j < cBounds[1]; j++)
          {
            /* Block at location in mat1's, mat2's and matAdd's matrices */
            size_t loc[2] = {i, j};

            mat_t *mat1Block = mat_fget_block(mat1, loc);
            mat_t *mat2Block = mat_fget_block(mat2, loc);
            mat_t *matAddBlock = mat_fget_block(matAdd, loc);

            if (matAddBlock == NULL)
              {
                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);

                continue;
              }

            if (mat1Block == NULL && mat2Block == NULL)
              {
                matAddBlock = mat_ffree(matAddBlock);

                continue;
              }

            if (mat1Block == NULL || mat2Block == NULL)
              {
                /* Add / Subtract a null Matrix to mat1Block or mat2Block */
                /* TODO: Can speed this up if mat_swap(mat1Block/mat2Block, matAddBlock) is called, and matrix element sign inversion is introduced as trafo */
                size_t fdimNull[2] = {0, 0};
                mat_get_fdim(matAddBlock, fdimNull);

                mat_t *matNull = mat_new(_matTypeNull.id, fdimNull, false);
                matAddBlock = mat_add_ip((mat1Block == NULL) ? job : 'a', matNull, (mat1Block == NULL) ? mat2Block : mat1Block, matAddBlock, NULL);
                matNull = mat_free(matNull);

                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);
                matAddBlock = mat_ffree(matAddBlock);

                continue;
              }

            /* Add the Blocks */
            matAddBlock = mat_add_ip(job, mat1Block, mat2Block, matAddBlock, NULL);

            /* Free memory */
            mat1Block = mat_ffree(mat1Block);
            mat2Block = mat_ffree(mat2Block);
            matAddBlock = mat_ffree(matAddBlock);
          }
      }

    return 0;
}


static int _mat_add_g(char job, mat_t *mat1, mat_t *mat2, mat_t *matAdd)
{
    /*

        Add / subtract two (general) matrices

    */

    /* Dimensions */
    size_t fdim[2] = {0, 0};
    mat_get_fdim(matAdd, fdim);

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < fdim[0]; i++)
      {
        for (size_t j = 0; j < fdim[1]; j++)
          {
            /* Value at location in mat1's nd mat2's matrices */
            size_t loc[2] = {i, j};
            double val1 = mat_get_value(mat1, loc);
            double val2 = mat_get_value(mat2, loc);

            /* Add values */
            double valAdd = (job == 'a') ? val1 + val2 : val1 - val2;

            /* Insert the new value into matAdd */
            mat_set_value(matAdd, loc, valAdd);
          }
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Multiply Matrices   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_mult_o(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value, bool ip);
static mat_t *_mat_mult_b(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value, bool ip);

static mat_t *_mat_mult_g(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_mult(mat_t *mat1, mat_t *mat2, double value, mat_reduce_t *reduce)
{
    /*

        Multiply two matrices (mat1, mat2) and a scalar

        TODO; Generalise this to n matrices

    */

    /* Can only multiply matrices if total row dimension of mat1 is equal to the total column dimension of mat2 */
    size_t fdim1[2], fdim2[2];
    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);

    if (fdim1[1] != fdim2[0])
      {
        printf("Cannot multiply two matrices with total dimensions (%ld,%ld) and (%ld,%ld).\n", fdim1[0], fdim1[1], fdim2[0], fdim2[1]);

        exit(1);

        return NULL;
      }


    /* Multiply the two matrices */

    mat_t *matMult = NULL;

    /* 'null' matrices */
    if (mat1 -> mtype == &_matTypeNull || mat2 -> mtype == &_matTypeNull)
      {
        size_t dimMult[2] = {fdim1[0], fdim2[1]};
        matMult = mat_new(_matTypeNull.id, dimMult, false);
      }

    /* Ordinary Matrices */
    else if (mat1 -> mblock == NULL && mat2 -> mblock == NULL)
      {
        matMult = _mat_mult_o('0', mat1, mat2, matMult, value, false);
      }

    /* Block Matrices with the same Block Structure of mat1's Rows and mat2's Columns */
    else if ((mat1 -> mblock != NULL && mat2 -> mblock != NULL) && mat_comp_bstruct_cr(mat1, mat2))
      {
        matMult = _mat_mult_b('0', mat1, mat2, matMult, value, false);
      }

    /* Any other Matrices */
    else
      {
        size_t fdim1[2] = {0, 0};
        size_t fdim2[2] = {0, 0};

        mat_get_fdim(mat1, fdim1);
        mat_get_fdim(mat2, fdim2);

        size_t dimMult[2] = {fdim1[0], fdim2[1]};

        matMult = mat_new(_matTypeF.id, dimMult, false);
        _mat_mult_g('0', mat1, mat2, matMult, value);
      }

    /* Try to reduce the matric */
    matMult = mat_reduce(matMult, reduce);

    return matMult;
}


/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_mult_ip(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value, mat_reduce_t *reduce)
{
    /*

        Multiply two matrices (mat1, mat2) in place (matMult) by adding (job == 'a') or substracting (job == 's') elements
        from matMult. If job == '0' == 0, elements in matMult are replaced by result of the matrix multiplication of mat1 and mat2.

        The final result of matMult is multiplied with value.

    */

    /* Check for NULL */
    if (mat1 == NULL || mat2 == NULL)
      {
        /* Overwrite matMult with zeros */
        if (job == '0')
          {
            size_t fdimMult[2];
            mat_get_fdim(matMult, fdimMult);

            mat_t *matNull = mat_new(_matTypeNull.id, fdimMult, false);
            mat_set_matrix(matMult, matNull);

            matNull = mat_free(matNull);

            return matMult;
          }

        /* Multiply matMult with result */
        matMult = mat_mult_val_ip(matMult, value, NULL);

        return matMult;
      }

    /* Check dimensions */
    size_t fdim1[2], fdim2[2], fdimMult[2];
    mat_get_fdim(mat1, fdim1);
    mat_get_fdim(mat2, fdim2);
    mat_get_fdim(matMult, fdimMult);

    if (fdim1[1] != fdim2[0])
      {
        printf("Cannot multiply two matrices with total dimensions (%ld,%ld) and (%ld,%ld).\n", fdim1[0], fdim1[1], fdim2[0], fdim2[1]);
        exit(1);

        return NULL;
      }

    if (fdim1[0] != fdimMult[0] || fdim2[1] != fdimMult[1])
      {
        printf("Cannot multiply two matrices with total dimensions (%ld,%ld) and (%ld,%ld) and store the result in a matrix with total dimensions (%ld,%ld).\n", fdim1[0], fdim1[1], fdim2[0], fdim2[1], fdimMult[0], fdimMult[1]);
        exit(1);

        return NULL;
      }

    /* Same addresses */
    if (matMult == mat1 && matMult == mat2)
      {
        /* Need new Matrix */
        mat_t *matMultRes = mat_cp(matMult);

        /* Multiply matrices */
        matMultRes = mat_mult_ip(job, mat1, mat2, matMultRes, value, NULL);

        /* Swap Matrices and free old Matrix */
        mat_swap(matMult, matMultRes);
        matMultRes = mat_free(matMultRes);
      }

    else
      {
        /* Ordinary Matrices */
        if (mat1 -> mblock == NULL && mat2 -> mblock == NULL)
          {
            matMult = _mat_mult_o(job, mat1, mat2, matMult, value, true);
          }

        /* Block Matrices with the same Block Structure of mat1's Rows and mat2's Columns */
        else if ((mat1 -> mblock != NULL && mat2 -> mblock != NULL) && mat_comp_bstruct_cr(mat1, mat2))
          {
            matMult = _mat_mult_b(job, mat1, mat2, matMult, value, true);
          }

        /* Any other Matrices */
        else
          {
            matMult = _mat_mult_g(job, mat1, mat2, matMult, value);
          }
      }

    /* Try to reduce the Matrix */
    matMult = mat_reduce(matMult, reduce);

    return matMult;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_mult_o(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value, bool ip)
{
    /*

        Multiply two ordinary matrices and a scalar

    */

    /* Dimensions */
    size_t dim1[2] = {mat1 -> dim[0], mat1 -> dim[1]};
    size_t dim2[2] = {mat2 -> dim[0], mat2 -> dim[1]};
    size_t dimMult[2] = {dim1[0], dim2[1]};


    /* Diagonal Matrix */

    if (mat1 -> mtype == &_matTypeD)
      {
        /* New matrix */
        matMult = (ip) ? matMult : mat_new((char*) ((mat2 -> mtype -> category == 0) ? mat2 -> mtype -> id : _matTypeF.id), dimMult, false);

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < dim2[0]; i++)
          {
            /* Value at location in mat1's matrix */
            size_t loc1[2] = {i, i};
            double val1 = mat_get_value(mat1, loc1);

            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(mat2, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Value at location in mat2's matrix */
                size_t loc2[2] = {i, j};
                double val2 = mat_get_value(mat2, loc2);

                /* Value does not exist */
                if (val2 == 0.)
                    continue;

                /* Insert the value into matMult */
                if (job == 'a')
                    mat_set_value(matMult, loc2, value * (mat_get_value(matMult, loc2) + val1 * val2));

                else if (job == 's')
                    mat_set_value(matMult, loc2, value * (mat_get_value(matMult, loc2) - val1 * val2));

                else
                    mat_set_value(matMult, loc2, value * val1 * val2);
              }
          }

        return matMult;
      }

    if (mat2 -> mtype == &_matTypeD)
      {
        /* New matrix */
        matMult = (ip) ? matMult : mat_new((char*) ((mat1 -> mtype -> category == 0) ? mat1 -> mtype -> id : _matTypeF.id), dimMult, false);

        #pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < dim1[1]; j++)
          {
            /* Value at location in mat2's matrix */
            size_t loc2[2] = {j, j};
            double val2 = mat_get_value(mat2, loc2);

            /* Row bounds */
            size_t rBounds[2] = {0, 0};
            mat_get_rbounds(mat1, j, rBounds);

            for (size_t i = rBounds[0]; i < rBounds[1]; i++)
              {
                /* Value at location in mat1's matrix */
                size_t loc1[2] = {i, j};
                double val1 = mat_get_value(mat1, loc1);

                /* Value does not exist */
                if (val1 == 0.)
                    continue;

                /* Insert the value into matMult */
                if (job == 'a')
                    mat_set_value(matMult, loc1, value * (mat_get_value(matMult, loc1) + val1 * val2));

                else if (job == 's')
                    mat_set_value(matMult, loc1, value * (mat_get_value(matMult, loc1) - val1 * val2));

                else
                    mat_set_value(matMult, loc1, value * val1 * val2);
              }
          }

        return matMult;
      }


    /* Other types */

    /* New matrix */
    matMult = (ip) ? matMult : mat_new(_matTypeF.id, dimMult, false);

    /* Array to store j'th row of matMult */
    double *valMult = malloc(sizeof(double) * dimMult[1]);

    /* Multiply the matrices */
    #pragma omp parallel
    #pragma omp single
    for (size_t i = 0; i < dimMult[0]; i++)
      {
        #pragma omp taskloop
        for (size_t j = 0; j < dimMult[1]; j++)
          {
            /* Value at location in matMult's matrix */
            valMult[j] = 0.;

            for (size_t k = 0; k < dim1[1]; k++)
              {
                /* Value at location in mat1's matrix */
                size_t loc1[2] = {i, k};
                double val1 = mat_get_value(mat1, loc1);

                /* Value does not exist */
                if (val1 == 0.)
                    continue;

                /* Value at location in mat2's matrix */
                size_t loc2[2] = {k, j};
                double val2 = mat_get_value(mat2, loc2);

                /* Value does not exist */
                if (val2 == 0.)
                    continue;

                /* Increment valMult */
                valMult[j] += val1 * val2;
              }
          }

        #pragma omp taskloop
        for (size_t j = 0; j < dimMult[1]; j++)
          {
            /* Location in matMult#s matrix */
            size_t locMult[2] = {i, j};

            /* Insert the value into matMult */
            if (job == 'a')
                mat_set_value(matMult, locMult, value * (mat_get_value(matMult, locMult) + valMult[j]));

            else if (job == 's')
                mat_set_value(matMult, locMult, value * (mat_get_value(matMult, locMult) - valMult[j]));

            else
                mat_set_value(matMult, locMult, value * valMult[j]);
          }
      }

    /* Free memory */
    free(valMult);

    return matMult;
}


static mat_t *_mat_mult_b(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value, bool ip)
{
    /*

        Multiply two block matrices

    */

    /* Dimensions */
    size_t dim1[2] = {mat1 -> dim[0], mat1 -> dim[1]};
    size_t dim2[2] = {mat2 -> dim[0], mat2 -> dim[1]};
    size_t dimMult[2] = {dim1[0], dim2[1]};


    /* Diagonal Matrix */

    if (mat1 -> mtype == &_matTypeD)
      {
        /* New matrix */
        matMult = (ip) ? matMult : mat_new((char*) ((mat2 -> mtype -> category == 0) ? mat2 -> mtype -> id : _matTypeF.id), dimMult, true);

        for (size_t i = 0; i < dim2[0]; i++)
          {
            /* Block at location in mat1's matrix */
            size_t loc1[2] = {i, i};
            mat_t *mat1Block = mat_fget_block(mat1, loc1);

            /* Block does not exist */
            if (mat1Block == NULL)
                continue;

            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(mat2, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Block at location in mat2's matrix */
                size_t loc2[2] = {i, j};
                mat_t *mat2Block = mat_fget_block(mat2, loc2);

                /* Block does not exist */
                if (mat2Block == NULL)
                    continue;

                if (ip)
                  {
                    /* New block */
                    mat_t *matMultBlock = mat_fget_block(matMult, loc2);

                    /* Multiply blocks */
                    matMultBlock = mat_mult_ip(job, mat1Block, mat2Block, matMultBlock, value, NULL);

                    /* Free memory */
                    matMultBlock = mat_ffree(matMultBlock);
                  }

                else
                  {
                    /* New block */
                    mat_t *matMultBlock = mat_mult(mat1Block, mat2Block, value, NULL);

                    /* Set the block in matMult */
                    mat_fset_block(matMult, loc2, matMultBlock);
                  }

                /* Free memory */
                mat2Block = mat_ffree(mat2Block);
              }

            /* Free memory */
            mat1Block = mat_ffree(mat1Block);
          }

        return matMult;
      }

    if (mat2 -> mtype == &_matTypeD)
      {
        /* New matrix */
        matMult = (ip) ? matMult : mat_new((char*) ((mat1 -> mtype -> category == 0) ? mat1 -> mtype -> id : _matTypeF.id), dimMult, true);

        for (size_t j = 0; j < dim1[1]; j++)
          {
            /* Block at location in mat2's matrix */
            size_t loc2[2] = {j, j};
            mat_t *mat2Block = mat_fget_block(mat2, loc2);

            /* Block does not exist */
            if (mat2Block == NULL)
                continue;

            /* Row bounds */
            size_t rBounds[2] = {0, 0};
            mat_get_rbounds(mat1, j, rBounds);

            for (size_t i = rBounds[0]; i < rBounds[1]; i++)
              {
                /* Block at location in mat1's matrix */
                size_t loc1[2] = {i, j};
                mat_t *mat1Block = mat_fget_block(mat1, loc1);

                /* Block does not exist */
                if (mat1Block == NULL)
                    continue;

                if (ip)
                  {
                    /* New block */
                    mat_t *matMultBlock = mat_fget_block(matMult, loc1);

                    /* Multiply blocks */
                    matMultBlock = mat_mult_ip(job, mat1Block, mat2Block, matMultBlock, value, NULL);

                    /* Free memory */
                    matMultBlock = mat_ffree(matMultBlock);
                  }

                else
                  {
                    /* New block */
                    mat_t *matMultBlock = mat_mult(mat1Block, mat2Block, value, NULL);

                    /* Set the block in matMult */
                    mat_fset_block(matMult, loc1, matMultBlock);
                  }

                /* Free memory */
                mat1Block = mat_ffree(mat1Block);
              }

            /* Free memory */
            mat2Block = mat_ffree(mat2Block);
          }

        return matMult;
      }


    /* Other types */

    /* New matrix */
    matMult = (ip) ? matMult : mat_new(_matTypeF.id, dimMult, true);

    for (size_t i = 0; i < dimMult[0]; i++)
      {
        for (size_t j = 0; j < dimMult[1]; j++)
          {
            /* Block at location in matMult's matrix */
            size_t locMult[2] = {i, j};
            mat_t *matMultBlock = (ip) ? mat_fget_block(matMult, locMult) : NULL;

            for (size_t k = 0; k < dim1[1]; k++)
              {
                /* Block at location in mat1's matrix */
                size_t loc1[2] = {i, k};
                mat_t *mat1Block = mat_fget_block(mat1, loc1);

                /* Block at location in mat2's matrix */
                size_t loc2[2] = {k, j};
                mat_t *mat2Block = mat_fget_block(mat2, loc2);

                /* Block does not exist */
                if (mat1Block == NULL || mat2Block == NULL)
                  {
                    mat1Block = mat_ffree(mat1Block);
                    mat2Block = mat_ffree(mat2Block);

                    continue;
                  }

                if (ip)
                  {
                    /* Multiply blocks */
                    matMultBlock = mat_mult_ip(job, mat1Block, mat2Block, matMultBlock, 1., NULL);
                  }

                else
                  {
                    /* Multiply blocks */
                    mat_t *matMultBlock_ = mat_mult(mat1Block, mat2Block, value, NULL);

                    /* Matrix not yet initialised */
                    if (matMultBlock == NULL)
                      {
                        matMultBlock = mat_cp_struct(matMultBlock_, _matTypeF.id, 0); // 0 : Set all ids to "f"
                        mat_set_matrix(matMultBlock, matMultBlock_);
                      }

                    /* Increment matMultBlock */
                    else
                      {
                        matMultBlock = mat_add_ip('a', matMultBlock, matMultBlock_, matMultBlock, NULL);
                      }

                    /* Free memory */
                    matMultBlock_ = mat_free(matMultBlock_);
                  }

                /* Free memory */
                mat1Block = mat_ffree(mat1Block);
                mat2Block = mat_ffree(mat2Block);
              }

            if (ip)
              {
                /* Multiply block with value */
                matMultBlock = mat_mult_val_ip(matMultBlock, value, NULL);
              }

            else
              {
                /* Matrix not yet initialised */
                if (matMultBlock == NULL)
                  {
                    /* Depth and location of block */
                    size_t depth = 0;

                    size_t **loc = malloc(sizeof(size_t*) * (depth + 1));
                    loc[0] = malloc(sizeof(size_t) * 2);

                    loc[0][0] = i;
                    loc[0][1] = j;

                    /* Get the block */
                    matMultBlock = mat_cp_struct_block_rc(mat1, mat2, depth, loc);

                    /* Free memory */
                    free(loc[0]);
                    free(loc);
                  }

                /* Insert matMultBlock into matMult */
                mat_fset_block(matMult, locMult, matMultBlock);
              }
          }
      }

    return matMult;
}


static mat_t *_mat_mult_g(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, double value)
{
    /*

        Multiply two (general) matrices

    */

    /* Dimensions */
    size_t fdim1[2] = {0, 0};
    mat_get_fdim(mat1, fdim1);

    size_t fdimMult[2] = {0, 0};
    mat_get_fdim(matMult, fdimMult);


    /* Write to matMult = mat1 or matMult != mat1, mat2 */
    if (matMult == mat1 || matMult != mat2)
      {
        /* Array to store j'th row of matMult */
        double *valMult = malloc(sizeof(double) * fdimMult[1]);

        /* Multiply the matrices */
        #pragma omp parallel
        for (size_t i = 0; i < fdimMult[0]; i++)
          {
            #pragma omp for schedule(dynamic)
            for (size_t j = 0; j < fdimMult[1]; j++)
              {
                valMult[j] = 0.;

                for (size_t k = 0; k < fdim1[1]; k++)
                  {
                    /* Value at location in mat1's matrix */
                    size_t loc1[2] = {i, k};
                    double val1 = mat_get_value(mat1, loc1);

                    /* Value does not exist */
                    if (val1 == 0.)
                        continue;

                    /* Value at location in mat2's matrix */
                    size_t loc2[2] = {k, j};
                    double val2 = mat_get_value(mat2, loc2);

                    /* Value does not exist */
                    if (val2 == 0.)
                        continue;

                    /* Increment valMult */
                    valMult[j] += val1 * val2;
                  }
              }

            #pragma omp for
            for (size_t j = 0; j < fdimMult[1]; j++)
              {
                /* Location in matMult#s matrix */
                size_t locMult[2] = {i, j};

                /* Insert the value into matMult */
                if (job == 'a')
                    mat_set_value(matMult, locMult, value * (mat_get_value(matMult, locMult) + valMult[j]));

                else if (job == 's')
                    mat_set_value(matMult, locMult, value * (mat_get_value(matMult, locMult) - valMult[j]));

                else
                    mat_set_value(matMult, locMult, value * valMult[j]);
              }
          }

        /* Free memory */
        free(valMult);
      }

    /* Write to matMult = mat2 */
    else
      {
        /* Array to store i'th column of matMult */
        double *valMult = malloc(sizeof(double) * fdimMult[0]);

        /* Multiply the matrices */
        #pragma omp parallel
        for (size_t j = 0; j < fdimMult[1]; j++)
          {
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < fdimMult[0]; i++)
              {
                valMult[i] = 0.;

                for (size_t k = 0; k < fdim1[1]; k++)
                  {
                    /* Value at location in mat1's matrix */
                    size_t loc1[2] = {i, k};
                    double val1 = mat_get_value(mat1, loc1);

                    /* Value does not exist */
                    if (val1 == 0.)
                        continue;

                    /* Value at location in mat2's matrix */
                    size_t loc2[2] = {k, j};
                    double val2 = mat_get_value(mat2, loc2);

                    /* Value does not exist */
                    if (val2 == 0.)
                        continue;

                    /* Increment valMult */
                    valMult[i] += val1 * val2;
                  }
              }

            #pragma omp for
            for (size_t i = 0; i < fdimMult[0]; i++)
              {
                /* Location in matMult#s matrix */
                size_t locMult[2] = {i, j};

                /* Insert the value into matMult */
                if (job == 'a')
                    mat_set_value(matMult, locMult, value * (mat_get_value(matMult, locMult) + valMult[i]));

                else if (job == 's')
                    mat_set_value(matMult, locMult, value * (mat_get_value(matMult, locMult) - valMult[i]));

                else
                    mat_set_value(matMult, locMult, value * valMult[i]);
              }
          }

        /* Free memory */
        free(valMult);
      }

    return matMult;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_mult_mat(mat_t *mat1, mat_t *mat2, mat_reduce_t *reduce)
{
    /*

        Multiply two matrices

    */

    /* Multiply the two matrices */
    mat_t *matMult = mat_mult(mat1, mat2, 1., reduce);

    return matMult;
}


mat_t *mat_mult_mat_ip(char job, mat_t *mat1, mat_t *mat2, mat_t *matMult, mat_reduce_t *reduce)
{
    /*

        Multiply two matrices in place

    */

    /* Multiply the two matrices */
    matMult = mat_mult_ip(job, mat1, mat2, matMult, 1., reduce);

    return matMult;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_mult_val_o(mat_t *mat, double val);
static mat_t *_mat_mult_val_b(mat_t *mat, double val);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_mult_val(mat_t *mat, double value, mat_reduce_t *reduce)
{
    /*

        Multiply a matrix with a scalar

    */

    mat_t *matMult = mat_mult_val_ip(mat_cp(mat), value, reduce);

    return matMult;
}


mat_t *mat_mult_val_ip(mat_t *mat, double value, mat_reduce_t *reduce)
{
    /*

        Multiply a matrix with a scalar by overwriting mat

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        mat = _mat_mult_val_o(mat, value);
      }

    /* Block Matrix */
    else
      {
        mat = _mat_mult_val_b(mat, value);
      }

    /* Attempt to reduce the Matrix */
    mat = mat_reduce(mat, reduce);

    return mat;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_mult_val_o(mat_t *mat, double value)
{
    /*

        Multiply an ordinary matrix with a scalar

    */

    for (size_t i = 0; i < mat -> dim[2]; i++)
      {
        ((double*) mat -> matrix)[i] *= value;
      }

    return mat;
}


static mat_t *_mat_mult_val_b(mat_t *mat, double value)
{
    /*

        Multiply a block matrix with a scalar

    */

    for (size_t i = 0; i < mat -> dim[2]; i++)
      {
        mat_mult_val(((mat_t**) mat -> matrix)[i], value, NULL);
      }

    return mat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Permutate Matrix   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_perm_o(char job, mat_t *mat, size_t *perm);
static mat_t *_mat_perm_b(char job, mat_t *mat, size_t *perm);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_perm(char job, mat_t *mat, size_t *perm, mat_reduce_t *reduce)
{
    /*

        Perform a permutation of mat's rows if job == 'r' or mat's columns if job == 'c'

    */

    /* Permutated Matrix */
    mat_t *matPerm = NULL;


    /* Permutate the Matrix */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        matPerm = _mat_perm_o(job, mat, perm);
      }

    /* Block Matrix */
    else
      {
        matPerm = _mat_perm_b(job, mat, perm);
      }

    /* Attempt to reduce the Matrix */
    matPerm = mat_reduce(matPerm, reduce);

    return matPerm;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_perm_o(char job, mat_t *mat, size_t *perm)
{
    /*

        Permutate an ordinary matrix depending on job

    */

    /* Permutated Matrix */
    mat_t *matPerm = NULL;

    /* "null" matrix */
    if (mat -> mtype == &_matTypeNull)
      {
        matPerm = mat_new(mat -> mtype -> id, mat -> dim, false);

        return matPerm;
      }

    /* Rows */
    if (job == 'r')
      {
        /* Allocate mmeory */
        const char *type = (mat -> mtype == &_matTypeD) ? _matTypeD.id : _matTypeF.id;
        matPerm = mat_new(type, mat -> dim, false);

        /* Permutate Matrix */
        for (size_t i = 0; i < matPerm -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matPerm, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Locations */
                size_t locPerm[2] = {i, j};
                size_t loc[2] = {perm[i], j};

                /* Set the new value */
                mat_set_value(matPerm, locPerm, mat_get_value(mat, loc));
              }
          }

        return matPerm;
      }

    /* Columns */
    if (job == 'c')
      {
        /* Allocate mmeory */
        const char *type = (mat -> mtype == &_matTypeD) ? _matTypeD.id : _matTypeF.id;
        matPerm = mat_new(type, mat -> dim, false);

        /* Permutate Matrix */
        for (size_t i = 0; i < matPerm -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matPerm, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Locations */
                size_t locPerm[2] = {i, j};
                size_t loc[2] = {i, perm[j]};

                /* Set the new value */
                mat_set_value(matPerm, locPerm, mat_get_value(mat, loc));
              }
          }

        return matPerm;
      }


    /* Wrong job */
    printf("Cannot permutate an ordinary matrix for unknown job '%c'.\n", job);
    exit(1);

    return NULL;
}


static mat_t *_mat_perm_b(char job, mat_t *mat, size_t *perm)
{
    /*

        Permutate a block matrix depending on job

    */

    /* Permutated Matrix */
    mat_t *matPerm = NULL;

    /* Rows */
    if (job == 'r')
      {
        /* Allocate mmeory */
        const char *type = (mat -> mtype == &_matTypeD) ? _matTypeD.id : _matTypeF.id;
        matPerm = mat_new(type, mat -> dim, true);

        /* Permutate Matrix */
        for (size_t i = 0; i < matPerm -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matPerm, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Locations */
                size_t locPerm[2] = {i, j};
                size_t loc[2] = {perm[i], j};

                /* Matrix at location */
                mat_t *matBlock = mat_get_block(mat, loc);

                if (matBlock == NULL)
                  {
                    size_t **locBlock = malloc(sizeof(size_t*) * 1);
                    locBlock[0] = loc;

                    matBlock = mat_cp_struct_block(mat, 0, locBlock, _matTypeF.id, 0); // TODO: 0 : Set all ids to "null"

                    free(locBlock);
                  }

                /* Set the new value */
                mat_fset_block(matPerm, locPerm, matBlock);
              }
          }

        return matPerm;
      }

    /* Columns */
    if (job == 'c')
      {
        /* Allocate mmeory */
        const char *type = (mat -> mtype == &_matTypeD) ? _matTypeD.id : _matTypeF.id;
        matPerm = mat_new(type, mat -> dim, true);

        /* Permutate Matrix */
        for (size_t i = 0; i < matPerm -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matPerm, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Locations */
                size_t locPerm[2] = {i, j};
                size_t loc[2] = {i, perm[j]};

                /* Matrix at location */
                mat_t *matBlock = mat_get_block(mat, loc);

                if (matBlock == NULL)
                  {
                    size_t **locBlock = malloc(sizeof(size_t*) * 1);
                    locBlock[0] = loc;

                    matBlock = mat_cp_struct_block(mat, 0, locBlock, _matTypeNull.id, 0); // 0 : Set all ids to "null"

                    free(locBlock);
                  }

                /* Set the new value */
                mat_fset_block(matPerm, locPerm, matBlock);
              }
          }

        return matPerm;
      }


    /* Wrong job */
    printf("Cannot permutate a block matrix for unknown job '%c'.\n", job);
    exit(1);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Correlation Matrix   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_corr_o(mat_t *mat, bool ip);
static mat_t *_mat_corr_b(mat_t *mat, bool ip);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_corr(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Get the correlation matrix from the covariance matrix mat

    */

    /* Correlation Matrix */
    mat_t *matCorr = NULL;

    /* Check dimensions */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    if (fdim[0] != fdim[1])
      {
        printf("Cannot get the correlation matrix from a non-square matrix.\n");
        exit(1);

        return NULL;
      }

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        matCorr = _mat_corr_o(mat, false);
      }

    /* Block Matrix */
    else
      {
        matCorr = _mat_corr_b(mat, false);
      }

    /* Attempt to reduce the Matrix */
    matCorr = mat_reduce(matCorr, reduce);

    return matCorr;
}


mat_t *mat_corr_ip(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Get the correlation matrix from the covariance matrix mat

    */

    /* Correlation Matrix */
    mat_t *matCorr = NULL;

    /* Check dimensions */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    if (fdim[0] != fdim[1])
      {
        printf("Cannot get the correlation matrix from a non-square matrix.\n");
        exit(1);

        return NULL;
      }

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        matCorr = _mat_corr_o(mat, true);
      }

    /* Block Matrix */
    else
      {
        matCorr = _mat_corr_b(mat, true);
      }

    /* Attempt to reduce the Matrix */
    matCorr = mat_reduce(matCorr, reduce);

    return matCorr;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_corr_o(mat_t *mat, bool ip)
{
    /*

        Get the correlation matrix from the ordinary covariance matrix mat

    */

    /* Correlation matrix */
    mat_t *matCorr = (ip) ? mat : mat_new(mat -> mtype -> id, mat -> dim, false);

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < mat -> dim[2]; i++)
      {
        /* Location at index */
        size_t loc[2];
        mat_get_loc(mat, i, loc);

        /* Do diagonal elements later */
        if (loc[0] == loc[1])
            continue;

        /* Get the diagonal elements */
        size_t locR[2] = {loc[0], loc[0]};
        double valR = sqrt(mat_get_value(mat, locR));

        size_t locC[2] = {loc[1], loc[1]};
        double valC = sqrt(mat_get_value(mat, locC));

        /* Set the correlation */
        mat_set_value(matCorr, loc, mat_get_value(mat, loc) / (valR * valC));
      }

    /* Set the diagonal elements */
    for (size_t i = 0; i < mat -> dim[0]; i++)
      {
        size_t loc[2] = {i, i};

        mat_set_value(matCorr, loc, 1.);
      }

    return matCorr;
}


static mat_t *_mat_corr_b(mat_t *mat, bool ip)
{
    /*

        Get the correlation matrix from the block covariance matrix mat

    */

    /* Correlation matrix */
    mat_t *matCorr = (ip) ? mat : mat_cp_struct(mat, NULL, 0); // 0 : Copy mat's ids

    /* (Inverse square root of) Diagonal matrix */
    mat_t *matDiag = mat_diag_func(mat_diag(mat), mat_diag_func_invsqrt);

    /* Correlation matrix in place */
    if (ip)
      {
        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            /* Correlation block at location */
            size_t loc[2];
            mat_get_loc(mat, i, loc);

            mat_t *matCorrBlock = mat_fget_block(matCorr, loc);

            /* Diagonal block */
            if (loc[0] == loc[1])
              {
                /* Get the correlation block */
                matCorrBlock = mat_corr_ip(matCorrBlock, NULL);

                /* Free memory */
                matCorrBlock = mat_ffree(matCorrBlock);

                continue;
              }

            /* Diagonal block at location */
            size_t locR[2] = {loc[0], loc[0]};
            mat_t *matDiagBlockR = mat_mget_block(matDiag, locR);

            size_t locC[2] = {loc[1], loc[1]};
            mat_t *matDiagBlockC = mat_mget_block(matDiag, locC);

            /* Get the correlation block */
            matCorrBlock = mat_mult_mat_ip('0', matDiagBlockR, matCorrBlock, matCorrBlock, NULL);
            matCorrBlock = mat_mult_mat_ip('0', matCorrBlock, matDiagBlockC, matCorrBlock, NULL);

            /* Free memory */
            matCorrBlock = mat_ffree(matCorrBlock);
          }
      }

    /* New correlation matrix */
    else
      {
        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            /* Correlation block at location */
            size_t loc[2];
            mat_get_loc(mat, i, loc);

            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Diagonal block */
            if (loc[0] == loc[1])
              {
                /* Get the correlation block */
                mat_t *matCorrBlock = mat_corr(matBlock, NULL);

                /* Set the block */
                mat_fset_block(matCorr, loc, matCorrBlock);

                /* Free memory */
                matBlock = mat_ffree(matBlock);

                continue;
              }

            /* Diagonal block at location */
            size_t locR[2] = {loc[0], loc[0]};
            mat_t *matDiagBlockR = mat_mget_block(matDiag, locR);

            size_t locC[2] = {loc[1], loc[1]};
            mat_t *matDiagBlockC = mat_mget_block(matDiag, locC);

            /* Get the correlation block */
            mat_t *matCorrBlock = mat_mult_mat(matDiagBlockR, matBlock, NULL);
            mat_mult_mat_ip('0', matCorrBlock, matDiagBlockC, matCorrBlock, NULL);

            /* Set the correlation block */
            mat_fset_block(matCorr, loc, matCorrBlock);
          }
      }

    /* Free memory */
    matDiag = mat_free(matDiag);

    return matCorr;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   QR Decomposition   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mats_t *_mat_qr_o(mat_t *mat);
static mats_t *_mat_qr_b(mat_t *mat);
static mats_t *_mat_qr_g(mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *mat_qr(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Decompose a (square) matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */

    /* Check for NULL */
    if (mat == NULL)
        return NULL;

    mats_t *matsQR;

    /* Ordinary matrix */
    if (mat -> mblock == NULL)
      {
        matsQR = _mat_qr_o(mat);
      }

    /* Block matrix */
    else
      {
        matsQR = _mat_qr_b(mat);
      }

    /* Reduce the matrices */
    mat_t *matQ = mats_get_mat(matsQR, 0);
    mat_t *matR = mats_get_mat(matsQR, 1);

    matQ = mat_reduce(matQ, reduce);
    matR = mat_reduce(matR, reduce);

    return matsQR;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int *_mat_qr_o_h(double *householderVector, double *matRArray, size_t col, size_t dim);
static int *_mat_qr_o_r(double *householderVector, double *matRArray, size_t col, size_t dim);
static int *_mat_qr_o_q(double *householderVector, double *matQArray, size_t col, size_t dim);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static mats_t *_mat_qr_o(mat_t *mat)
{
    /*

        Decompose a (square) ordinary matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */

    /* Diagonal or upper triangular matrix */
    if (mat -> mtype == &_matTypeD || !strcmp(mat -> mtype -> id, "ut"))
      {
        mat_t *matR = mat_cp(mat);
        mat_t *matQ = mat_new_id(mat);

        /* Return matrices as collection */
        size_t matsSize = 2;
        mats_t *matsQR = mats_new(matsSize);

        /* Insert the matrices */
        mats_set_mat(matsQR, 0, matQ);
        mats_set_mat(matsQR, 1, matR);

        return matsQR;
      }


    /* Other matrix type */

    /* Matrix dimensions */
    size_t *dim = mat -> dim;

    /* Matrices R and Q */
    mat_t *matR = mat_trafo(mat, _matTypeF.id, mat -> dim, false, NULL);
    mat_t *matQ = mat_new(_matTypeF.id, dim, false);

    /* Initialise Q */
    for (size_t i = 0; i < *dim; i++)
      {
        size_t loc[2] = {i, i};
        mat_set_value(matQ, loc, 1.);
      }

    /* Householder vector */
    double *householderVector = malloc(sizeof(double) * *dim);

    #pragma omp parallel shared(householderVector, matR, matQ)
    #pragma omp single
    for (size_t i = 0; i < *dim - 1; i++)
      {
        /* Householder Vector */
        #pragma omp task depend(out: householderVector)
        _mat_qr_o_h(householderVector, matR -> matrix, i, *dim);

        /* R matrix */
        #pragma omp task depend(in: householderVector) depend(inout: matR -> matrix)
        _mat_qr_o_r(householderVector, matR -> matrix, i, *dim);

        /* Q matrix */
        #pragma omp task depend(in: householderVector)
        _mat_qr_o_q(householderVector, matQ -> matrix, i, *dim);
      }

    /* Return matrices as collection */
    size_t matsSize = 2;

    mats_t *matsQR = mats_new(matsSize);
    mats_set_mat(matsQR, 0, matQ);
    mats_set_mat(matsQR, 1, matR);

    /* Free memory */
    free(householderVector);

    return matsQR;
}


static mats_t *_mat_qr_g(mat_t *mat)
{
    /*

        Decompose a (square) ordinary matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */

    /* Other matrix type */

    /* Matrix dimensions */
    size_t *dim = mat -> dim;

    /* Matrices R and Q */
    mat_t *matR = mat;
    mat_t *matQ = mat_new(_matTypeF.id, dim, false);

    /* Initialise Q */
    for (size_t i = 0; i < *dim; i++)
      {
        size_t loc[2] = {i, i};
        mat_set_value(matQ, loc, 1.);
      }

    /* Householder vector */
    double *householderVector = malloc(sizeof(double) * *dim);

    #pragma omp parallel shared(householderVector, matR, matQ)
    #pragma omp single
    for (size_t i = 0; i < *dim - 1; i++)
      {
        /* Householder Vector */
        #pragma omp task depend(out: householderVector)
        _mat_qr_o_h(householderVector, matR -> matrix, i, *dim);

        /* R matrix */
        #pragma omp task depend(in: householderVector) depend(inout: matR -> matrix)
        _mat_qr_o_r(householderVector, matR -> matrix, i, *dim);

        /* Q matrix */
        #pragma omp task depend(in: householderVector)
        _mat_qr_o_q(householderVector, matQ -> matrix, i, *dim);
      }

    /* Return matrices as collection */
    size_t matsSize = 2;

    mats_t *matsQR = mats_new(matsSize);
    mats_set_mat(matsQR, 0, matQ);
    mats_set_mat(matsQR, 1, matR);

    /* Free memory */
    free(householderVector);

    return matsQR;
}


static int *_mat_qr_o_h(double *householderVector, double *matRArray, size_t col, size_t dim)
{
    /*

        Decompose a (square) ordinary matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */

    /* Norm of the first column of each submatrix (i.e. R_ii of the upper triangular matrix R) */
    double rColNorm = 0.;

    /* Norm of the (raw) Householder vector */
    double householderNorm = 0.;

    /* Work in dim - col dimensions for the submatrices */
    for (size_t row = col; row < dim; row++)
      {
        /* Norm of the first column */
        rColNorm += matRArray[row * dim + col] * matRArray[row * dim + col];

        /* Householder vector */
        householderVector[row] = matRArray[row * dim + col];
      }

    /* Fix the sign of colNorm (i.e. the sign of R_ii) */
    if (matRArray[(dim + 1) * col] > 0.)
        rColNorm = -sqrt(rColNorm);
    else
        rColNorm = sqrt(rColNorm);

    /* Raw Householder vector */
    householderVector[col] -= rColNorm;

    /* Norm of the raw Householder vector */
    householderNorm = sqrt(rColNorm * (rColNorm - matRArray[(dim + 1) * col]));

    /* Normalise the Householder vector */
    for (size_t i = col; i < dim; i++)
        householderVector[i] /= householderNorm;

    return 0;
}


static int *_mat_qr_o_r(double *householderVector, double *matRArray, size_t col, size_t dim)
{
    /*

        Decompose a (square) ordinary matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */

    /* Work in dim - col dimensions for the submatrices */
    #pragma omp taskloop shared(householderVector, matRArray)
    for (size_t i = col; i < dim; i++)
      {
        /* Scalar product of the j'th column of the submatrix and the Householder vector */
        double rColHouseholderProduct = 0.;

        for (size_t j = col; j < dim; j++)
            rColHouseholderProduct += matRArray[dim * j + i] * householderVector[j];

        /* Next submatrix */
        for (size_t j = col; j < dim; j++)
            matRArray[dim * j + i] -= householderVector[j] * rColHouseholderProduct;
      }

    return 0;
}


static int *_mat_qr_o_q(double *householderVector, double *matQArray, size_t col, size_t dim)
{
    /*

        Decompose a (square) ordinary matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */

    /* Product between the Q matrix and the Householder vector */
    #pragma omp taskloop shared(matQArray, householderVector)
    for (size_t i = 0; i < dim; i++)
      {
        double qHouseholderProduct = 0.;

        /* Work in dim - col dimensions for the submatrices */
        for (size_t j = col; j < dim; j++)
          {
            qHouseholderProduct += matQArray[dim * i + j] * householderVector[j];
          }

        /* Work in dim - col dimensions for the submatrices */
        for (size_t j = col; j < dim; j++)
          {
            matQArray[dim * i + j] -= householderVector[j] * qHouseholderProduct;
          }
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mats_t *_mat_qr_b(mat_t *mat)
{
    /*

        Decompose a (square) block matrix into an orthogonal matrix Q and and upper triangular
        matrix R using Householder decompositions.

    */


    /* Diagonal matrix */

    if (mat -> mtype == &_matTypeD)
      {
        mat_t *matQ = mat_new(_matTypeD.id, mat -> dim, true);
        mat_t *matR = mat_new(_matTypeD.id, mat -> dim, true);

        for (size_t i = 0; i < *(mat -> dim); i++)
          {
            /* Block in mat's memory */
            size_t loc[2] = {i, i};
            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Check for NULL */
            if (matBlock == NULL)
                continue;

            /* Get the Q and R matrices of the i'th block */
            mats_t *matsQRBlock = mat_qr(matBlock, NULL);
            mat_t *matQBlock = mats_get_mat(matsQRBlock, 0);
            mat_t *matRBlock = mats_get_mat(matsQRBlock, 1);

            /* Insert the Q and R blocks */
            mat_fset_block(matQ, loc, matQBlock);
            mat_fset_block(matR, loc, matRBlock);

            /* Free memory */
            matBlock = mat_ffree(matBlock);
            matsQRBlock = mats_free(matsQRBlock);
          }

        /* Return matrices as collection */
        size_t matsSize = 2;

        mats_t *matsQR = mats_new(matsSize);
        mats_set_mat(matsQR, 0, matQ);
        mats_set_mat(matsQR, 1, matR);

        return matsQR;
      }


    /* Other type */

    /* Use a full ordinary matrix as input */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    mat_t *matFull = mat_trafo(mat, _matTypeF.id, fdim, false, NULL);

    /* Get the QR decomposition */
    mats_t *matsQRFull = _mat_qr_g(matFull);

    /* Transform Q and R into block matrices */
    mat_t *matQ = mat_cp_struct(mat, _matTypeF.id, 0); // 0 : Set all ids to "f"
    mat_t *matR = mat_cp_struct(mat, "ut", -1); // -1 : Only set the id of the highest block (ie of matR) to "ut" while setting all lower blocks' ids to "f"

    mat_set_matrix(matQ, mats_get_mat(matsQRFull, 0));
    mat_set_matrix(matR, mats_get_mat(matsQRFull, 1));

    matsQRFull = mats_free_full(matsQRFull);


    /* Return matrices */

    size_t matsSize = 2;

    mats_t *matsQR = mats_new(matsSize);
    mats_set_mat(matsQR, 0, matQ);
    mats_set_mat(matsQR, 1, matR);

    return matsQR;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Eigenvalue Decomposition   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mats_t *_mat_eig_o(mat_t *mat);
static mats_t *_mat_eig_b(mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *mat_eig(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Decompose a (square) symmetric matrix into an orthogonal matrix O of its eigenvectors, and into a
        diagonal matrix D of its eigenvalues.

    */

    /* Check for NULL */
    if (mat == NULL)
        return NULL;

    mats_t *matsEig = NULL;

    /* Ordinary matrix */
    if (mat -> mblock == NULL)
      {
        matsEig = _mat_eig_o(mat);
      }

    /* Block matrix */
    else
      {
        matsEig = _mat_eig_b(mat);
      }

    /* Reduce the matrices */
    mat_t *matD = mats_get_mat(matsEig, 0);
    mat_t *matO = mats_get_mat(matsEig, 1);

    matD = mat_reduce(matD, reduce);
    matO = mat_reduce(matO, reduce);

    return matsEig;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static mats_t *_mat_eig_o(mat_t *mat)
{
    /*

        Eigendecomposition of an ordinary matrix (mat) by using LAPACK routines

    */

    /* Eigendecomposition */
    mats_t *matsEig = mats_new(2);


    /* LAPACK VARIABLES */

    /* Dimensions of matrix */
    int laN = (int) mat -> dim[0];


    /* Leading dimension of A */
    int laLDA = laN;


    /* Size of WORK (should be larger than 3*N-1) */
    int laLWORK = 6 * laN;

    /* WORK */
    double *laWORK = malloc(sizeof(double) * ((size_t) laLWORK));


    /* INFO == 0 (success), == -i (i-th argument had illegal value, == i (algorithm failed...)) */
    int laINFO;


    /* Diagonal Matrix */
    if (mat -> mtype == &_matTypeD)
      {
        /* Eigenvalues */
        mat_t *matD = mat_cp(mat);

        /* Eigenvectors */
        mat_t *matO = mat_new_id(mat);

        /* Store matrices in matsEig */
        mats_set_mat(matsEig, 0, matD);
        mats_set_mat(matsEig, 1, matO);
      }

    /* Symmetric Matrix */
    else if (mat -> mtype == &_matTypeS)
      {
        /* Eigenvalues */
        mat_t *matD = mat_new(_matTypeD.id, mat -> dim, false);

        /* Eigenvectors */
        mat_t *matO = mat_trafo(mat, _matTypeF.id, mat -> dim, false, NULL);


        /* Job (eigenvalues and eigenvectors) */
        char laJOBZ = 'V';

        /* Storage of matrix (can choose either U or L) */
        char laUPLO = 'U';

        /* Input matrix + Output eigenvectors (LDA, N) */
        double *laA = matO -> matrix;

        /* Output eigenvalues */
        double *laW = matD -> matrix;

        /* Calculate the eigendexomoposition */
        dsyev_(&laJOBZ, &laUPLO, &laN, laA, &laLDA, laW, laWORK, &laLWORK, &laINFO);

        if (laINFO != 0)
          {
            matO = mat_free(matO);

            goto freeMemory;
          }

        /* Store matrices in matsEig */
        mats_set_mat(matsEig, 0, matD);
        mats_set_mat(matsEig, 1, matO);
      }

    /* Other Matrix Type */
    else
      {
        printf("Eigendecomposition of a matrix that is not diagonal or symmetric is not yet supported.\n");
        exit(1);

        return NULL;
      }


    freeMemory:

    /* Free memory */
    free(laWORK);

    return matsEig;
}


static mats_t *_mat_eig_b(mat_t *mat)
{
    /*

        Eigendecomposition of a block matrix (mat) by using LAPACK routines

    */

    /* Diagonal matrix */
    if (mat -> mtype == &_matTypeD)
      {
        mats_t *matsEig = mats_new(2);

        mat_t *matD = mat_new(_matTypeD.id, mat -> dim, true);
        mat_t *matO = mat_new(_matTypeD.id, mat -> dim, true);

        mats_set_mat(matsEig, 0, matD);
        mats_set_mat(matsEig, 1, matO);

        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            /* Block at location */
            size_t loc[2] = {i, i};
            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Eigendecomoposition of the block */
            mats_t *matsBlockEig = mat_eig(matBlock, NULL);

            if (matsBlockEig == NULL)
              {
                matO = mat_free(matO);
                matD = mat_free(matD);

                matBlock = mat_ffree(matBlock);

                return NULL;
              }

            /* Set the blocks */
            mat_fset_block(matD, loc, mats_get_mat(matsBlockEig, 0));
            mat_fset_block(matO, loc, mats_get_mat(matsBlockEig, 1));

            /* Free memory */
            matBlock = mat_ffree(matBlock);
            matsBlockEig = mats_free(matsBlockEig);
          }

        return matsEig;
      }


    /* Other type */

    /* Full dimensions */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    /* Full matrix */
    mat_t *matFull = mat_new(mat -> mtype -> id, fdim, false);
    mat_set_matrix(matFull, mat);

//    mat_trafo(mat, "s", fdim, false, NULL); // TODO: This trafo copies mat's type -> should however also consider blocks' types

    mats_t *matsEig = mat_eig(matFull, NULL);
    matFull = mat_free(matFull);

    if (matsEig == NULL)
      {
        return NULL;
      }

    /* Copy the struct of mat */
    mat_t *matD = mat_cp_struct(mat, _matTypeD.id, 0); // 0 : Set all ids to "d"
    mat_t *matO = mat_cp_struct(mat, _matTypeF.id, 0); // 0 : Set all ids to "f"

    /* Set the matrices */
    mat_set_matrix(matD, mats_get_mat(matsEig, 0));
    mat_free(mats_get_mat(matsEig, 0));

    mat_set_matrix(matO, mats_get_mat(matsEig, 1));
    mat_free(mats_get_mat(matsEig, 1));

    mats_set_mat(matsEig, 0, matD);
    mats_set_mat(matsEig, 1, matO);

    return matsEig;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%     LINEAR EQUATION SOLVER + MATRIX INVERSION     %%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Back Substitutiom   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_backsub_o(mat_t *matA, mat_t *matB);
static mat_t *_mat_backsub_b(mat_t *matA, mat_t *matB);
static mat_t *_mat_backsub_g(mat_t *matA, mat_t *matB);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_backsub(mat_t *matA, mat_t *matB, mat_reduce_t *reduce)
{
    /*

        Solve the linear equation

            A X = B

        for X, where A is an upper triangular matrix and B is a general matrix.

    */

    /* Check dimensions */
    size_t dimA[2], dimB[2];

    mat_get_fdim(matA, dimA);
    mat_get_fdim(matB, dimB);

    if (dimA[0] != dimB[0])
      {
        printf("Cannot solve the linear equation A X = B if A's dimensions are (%ld,%ld) and B's dimensions are (%ld,%ld).\n", dimA[0], dimA[1], dimB[0], dimB[1]);
        exit(1);

        return NULL;
      }

    /* Solution to the linear equation */
    mat_t *matX = NULL;

    /* Ordinary Matrices */
    if (matA -> mblock == NULL && matB -> mblock == NULL)
      {
        matX = _mat_backsub_o(matA, matB);
      }

    /* Block Matrices with matching Dimensions */
    else if ((matA -> mblock != NULL && matB -> mblock != NULL) && (matA -> dim[0] == matB -> dim[0]))
      {
        matX = _mat_backsub_b(matA, matB);
      }

    /* Any other Matrices */
    else
      {
        matX = _mat_backsub_g(matA, matB);
      }

    /* Try to reduce the matrix */
    matX = mat_reduce(matX, reduce);

    return matX;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_backsub_o(mat_t *matA, mat_t *matB)
{
    /*

        Solve the (ordinary) linear equation

            A X = B

        where A is an upper triangular matrix with back substitution.

        This function assumes A and B to be ordinary matrices.

    */

    /* Dimensions */
//    size_t dimA[2] = {matA -> dim[0], matA -> dim[1]};
    size_t dimB[2] = {matB -> dim[0], matB -> dim[1]};


    /* Diagonal matrix */

    if (matA -> mtype == &_matTypeD)
      {
        /* Only need to copy matB's struct */
        mat_t *matX = mat_new((char*) matB -> mtype -> id, matB -> dim, false);


        /* Diagonal matrix */

        if (matB -> mtype == &_matTypeD)
          {
            for (size_t i = 0; i < dimB[0]; i++)
              {
                /* Value at location in matA's and matX's matrices */
                size_t loc[2] = {i, i};
                double valA = mat_get_value(matA, loc);
                double valX = mat_get_value(matB, loc);

                /* Matrix is not invertible */
                if (valA == 0.)
                  {
                    matX = mat_free(matX);

                    return NULL;
                  }

                /* Insert the new value */
                mat_set_value(matX, loc, valX / valA);
              }

            return matX;
          }


        /* Other matrix type */

        #pragma omp parallel for schedule(dynamic) shared(matX, matA, matB)
        for (size_t i = 0; i < dimB[0]; i++)
          {
            /* Matrix not invertible (freed matX) */
            if (matX == NULL)
                continue;

            /* Value at location in matA's matrix */
            size_t locA[2] = {i, i};
            double valA = mat_get_value(matA, locA);

            /* Matrix is not invertible */
            if (valA == 0.)
              {
                matX = mat_free(matX);

                continue;
              }

            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matB, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Value at location in matB's matrix */
                size_t locX[2] = {i, j};
                double valX = mat_get_value(matB, locX);

                /* Insert the new value */
                mat_set_value(matX, locX, valX / valA);
              }
          }

        return matX;
      }


    /* Other type */

    /* Create general matrix with matB's struct */
    mat_t *matX = mat_cp_struct(matB, _matTypeF.id, 0); // 0 : Set all ids to "f"
    mat_set_matrix(matX, matB);

    #pragma omp parallel shared(matX, matA)
    for (int iInt = (int) dimB[0] - 1; iInt > -1; iInt--)
      {
        size_t i = (size_t) iInt;

        #pragma omp for schedule(dynamic)
        for (size_t j = 0; j < dimB[1]; j++)
          {
            /* Matrix not invertible (freed matX) */
            if (matX == NULL)
                continue;

            /* Value at location in matX' matrix */
            size_t locX[2] = {i, j};
            double valX = mat_get_value(matX, locX);

            /* Back substitution */
            for (size_t k = i + 1; k < dimB[0]; k++)
              {
                /* Value at location in matA's matrix */
                size_t locA_[2] = {i, k};
                double valA_ = mat_get_value(matA, locA_);

                /* Value at location in matX' matrix */
                size_t locX_[2] = {k, j};
                double valX_ = mat_get_value(matX, locX_);

                /* Adjust valX */
                valX -= valA_ * valX_;
              }

            /* Value at location in matA's matrix */
            size_t locA[2] = {i, i};
            double valA = mat_get_value(matA, locA);

            /* Matrix is not invertible */
            if (valA == 0.)
              {
                matX = mat_free(matX);

                continue;
              }

            /* Set the value in matX */
            mat_set_value(matX, locX, valX / valA);
          }
      }

    return matX;
}


static mat_t *_mat_backsub_b(mat_t *matA, mat_t *matB)
{
    /*

        Solve the (block) linear equation

            A X = B

        where A is an upper triangular matrix with back substitution.

        This function assumes A and B to be block matrices.

    */

    /* Dimensions */
//    size_t dimA[2] = {matA -> dim[0], matA -> dim[1]};
    size_t dimB[2] = {matB -> dim[0], matB -> dim[1]};


    /* Diagonal matrix */

    if (matA -> mtype == &_matTypeD)
      {
        /* Only need to copy matB's struct (not the elements) */
        mat_t *matX = mat_new((char*) matB -> mtype -> id, matB -> dim, true);

        /* Keep track of visited indices */
        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        for (size_t i = 0; i < dimB[0]; i++)
          {
            /* Block at location in matA's matrix */
            size_t locA[2] = {i, i};
            mat_t *matABlock = mat_fget_block(matA, locA);

            /* Matrix is not invertible */
            if (matABlock == NULL)
              {
                matX = mat_free(matX);

                return NULL;
              }

            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(matB, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Index at location in matB's memory */
                size_t locX[2] = {i, j};
                size_t indexMem = 0;
                int success = mat_mget_index(matX, locX, &indexMem);

                /* Block does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) indexMem, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                /* Get the old block */
                mat_t *matXBlock = mat_fget_block(matB, locX);

                /* Get the new block */
                mat_t *matXBlockNew = mat_backsub(matABlock, matXBlock, NULL);

                /* Block is not invertible */
                if (matXBlockNew == NULL)
                  {
                    matX = mat_free(matX);

                    return NULL;
                  }

                /* Insert the new block */
                mat_fset_block(matX, locX, matXBlockNew);

                /* Free memory */
                matXBlock = mat_ffree(matXBlock);
              }

            /* Free memory */
            matABlock = mat_ffree(matABlock);
          }

        /* Free memory */
        free(visitedIndices);

        return matX;
      }


    /* Other type */

    /* Create general matrix with matB's block struct */
    mat_t *matX = mat_cp_struct(matB, _matTypeF.id, 0); // 0 : Set all ids to "f"
    mat_set_matrix(matX, matB);

    for (size_t j = 0; j < dimB[1]; j++)
      {
        for (int iInt = (int) dimB[0] - 1; iInt > -1; iInt--)
          {
            size_t i = (size_t) iInt;

            /* Block at location in matX' matrix */
            size_t locX[2] = {i, j};
            mat_t *matXBlock = mat_fget_block(matX, locX);

            /* Skip if a block does not exist */
            if (matXBlock == NULL)
                continue;

            /* Back substitution */
            for (size_t k = i + 1; k < dimB[0]; k++)
              {
                /* Block at location in matA's matrix */
                size_t locA_[2] = {i, k};
                mat_t *matABlock_ = mat_fget_block(matA, locA_);

                /* Block at location in matX' matrix */
                size_t locX_[2] = {k, j};
                mat_t *matXBlock_ = mat_fget_block(matX, locX_);

                /* Block does not exist */
                if (matABlock_ == NULL || matXBlock_ == NULL)
                  {
                    matABlock_ = mat_ffree(matABlock_);
                    matXBlock_ = mat_ffree(matXBlock_);

                    continue;
                  }

                /* Get the new block (generelisation of the ordinary case) */
                matXBlock = mat_mult_mat_ip('s', matABlock_, matXBlock_, matXBlock, NULL);

                /* Free memory */
                matABlock_ = mat_ffree(matABlock_);
                matXBlock_ = mat_ffree(matXBlock_);
              }

            /* Block at location in matA's matrix */
            size_t locA[2] = {i, i};
            mat_t *matABlock = mat_fget_block(matA, locA);

            /* Solve A' X' = B' by back substitution */
            mat_t *matXBlock_ = mat_backsub(matABlock, matXBlock, NULL);

            /* Matrix is not invertible */
            if (matXBlock_ == NULL)
              {
                matX = mat_free(matX);

                return NULL;
              }

            /* Swap matrices */
            mat_swap(matXBlock, matXBlock_);

            /* Free memory */
            matABlock = mat_ffree(matABlock);
            matXBlock_ = mat_free(matXBlock_);
            mat_ffree(mat_mget_block(matX, locX)); // Need to free everything but matX -> matrix (not been freed previously)

            /* Reinsert the block into matX */
            mat_fset_block(matX, locX, matXBlock);
          }
      }

    return matX;
}


static mat_t *_mat_backsub_g(mat_t *matA, mat_t *matB)
{
    /*

        Solve the (general) linear equation

            A X = B

        where A is an upper triangular matrix and B is a general matrix.

        This function does not assume A and B to have any particular block structure.

    */

    /* Dimensions */
//    size_t dimA[2] = {0, 0};
//    mat_get_fdim(matA, dimA);

    size_t dimB[2] = {0, 0};
    mat_get_fdim(matB, dimB);


    /* Initialise matX as the full version of matB */
    mat_t *matX = mat_cp_struct(matB, _matTypeF.id, 0); // 0 : Set all ids to "f"
    mat_set_matrix(matX, matB);

    /* Perform the back substitution */
    #pragma omp parallel shared(matX, matA)
      {
        for (int iInt = (int) dimB[0] - 1; iInt > -1; iInt--)
          {
            size_t i = (size_t) iInt;

            #pragma omp for schedule(dynamic)
            for (size_t j = 0; j < dimB[1]; j++)
              {
                /* Matrix not invertible (freed matX) */
                if (matX == NULL)
                    continue;

                /* Value at location in matX' matrix */
                size_t locX[2] = {i, j};
                double valX = mat_get_value(matX, locX);

                /* Back substitution */
                for (size_t k = i + 1; k < dimB[0]; k++)
                  {
                    /* Value at location in matA's matrix */
                    size_t locA_[2] = {i, k};
                    double valA_ = mat_get_value(matA, locA_);

                    /* Value at location in matX' matrix */
                    size_t locX_[2] = {k, j};
                    double valX_ = mat_get_value(matX, locX_);

                    /* Adjust valX */
                    valX -= valA_ * valX_;
                  }

                /* Value at location in matA's matrix */
                size_t locA[2] = {i, i};
                double valA = mat_get_value(matA, locA);

                /* Matrix is not invertible */
                if (valA == 0.)
                  {
                    matX = mat_free(matX);

                    continue;
                  }

                /* Set the value in matX */
                mat_set_value(matX, locX, valX / valA);
              }
          }
      }

    return matX;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Matrix Inversion   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_inv(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Invert a matrix (mat) by first performing a QR decomposing of the matrix

    */

    /* Get the QR decomposition */
    mats_t *matsQR = mat_qr(mat, reduce);

    /* Q and R matrices */
    mat_t *matQ = matsQR -> mat[0];
    mat_t *matR = matsQR -> mat[1];

    /* Q^T matrix */
    matQ = mat_ftp(matQ, NULL);

    /* Solve R A^-1 = Q^T using back substitution (allow for reduction here) */
    mat_t *matInv = mat_backsub(matR, matQ, reduce);

    /* Free memory */
    matsQR = mats_free_full(matsQR);

    return matInv;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_inv_lapack_o(mat_t *mat);
static mat_t *_mat_inv_lapack_b(mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_inv_lapack(mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Invert a matrix (mat) by using LAPACK routines

    */

    /* Inverse matrix */
    mat_t *matInv = NULL;

    /* Ordinary matrix */
    if (mat -> mblock == NULL)
      {
        matInv = _mat_inv_lapack_o(mat);
      }

    /* Block matrix */
    else
      {
        matInv = _mat_inv_lapack_b(mat);
      }

    /* Reduce the matrix */
    matInv = mat_reduce(matInv, reduce);

    return matInv;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_inv_lapack_o(mat_t *mat)
{
    /*

        Invert an ordinary matrix (mat) by using LAPACK routines

    */

    /* Inverse matrix */
    mat_t *matInv = NULL;


    /* LAPACK VARIABLES */

    /* Dimensions of matrix */
    int laM = (int) mat -> dim[0];
    int laN = (int) mat -> dim[1];


    /* Leading dimension of A */
    int laLDA = laN;


    /* Pivot matrix */
    int *laIPIV = malloc(sizeof(int) * ((size_t) laN));


    /* Size of WORK (should be larger than 4*N) */
    int laLWORK = 6 * laN;

    /* WORK */
    double *laWORK = malloc(sizeof(double) * ((size_t) laLWORK));


    /* INFO == 0 (success), == -i (i-th argument had illegal value, == i (algorithm failed...)) */
    int laINFO;


    /* Diagonal Matrix */
    if (mat -> mtype == &_matTypeD)
      {
        matInv = mat_inv(mat, NULL);
      }

    /* Symmetric Matrix */
    else if (mat -> mtype == &_matTypeS)
      {
        /* Inverse matrix */
        matInv = mat_trafo(mat, _matTypeF.id, mat -> dim, false, NULL);

        /* UDU decomposition (Should form upper triangular part of inverse but actually forms lower triangular part (?)) */
        char laUPLO = 'U';

        /* Input + Output matrix (LDA, N) */
        double *laA = matInv -> matrix;

        /* Calculate the LDL decomposition */
        dsytrf_(&laUPLO, &laN, laA, &laLDA, laIPIV, laWORK, &laLWORK, &laINFO);

        if (laINFO != 0)
          {
            printf("UDU decomposition failed.\n");
            printf("LAPACK output : %d\n", laINFO);

            matInv = mat_free(matInv);

            goto freeMemory;
          }

        /* Calculate the inverse */
        dsytri_(&laUPLO, &laN, laA, &laLDA, laIPIV, laWORK, &laLWORK, &laINFO);

        if (laINFO != 0)
          {
            printf("Inversion failed.\n");
            printf("LAPACK output : %d\n", laINFO);

            matInv = mat_free(matInv);

            goto freeMemory;
          }

        /* Trafo to symmetric matrix */
        /* TODO: Let mat_trafo be able to handle this or create new function... */
        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            for (size_t j = i + 1; j < mat -> dim[1]; j++)
              {
                size_t loc[2] = {i, j};
                size_t locTp[2] = {j, i};

                mat_set_value(matInv, loc, mat_get_value(matInv, locTp));
              }
          }

        matInv = mat_trafo_ip(matInv, "s", mat -> dim, false, NULL);
      }

    /* Other Matrix Type */
    else
      {
        /* Inverse matrix */
        matInv = mat_trafo(mat, _matTypeF.id, mat -> dim, false, NULL);


        /* Input + Output matrix (LDA, N) */
        double *laA = matInv -> matrix;

        /* Calculate the LU decomposition */
        dgetrf_(&laM, &laN, laA, &laLDA, laIPIV, &laINFO);

        if (laINFO != 0)
          {
            printf("LU decomposition failed.\n");
            printf("LAPACK output : %d\n", laINFO);

            matInv = mat_free(matInv);

            goto freeMemory;
          }

        /* Calculate the inverse */
        dgetri_(&laN, laA, &laLDA, laIPIV, laWORK, &laLWORK, &laINFO);

        if (laINFO != 0)
          {
            printf("Inversion failed.\n");
            printf("LAPACK output : %d\n", laINFO);

            matInv = mat_free(matInv);

            goto freeMemory;
          }
      }


    freeMemory:

    /* Free memory */
    free(laIPIV);
    free(laWORK);

    return matInv;
}


static mat_t *_mat_inv_lapack_b(mat_t *mat)
{
    /*

        Invert a block matrix (mat) by using LAPACK routines

    */

    /* Diagonal matrix */
    if (mat -> mtype == &_matTypeD)
      {
        mat_t *matInv = mat_new(_matTypeD.id, mat -> dim, true);

        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            /* Block at location */
            size_t loc[2] = {i, i};
            mat_t *matBlock = mat_fget_block(mat, loc);

            /* Inverse block */
            mat_t *matInvBlock = mat_inv_lapack(matBlock, NULL);

            /* Set the inverse block */
            mat_fset_block(matInv, loc, matInvBlock);

            /* Free memory */
            matBlock = mat_ffree(matBlock);
          }

        return matInv;
      }


    /* 2x2 Block Matrix */
    if (mat -> dim[0] == 2 && mat -> dim[1] == 2)
      {
        mat_t *matInv = mat_new(_matTypeF.id, mat -> dim, true);

        /* D^-1 Matrix */
        size_t locD[2] = {1, 1};
        mat_t *matD = mat_fget_block(mat, locD);
        mat_t *matDInv = mat_inv_lapack(matD, NULL);

        if (matDInv == NULL)
          {
            matD = mat_ffree(matD);

            goto fullInv;
          }

        /* -D^-1 C Matrix */
        size_t locC[2] = {1, 0};
        mat_t *matC = mat_fget_block(mat, locC);
        mat_t *matDInvC = mat_mult(matDInv, matC, -1., NULL);

        /* -B D^-1 Matrix */
        size_t locB[2] = {0, 1};
        mat_t *matB = mat_fget_block(mat, locB);
        mat_t *matBDInv = mat_mult(matB, matDInv, -1., NULL);

        /* (A - B D^-1 C)^-1 Matrix */
        size_t locA[2] = {0, 0};
        mat_t *matA = mat_fget_block(mat, locA);
        mat_t *matAInv = mat_mult_mat(matBDInv, matC, NULL);

        if (!mat_comp_bstruct(matA, matAInv) || matA -> mtype != matAInv -> mtype)
          {
            mat_t *matAInv_ = mat_cp_struct(matA, _matTypeF.id, 0);

            mat_set_matrix(matAInv_, matAInv);
            mat_swap(matAInv_, matAInv);

            matAInv_ = mat_free(matAInv_);
          }

        matAInv = mat_add_ip('a', matA, matAInv, matAInv, NULL); // A - B D^-1 C

//        matAInv = mat_inv_lapack_ip(matAInv, NULL); // TODO
        mat_t *matAInv_ = mat_inv_lapack(matAInv, NULL); // (A - B D^-1 C)^-1

        mat_swap(matAInv, matAInv_);
        matAInv_ = mat_free(matAInv_);

        if (matAInv == NULL)
          {
            matA = mat_ffree(matA);
            matB = mat_ffree(matB);
            matC = mat_ffree(matC);
            matD = mat_ffree(matD);

            matDInvC = mat_free(matDInvC);
            matBDInv = mat_free(matBDInv);

            matDInv = mat_free(matDInv);

            goto fullInv;
          }

        mat_fset_block(matInv, locA, matAInv);

        /* -(A - B D^-1 C)^-1 B D^-1 Matrix */
        mat_t *matBInv = mat_mult_mat(matAInv, matBDInv, NULL);
        mat_fset_block(matInv, locB, matBInv);

        /* -D^-1 C (A - B D^-1 C)^-1 Matrix */
        mat_t *matCInv = mat_mult_mat(matDInvC, matAInv, NULL);
        mat_fset_block(matInv, locC, matCInv);

        matDInvC = mat_free(matDInvC);

        /* D^-1 + D^-1 C (A - B D^-1 C)^-1 B D^-1 Matrix */
        matDInvC = mat_mult_mat(matCInv, matBDInv, NULL); // D^-1 C (A - B D^-1 C)^-1 B D^-1

        if (matDInvC -> mtype != &_matTypeNull && matDInv -> mtype != matDInvC -> mtype)
          {
            mat_t *matDInv_ = mat_cp_struct(matDInv, _matTypeF.id, 0);

            mat_set_matrix(matDInv_, matDInv);
            mat_swap(matDInv_, matDInv);

            matDInv_ = mat_free(matDInv_);
          }

        matDInv = mat_add_ip('a', matDInv, matDInvC, matDInv, NULL); // D^-1 + D^-1 C (A - B D^-1 C)^-1 B D^-1

        matDInvC = mat_free(matDInvC);
        matBDInv = mat_free(matBDInv);

        mat_fset_block(matInv, locD, matDInv);


        /* Free memory */
        matA = mat_ffree(matA);
        matB = mat_ffree(matB);
        matC = mat_ffree(matC);
        matD = mat_ffree(matD);

        return matInv;
      }


    fullInv: ;

    /* Other type */

    /* Full dimensions */
    size_t fdim[2];
    mat_get_fdim(mat, fdim);

    /* Full matrix */
    mat_t *matFull = mat_trafo(mat, _matTypeF.id, fdim, false, NULL);
    mat_t *matInv_ = mat_inv_lapack(matFull, NULL);
    matFull = mat_free(matFull);

    if (matInv_ == NULL)
      {
        return NULL;
      }

    /* Copy the struct */
    mat_t *matInv = mat_cp_struct(mat, _matTypeF.id, 0); // 0 : Set all ids to "f"
    mat_set_matrix(matInv, matInv_);
    matInv_ = mat_free(matInv_);

    return matInv;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mat_t *_mat_inv_cov_o(mat_t *inv(mat_t*, mat_reduce_t*), mat_t *mat);
static mat_t *_mat_inv_cov_b(mat_t *inv(mat_t*, mat_reduce_t*), mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mat_inv_cov(mat_t *inv(mat_t*, mat_reduce_t*), mat_t *mat, mat_reduce_t *reduce)
{
    /*

        Invert a covariance matrix

    */

    /* Inverse Matrix */
    mat_t *matInv = NULL;


    /* Invert the Matrix */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        matInv = _mat_inv_cov_o(inv, mat);
      }

    /* Block Matrix */
    else
      {
        matInv = _mat_inv_cov_b(inv, mat);
      }

    /* Attempt to reduce the Matrix */
    matInv = mat_reduce(matInv, reduce);

    return matInv;
}


/*  ------------------------------------------------------------------------------------------------------  */


static mat_t *_mat_inv_cov_o(mat_t *inv(mat_t*, mat_reduce_t*), mat_t *mat)
{
    /*

        Invert an ordinary covariance matrix

    */

    /* Inverse Matrix */
    mat_t *matInv = NULL;

    /* Diagonal Matrix */
    if (mat -> mtype == &_matTypeD)
      {
        matInv = inv(mat, NULL);
      }

    /* Other Matrix Type */
    else
      {
        /* Correlation matrix */
        mat_t *matCorr = mat_corr(mat, NULL);

        /* Inverse of Correlation Matrix */
        matInv = inv(matCorr, NULL);

        /* Free memory */
        matCorr = mat_free(matCorr);

        /* Inverse of Covariance Matrix */
        #pragma omp parallel for
        for (size_t i = 0; i < matInv -> dim[2]; i++)
          {
            /* Location and value in matInv's matrix */
            size_t loc[2];
            mat_get_loc(matInv, i, loc);

            double val = mat_get_value(matInv, loc);

            /* Location and values mat's matrix */
            size_t locR[2] = {loc[0], loc[0]};
            double valR = sqrt(mat_get_value(mat, locR));

            size_t locC[2] = {loc[1], loc[1]};
            double valC = sqrt(mat_get_value(mat, locC));

            mat_set_value(matInv, loc, val / (valR * valC));
          }
      }

    return matInv;
}


static mat_t *_mat_inv_cov_b(mat_t *inv(mat_t*, mat_reduce_t*), mat_t *mat)
{
    /*

        Invert a block covariance matrix

    */

    /* Inverse Matrix */
    mat_t *matInv = NULL;

    /* Diagonal Matrix */
    if (mat -> mtype == &_matTypeD)
      {
        matInv = mat_new(_matTypeD.id, mat -> dim, true);

        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            size_t loc[2] = {i, i};
            mat_t *matBlock = mat_fget_block(mat, loc);

            mat_t *matInvBlock = mat_inv_cov(inv, matBlock, NULL);
            mat_fset_block(matInv, loc, matInvBlock);

            matBlock = mat_ffree(matBlock);
          }
      }

    /* Other Matrix Type */
    else
      {
        /* Correlation matrix */
        mat_t *matCorr = mat_corr(mat, NULL);

        /* Inverse of Correlation Matrix */
        matInv = inv(matCorr, NULL);

        /* Free memory */
        matCorr = mat_free(matCorr);

        /* (Inverse square root of) Diagonal matrix */
        mat_t *matDiag = mat_diag_func(mat_diag(mat), mat_diag_func_invsqrt);

        /* Get the correlation block */
        matInv = mat_mult_mat_ip('0', matDiag, matInv, matInv, NULL);
        matInv = mat_mult_mat_ip('0', matInv, matDiag, matInv, NULL);

        /* Free memory */
        matDiag = mat_free(matDiag);

//        /* Full dimensions */
//        size_t fdim[2] = {0, 0};
//        mat_get_fdim(matInv, fdim);
//
//        printf("\t%ld %ld - %ld %ld - %d\n", matInv -> dim[0], matInv -> dim[1], fdim[0], fdim[1], matInv -> mblock == NULL);
//
//        /* Inverse of Covariance Matrix */
//        /* TODO: This will be slow for large matrices (10,000 x 10,000) -> improve this */
//        size_t visitedIndicesSize = 0;
//        double *visitedIndices = NULL;
//        bool successInsert;
//
////        #pragma omp parallel for
//        for (size_t i = 0; i < fdim[0]; i++)
//          {
//            size_t locR[2] = {i, i};
//            double valR = sqrt(mat_get_value(mat, locR));
//
//            for (size_t j = 0; j < fdim[1]; j++)
//              {
//                /* Index of location in mat's memory */
//                size_t loc[2] = {i, j};
//                size_t index = 0;
//                int success = mat_mget_index(matInv, loc, &index);
//
//                /* Element does not exist */
//                if (success == -1)
//                    continue;
//
//                /* Search for and insert the index (as double) in visitedIndices */
//                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);
//
//                /* Skip if an index has already been visited */
//                if (!successInsert)
//                    continue;
//
//                /* Location and value in mat's matrix */
//                size_t locC[2] = {j, j};
//                double valC = sqrt(mat_get_value(mat, locC));
//
//                /* Value in matInv's matrix */
//                double val = mat_get_value(matInv, loc);
//
//                mat_set_value(matInv, loc, val / (valR * valC));
//              }
//          }
//
//        /* Free memory */
//        free(visitedIndices);
      }

    return matInv;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     OUTPUT MATRICES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------------   Output Matrix   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static char *_mats_output_fprintf_get_header(mats_t *mats, size_t offset);
static char *_mats_output_fwrite_get_header(mats_t *mats, size_t offset);

static int _mats_output_fprintf_header(FILE *stream, char *header);
static int _mats_output_fwrite_header(FILE *stream, char *header);

static int _mats_output_fprintf(FILE *stream, mat_t *mat, int precision);
static int _mats_output_fwrite(FILE *stream, mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


int mats_output(char *fileName, mats_t *mats, int *precision)
{
    /*

        Output a matrix (mat) to a file (fileName) with an optional additional header

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

    /* Write the data to binary file */
    if (precision == NULL)
      {
        header = _mats_output_fwrite_get_header(mats, 0);

        _mats_output_fwrite_header(stream, header);

        for (size_t i = 0; i < mats -> size; i++)
          {
            _mats_output_fwrite(stream, mats -> mat[i]);
          }
      }

    /* Write the data to text file */
    else
      {
        header = _mats_output_fprintf_get_header(mats, 0);

        _mats_output_fprintf_header(stream, header);

        for (size_t i = 0; i < mats -> size; i++)
          {
            _mats_output_fprintf(stream, mats -> mat[i], *precision);

            /* Add a new line to separate the matrices */
            if (i < mats -> size - 1)
              {
                fprintf(stream, "\n\n");
              }
          }
      }

    /* Free memory */
    free(header);

    /* Close the file */
    fclose(stream);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fprintf_header(FILE *stream, char *header)
{
    /*

        Write the total header to stream in text mode

    */

    /* Write header to file */
    fprintf(stream, "%s", header);

    return 0;
}


static int _mats_output_fwrite_header(FILE *stream, char *header)
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


static int _mats_output_fprintf(FILE *stream, mat_t *mat, int precision)
{
    /*

        Write a mat to stream in text mode

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        bool newline;

        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            newline = true;

            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(mat, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Index of location in mat's memory */
                size_t loc[2] = {i, j};
                size_t index = 0;
                int success = mat_mget_index(mat, loc, &index);

                /* Element does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                /* Add one white space char if newline is false */
                if (!newline)
                    fprintf(stream, " ");

                else
                    newline = false;

                /* Write the element to file */
                fprintf(stream, "%.*e", precision, ((double*) mat -> matrix)[index]);
              }

            /* Line break */
            fprintf(stream, "\n");
          }

        /* Free memory */
        free(visitedIndices);
      }

    /* Block Matrix */
    else
      {
        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(mat, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Index at location in mat's memory */
                size_t loc[2] = {i, j};
                size_t index = 0;
                int success = mat_mget_index(mat, loc, &index);

                /* Block does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                /* Block at location in mat's matrix */
                mat_t *matBlock = mat_fget_block(mat, loc);

                /* Write the block to file */
                _mats_output_fprintf(stream, matBlock, precision);

                /* Line break */
                if (i < mat -> dim[0] - 1 || j < mat -> dim[1] - 1)
                    fprintf(stream, "\n");

                /* Free memory */
                matBlock = mat_ffree(matBlock);
              }
          }

        /* Free memory */
        free(visitedIndices);
      }

    return 0;
}


static int _mats_output_fwrite(FILE *stream, mat_t *mat)
{
    /*

        Write a mat to stream in binary mode

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* No location trafos */
        if (mat -> mltrfs -> size == 0)
            fwrite(mat -> matrix, sizeof(double), mat -> dim[2], stream);

        /* TODO: This is slow for large matrices (~ 10,000 x 10,000) -> undo location trafo and output as above */
        else
          {
            size_t visitedIndicesSize = 0;
            double *visitedIndices = NULL;
            bool successInsert;

            for (size_t i = 0; i < mat -> dim[0]; i++)
              {
                /* Column bounds */
                size_t cBounds[2] = {0, 0};
                mat_get_cbounds(mat, i, cBounds);

                for (size_t j = cBounds[0]; j < cBounds[1]; j++)
                  {
                    /* Index of location in mat's memory */
                    size_t loc[2] = {i, j};
                    size_t index = 0;
                    int success = mat_mget_index(mat, loc, &index);

                    /* Element does not exist */
                    if (success == -1)
                        continue;

                    /* Search for and insert the index (as double) in visitedIndices */
                    misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                    /* Skip if an index has already been visited */
                    if (!successInsert)
                        continue;

                    /* Write the element to file */
                    double val = ((double*) mat -> matrix)[index];
                    fwrite(&val, sizeof(double), 1, stream);
                  }
              }

            /* Free memory */
            free(visitedIndices);
          }
      }

    /* Block Matrix */
    else
      {
        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            /* Column bounds */
            size_t cBounds[2] = {0, 0};
            mat_get_cbounds(mat, i, cBounds);

            for (size_t j = cBounds[0]; j < cBounds[1]; j++)
              {
                /* Index of location in mat's memory */
                size_t loc[2] = {i, j};
                size_t index = 0;
                int success = mat_mget_index(mat, loc, &index);

                /* Block does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                /* Block at location in mat's matrix */
                mat_t *matBlock = mat_fget_block(mat, loc);

                /* Write the block to file */
                _mats_output_fwrite(stream, matBlock);

                /* Free memory */
                matBlock = mat_ffree(matBlock);
              }
          }

        /* Free memory */
        free(visitedIndices);
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------   Construct a Header for a Text File   --------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mats_output_fprintf_get_header_guide(char *header);
static int _mats_output_fprintf_get_header_section(mat_t *mat, size_t depth, size_t index, size_t *matline, char *header);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static char *_mats_output_fprintf_get_header(mats_t *mats, size_t offset)
{
    /*

        Get the header for a matrix (mat)

    */

    char *header = malloc(HEADER_LEN);
    header[0] = '\0';

    /* Add guide about displayed information */
    _mats_output_fprintf_get_header_guide(header);

    /* Add a new line */
    strcat(header, "\n");

    /* Line of the first matrix */
    size_t matline = offset + 1;

    for (size_t i = 0; i < mats -> size; i++)
      {
        _mats_output_fprintf_get_header_section(mats -> mat[i], 0, i, &matline, header);

        /* Add a line break */
        strcat(header, "\n");
        matline++;
      }

    /* Separate the header from the data */
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
        char *headerNew = _mats_output_fprintf_get_header(mats, count);

        free(header);
        header = headerNew;
      }

    return header;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fprintf_get_header_guide(char *header)
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

    snprintf(line, HEADER_LINE_LEN, "* Information about the matrices is displayed as:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     matrix at line <line> : <label>, <type>, (<row dim>, <column dim>), (<depth in block>, <index at depth>)\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     block matrix at line <line> : <label>, <type>, (<row dim>, <column dim>), (<depth in block>, <index at depth>)\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*       matrix at line <line> : <label>, <type>, (<row dim>, <column dim>), (<depth in block>, <index at depth>)\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "* Can have types:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - f : full\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - d : diagonal\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - s : symmetric\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - ut : upper-triangular\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - lt : lower-triangular\n");
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
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mats_output_fprintf_get_header_line(mat_t *mat, size_t depth, size_t index, size_t *matline, char *header);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fprintf_get_header_section(mat_t *mat, size_t depth, size_t index, size_t *matline, char *header)
{
    /*

        Get the header section for a matrix (mat) located at some depth and some index.

    */

    /* Get the header line for mat */
    _mats_output_fprintf_get_header_line(mat, depth, index, matline, header);

    /* Get the header lines for sub-matrices */
    if (mat -> mblock != NULL)
      {
        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            _mats_output_fprintf_get_header_section(((mat_t**) mat -> matrix)[i], depth + 1, i, matline, header);
          }
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fprintf_get_header_line(mat_t *mat, size_t depth, size_t index, size_t *matline, char *header)
{
    /*

        Get the header line for a matrix (mat) located at some depth and some index.

    */

    char *headerLine = malloc(HEADER_LINE_LEN);

    /* Without label */
    if (mat -> label == NULL)
      {
        if (mat -> mblock != NULL)
            snprintf(headerLine, HEADER_LINE_LEN, "%*sblock matrix at line %ld : --, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matline, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);

        else
            snprintf(headerLine, HEADER_LINE_LEN, "%*smatrix at line %ld : --, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matline, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);
      }

    /* With label */
    else
      {
        if (mat -> mblock != NULL)
            snprintf(headerLine, HEADER_LINE_LEN, "%*sblock matrix at line %ld : %s, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matline, mat -> label, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);

        else
            snprintf(headerLine, HEADER_LINE_LEN, "%*smatrix at line %ld : %s, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matline, mat -> label, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);
      }

    /* Increment matline for ordinary matrices (+1 since matrices are always separated by \n) */
    if (mat -> mblock == NULL)
        *matline += mat -> dim[0] + 1;

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

static int _mats_output_fwrite_get_header_guide(char *header);
static int _mats_output_fwrite_get_header_section(mat_t *mat, size_t depth, size_t index, size_t *matbyte, char *header);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static char *_mats_output_fwrite_get_header(mats_t *mats, size_t offset)
{
    /*

        Get the header for a matrix (mat)

    */

    char *header = malloc(HEADER_LEN);
    header[0] = '\0';

    /* Add guide about displayed information */
    _mats_output_fwrite_get_header_guide(header);

    /* Add a new line */
    strcat(header, "\n");

    /* Line of the first matrix */
    size_t matbyte = offset + 1;

    for (size_t i = 0; i < mats -> size; i++)
      {
        _mats_output_fwrite_get_header_section(mats -> mat[i], 0, i, &matbyte, header);

        /* Add a line break */
        strcat(header, "\n");
      }

    /* Separate the header from the data */
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
        char *headerNew = _mats_output_fwrite_get_header(mats, strlen(header));

        free(header);
        header = headerNew;
      }

    return header;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fwrite_get_header_guide(char *header)
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

    snprintf(line, HEADER_LINE_LEN, "* Information about the matrices is displayed as:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     matrix at byte <byte> : <label>, <type>, (<row dim>, <column dim>), (<depth in block>, <index at depth>)\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     block matrix at byte <byte> : <label>, <type>, (<row dim>, <column dim>), (<depth in block>, <index at depth>)\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*       matrix at byte <byte> : <label>, <type>, (<row dim>, <column dim>), (<depth in block>, <index at depth>)\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*                 :\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "* Can have types:\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - f : full\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - d : diagonal\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - s : symmetric\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - ut : upper-triangular\n");
    strcat(header, line);

    snprintf(line, HEADER_LINE_LEN, "*     - lt : lower-triangular\n");
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
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mats_output_fwrite_get_header_line(mat_t *mat, size_t depth, size_t index, size_t *matbyte, char *header);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fwrite_get_header_section(mat_t *mat, size_t depth, size_t index, size_t *matbyte, char *header)
{
    /*

        Get the header section for a matrix (mat) located at some depth and some index.

    */

    /* Get the header line for mat */
    _mats_output_fwrite_get_header_line(mat, depth, index, matbyte, header);

    /* Get the header lines for sub-matrices */
    if (mat -> mblock != NULL)
      {
        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            _mats_output_fwrite_get_header_section(((mat_t**) mat -> matrix)[i], depth + 1, i, matbyte, header);
          }
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_output_fwrite_get_header_line(mat_t *mat, size_t depth, size_t index, size_t *matbyte, char *header)
{
    /*

        Get the header line for a matrix (mat) located at some depth and some index.

    */

    char *headerLine = malloc(HEADER_LINE_LEN);

    /* Without label */
    if (mat -> label == NULL)
      {
        if (mat -> mblock != NULL)
            snprintf(headerLine, HEADER_LINE_LEN, "%*sblock matrix at byte %ld : --, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matbyte, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);

        else
            snprintf(headerLine, HEADER_LINE_LEN, "%*smatrix at byte %ld : --, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matbyte, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);
      }

    /* With label */
    else
      {
        if (mat -> mblock != NULL)
            snprintf(headerLine, HEADER_LINE_LEN, "%*sblock matrix at byte %ld : %s, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matbyte, mat -> label, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);

        else
            snprintf(headerLine, HEADER_LINE_LEN, "%*smatrix at byte %ld : %s, %s, (%ld, %ld), (%ld, %ld)\n",
                     (int) depth, "", *matbyte, mat -> label, mat -> mtype -> id, mat -> dim[0], mat -> dim[1], depth, index);
      }

    /* Increment matbyte for ordinary matrices */
    if (mat -> mblock == NULL)
        *matbyte += sizeof(double) * mat -> dim[2];

    /* Add the line to the header */
    strcat(header, headerLine);

    /* Free memory */
    free(headerLine);

    return 0;
}






/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INPUT MATRIX     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
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

    bool block;

    bool binary;
    char *loc;

    char *label;
    char *type;

    char **dim;

    char *depth;
    char *index;

} _mats_input_header_line_t;


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static _mats_input_header_line_t *_mats_input_get_header_line_free(_mats_input_header_line_t *headerLine)
{
    /*

        Free a _mats_input_header_line_t struct

    */

    /* If headerLine is NULL simply return NULL */
    if (headerLine == NULL)
        return NULL;

    /* Free headerLine's contents */
    free(headerLine -> loc);

    free(headerLine -> label);
    free(headerLine -> type);

    free(headerLine -> dim[0]);
    free(headerLine -> dim[1]);
    free(headerLine -> dim);

    free(headerLine -> depth);
    free(headerLine -> index);

    /* Free headerLine itself */
    free(headerLine);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Input Matrix from File   -------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static mats_t *_mats_input_get_header(char *fileName, char *label);

static int _mats_input_fscanf_section(char *fileName, mat_t *mat);
static int _mats_input_fread_section(char *fileName, mat_t *mat);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *mats_input(char *fileName, char *label, bool fill)
{
    /*

        Read matrices from file.

        If fill is false only get the matrices as null matrices (i.e. dimensions, labels...)

    */

    /* Get the matrix structs and locations in the file */
    mats_t *mats = _mats_input_get_header(fileName, label);

    /* No matrix in the file */
    if (mats == NULL)
      {
        printf("Could not find a matrix with label '%s' in the file '%s'.\n", label, fileName);
        exit(1);

        return NULL;
      }

    /* Fill the matrices */
    if (fill)
      {
        for (size_t i = 0; i < mats -> size; i++)
          {
            if (mats -> mat[i] -> mfile -> binary)
                _mats_input_fread_section(fileName, mats -> mat[i]);

            else
                _mats_input_fscanf_section(fileName, mats -> mat[i]);
          }
      }

    return mats;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_input_fscanf_section(char *fileName, mat_t *mat)
{
    /*

        Read the matrix from a text file

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Open the file */
        FILE *stream = fopen(fileName, "r");

        /* Skip to the correct line */
        misc_stream_skip(stream, mat -> mfile -> loc - 1);

        /* Must get the correct type */
        mat_t *matDum = mat_new(mat -> mfile -> type, mat -> dim, false);
        mat_set_label(matDum, mat -> label);
        matDum -> mfile = _mat_file_cp(mat -> mfile);

        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        for (size_t i = 0; i < matDum -> dim[0]; i++)
          {
            for (size_t j = 0; j < matDum -> dim[1]; j++)
              {
                /* Index of location in mat's memory */
                size_t loc[2] = {i, j};
                size_t index = 0;
                int success = mat_mget_index(matDum, loc, &index);

                /* Element does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                int numScan = fscanf(stream, "%lf", &((double*) matDum -> matrix)[index]);
                (void) numScan; // Ignore number of variables scanned and assume everything went smoothly :)
              }
          }

        /* Swap matrices */
        mat_swap(mat, matDum);

        /* Free memory */
        matDum = mat_free(matDum);
        free(visitedIndices);

        /* Close the file */
        fclose(stream);
      }

    /* Block Matrix */
    else
      {
        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            for (size_t j = 0; j < mat -> dim[1]; j++)
              {
                /* Index of location in mat's memory */
                size_t loc[2] = {i, j};
                size_t index = 0;
                int success = mat_mget_index(mat, loc, &index);

                /* Block does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                /* Block */
                mat_t *matBlock = mat_mget_block(mat, loc);

                /* Get the block */
                _mats_input_fscanf_section(fileName, matBlock);
              }
          }

        /* Free memory */
        free(visitedIndices);
      }

    return 0;
}


static int _mats_input_fread_section(char *fileName, mat_t *mat)
{
    /*

        Read the matrix from a binary file

    */

    /* Ordinary Matrix */
    if (mat -> mblock == NULL)
      {
        /* Open the file */
        FILE *stream = fopen(fileName, "r");

        /* Skip to the correct byte */
        fseek(stream, (long int) mat -> mfile -> loc, SEEK_SET);

        /* Must get the correct type */
        mat_t *matDum = mat_new(mat -> mfile -> type, mat -> dim, false);
        mat_set_label(matDum, mat -> label);
        matDum -> mfile = _mat_file_cp(mat -> mfile);

//        size_t visitedIndicesSize = 0;
//        double *visitedIndices = NULL;
//        bool successInsert;
//
//        for (size_t i = 0; i < matDum -> dim[0]; i++)
//          {
//            for (size_t j = 0; j < matDum -> dim[1]; j++)
//              {
//                /* Index of location in matDum's memory */
//                size_t loc[2] = {i, j};
//                size_t index = 0;
//                int success = mat_mget_index(matDum, loc, &index);
//
//                /* Element does not exist */
//                if (success == -1)
//                    continue;
//
//                /* Search for and insert the index (as double) in visitedIndices */
//                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);
//
//                /* Skip if an index has already been visited */
//                if (!successInsert)
//                    continue;
//
//                size_t numRead = fread(&((double*) matDum -> matrix)[index], sizeof(double), 1, stream);
//                (void) numRead; // Ignore number of variables read and assume everything went smoothly :)
//              }
//          }

        size_t numRead = fread(matDum -> matrix, sizeof(double), matDum -> dim[2], stream); // This way works, as opposed to the 'output' function, since index trafos are not assumed
        (void) numRead; // Ignore number of variables read and assume everything went smoothly :)

        /* Swap matrices */
        mat_swap(mat, matDum);

        /* Free memory */
        matDum = mat_free(matDum);
//        free(visitedIndices);

        /* Close the file */
        fclose(stream);
      }

    /* Block Matrix */
    else
      {
        size_t visitedIndicesSize = 0;
        double *visitedIndices = NULL;
        bool successInsert;

        for (size_t i = 0; i < mat -> dim[0]; i++)
          {
            for (size_t j = 0; j < mat -> dim[1]; j++)
              {
                /* Index of location in mat's memory */
                size_t loc[2] = {i, j};
                size_t index = 0;
                int success = mat_get_index(mat, loc, &index);

                /* Block does not exist */
                if (success == -1)
                    continue;

                /* Search for and insert the index (as double) in visitedIndices */
                misc_insort(&visitedIndices, &visitedIndicesSize, (double) index, __ABSTOL__, &successInsert);

                /* Skip if an index has already been visited */
                if (!successInsert)
                    continue;

                /* Block at location in mat's matrix */
                mat_t *matBlock = mat_mget_block(mat, loc);

                /* Get the block */
                _mats_input_fread_section(fileName, matBlock);
              }
          }

        /* Free memory */
        free(visitedIndices);
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Matrix Struct from Header   ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _mats_input_header_skip_guide(FILE *stream);

static bool _mats_input_get_header_section(FILE *stream, char *label, mats_t *mats);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static mats_t *_mats_input_get_header(char *fileName, char *label)
{
    /*

        Read a file (fileName) and return the matrices and locations of each block (label == NULL) or for only
        one label.

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
    _mats_input_header_skip_guide(stream);

    /* Read blocks */
    mats_t *mats = mats_new(0);

    bool headerEnd = false;

    do
      {
        headerEnd = _mats_input_get_header_section(stream, label, mats);
      }

    while (!headerEnd);

    /* Close the file */
    fclose(stream);

    /* Header ended and correct label was not found */
    if (mats -> size == 0)
      {
        mats = mats_free(mats);
      }

    return mats;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _mats_input_header_skip_guide(FILE *stream)
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
        printf("File ended before finding the beginning of the matrix informations.\n");
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

static _mats_input_header_line_t *_mats_input_get_header_line_dissect(FILE *stream, bool *end);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static bool _mats_input_get_header_section(FILE *stream, char *label, mats_t *mats)
{
    /*

        Get the matrix struct of the block and the location of the matrix in the file.

    */

    /* Header line */
    bool headerEnd = false;
    _mats_input_header_line_t *headerLine = _mats_input_get_header_line_dissect(stream, &headerEnd);

    /* Header has ended */
    if (headerEnd)
        return true;

    /* Empty line */
    if (headerLine == NULL)
        return false;

    /* Matrix has the correct label */
    bool correctLabel = true;

    /* label is provided but matrix does not have the same label */
    if (label != NULL && strcmp(headerLine -> label, label))
      {
        correctLabel = false;
      }

    /* Matrix */
    mat_t *mat = NULL;

    /* Dimensions */
    size_t dim[2] = {(size_t) atoi(headerLine -> dim[0]), (size_t) atoi(headerLine -> dim[1])};

    /* Ordinary Matrix */
    if (!headerLine -> block)
      {
        /* Empty matrix only if the labels match */
        if (correctLabel)
          {
            mat = mat_new(_matTypeNull.id, dim, false);
          }
      }

    /* Block Matrix */
    if (headerLine -> block)
      {
        mat = mat_new(headerLine -> type, dim, true);

        /* Current size of mats */
        size_t currentSize = mats -> size;

        for (size_t i = 0; i < mat -> dim[2]; i++)
          {
            /* (i, j) location in matrix */
            size_t loc[2];
            mat -> mtype -> loc(mat -> dim, i, loc);

            /* Block at index i */
            mat_t *matBlock;

            if (correctLabel)
              {
                /* Already found the matrix with the correct label -> must use NULL */
                _mats_input_get_header_section(stream, NULL, mats);

                /* Get the block (last entry of mats -> mat) - TODO: If some matrix is missing in block -> undefined behaviour! */
                matBlock = mats -> mat[mats -> size - 1];

                /* Insert the block */
                mat_fset_block(mat, loc, matBlock);

                /* Remove the last entry - TODO: Can have general function to remove i'th matrix in mats */
                mats -> mat = realloc(mats -> mat, sizeof(mat_t*) * (mats -> size - 1));
                mats -> size -= 1;
              }

            else
              {
                /* Still searching for the correct matrix */
                _mats_input_get_header_section(stream, label, mats);

                /* If mats -> size was increased (ie correct label was found) size has changed -> break out of loop */
                if (currentSize != mats -> size)
                  {
                    break;
                  }
              }
          }

        /* Must free mat for the wrong label */
        if (!correctLabel)
          {
            mat = mat_free(mat);
          }
      }

    /* Set label and insert the matrix into mats */
    if (correctLabel)
      {
        /* Set the label */
        mat_set_label(mat, headerLine -> label);

        /* Set the matrix file parameters */
        mat -> mfile = _mat_file_new();

        /* Matrix Type */
        mat -> mfile -> type = misc_scat(1, headerLine -> type);

        /* Matrix location */
        sscanf(headerLine -> loc, "%ld", &(mat -> mfile -> loc));

        /* Binary flag */
        mat -> mfile -> binary = headerLine -> binary;

        /* Insert the matrix into mats */
        mats_set_mat(mats, mats -> size, mat);
      }

    /* Free memory */
    headerLine = _mats_input_get_header_line_free(headerLine);

    return false;
}


/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static char *_mats_input_get_header_line(FILE *stream);

static bool _mats_input_get_header_line_empty(char *line);
static bool _mats_input_get_header_line_end(char *line);

static char *_mats_input_get_header_line_type(char *line);

static char *_mats_input_get_header_line_loc(char *line);

static char *_mats_input_get_header_line_mtype(char *line);
static char **_mats_input_get_header_line_dim(char *line);
static char **_mats_input_get_header_line_depth(char *line);

static char *_mats_input_get_header_line_label(char *line);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static _mats_input_header_line_t *_mats_input_get_header_line_dissect(FILE *stream, bool *end)
{
    /*

        Read a line from a mat header and store the relevant information in a _mats_input_header_line_t struct.

    */

    /* Read the next line */
    char *line = _mats_input_get_header_line(stream);

    /* Empty line (copy and pasting stuff in file automatically adds \r at the end of every line... interesting stuff :)) */
    if (_mats_input_get_header_line_empty(line))
      {
        free(line);

        return NULL;
      }

    /* Header has ended */
    if (_mats_input_get_header_line_end(line))
      {
        *end = true;
        free(line);

        return NULL;
      }

    /* headerLine */
    _mats_input_header_line_t *headerLine = malloc(sizeof(_mats_input_header_line_t));
    headerLine -> dim = malloc(sizeof(char*) * 2);

    /* Type */
    char *type = _mats_input_get_header_line_type(line);

    headerLine -> block = (type[0] == 'b') ? true : false;
    headerLine -> binary = (type[1] == 'b') ? true : false;

    /* Location */
    headerLine -> loc = _mats_input_get_header_line_loc(line);

    /* Depth */
    char **depth = _mats_input_get_header_line_depth(line);

    headerLine -> depth = depth[0];
    headerLine -> index = depth[1];

    /* Dimensions */
    char **dim = _mats_input_get_header_line_dim(line);

    headerLine -> dim[0] = dim[0];
    headerLine -> dim[1] = dim[1];

    /* Matrix type */
    headerLine -> type = _mats_input_get_header_line_mtype(line);

    /* Label */
    headerLine -> label = _mats_input_get_header_line_label(line);

    /* Free memory */
    free(line);

    free(type);
    free(depth);
    free(dim);

    return headerLine;
}


/*  ------------------------------------------------------------------------------------------------------  */


static char *_mats_input_get_header_line(FILE *stream)
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


static bool _mats_input_get_header_line_empty(char *line)
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


static bool _mats_input_get_header_line_end(char *line)
{
    /*

        CHeck if the line signals the end of the header (128 x '=')

    */

    size_t headerBarrierLen = 0;
    size_t lineLen = strlen(line);

    /* Line is smaller than max header line length -> can return false */
    if (lineLen < HEADER_LINE_LEN)
        return false;

    /* Check how many '=' are in line */
    for (size_t i = 0; i < lineLen; i++)
      {
        if (line[i] == '=')
            headerBarrierLen ++;
      }

    /* If entire line is made up of '=' return true */
    if (headerBarrierLen == HEADER_LINE_LEN)
        return true;

    return false;
}


/*  ------------------------------------------------------------------------------------------------------  */


static char *_mats_input_get_header_line_type(char *line)
{
    /*

        Check which type of information the line contains

        TODO: Problems arise if the labels contain those keywords...

    */

    char *type = malloc(2);

    /* Böock matrix */
    if (misc_sin_rm(line, "block matrix at "))
      {
        type[0] = 'b';

        /* Binary file */
        if (misc_sin_rm(line, "byte "))
          {
            type[1] = 'b';
          }

        /* Binary file */
        else if (misc_sin_rm(line, "line "))
          {
            type[1] = 'l';
          }

        else
            goto lineError;

        return type;
      }

    /* Ordinary matrix */
    if (misc_sin_rm(line, "matrix at "))
      {
        type[0] = 'm';

        /* Binary file */
        if (misc_sin_rm(line, "byte "))
          {
            type[1] = 'b';
          }

        /* Binary file */
        else if (misc_sin_rm(line, "line "))
          {
            type[1] = 'l';
          }

        else
            goto lineError;

        return type;
      }

    lineError:

    printf("Expected line to be of shape '%s', '%s', '%s' or '%s' but instead got '%s'!\n"
           , "block matrix at line <> : <>", "block matrix at byte <> : <>", "matrix at line <> : <>", "matrix at byte <> : <>", line);

    free(type);
    exit(1);

    return NULL;
}


static char *_mats_input_get_header_line_loc(char *line)
{
    /*

        Get the location of the matrix

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


static char **_mats_input_get_header_line_depth(char *line)
{
    /*

        Get the depth and the index of the matrix

    */

    char c;

    bool getDepth = false;

    char *depth = NULL;
    size_t depthLen = 0;

    char *indexAtDepth = NULL;
    size_t indexAtDepthLen = 0;

    for (int i = (int) strlen(line) - 1; i > -1; i--)
      {
        /* Current character */
        c = line[i];

        /* Switch the flag at the first delimiter */
        if (c == ',' && !getDepth)
            getDepth = true;

        /* Break at the second delimiter */
        else if (c == ',' && getDepth)
          {
            /* Shorten the string */
            line[i] = '\0';
            break;
          }

        /* Only keep digits */
        if (isdigit(c))
          {
            if (getDepth)
              {
                depth = realloc(depth, ++depthLen + 1);
                depth[depthLen - 1] = c;
                depth[depthLen] = '\0';
              }

            else
              {
                indexAtDepth = realloc(indexAtDepth, ++indexAtDepthLen + 1);
                indexAtDepth[indexAtDepthLen - 1] = c;
                indexAtDepth[indexAtDepthLen] = '\0';
              }
          }
      }

    /* Reverse the strings */
    for (int i = 0, j = (int) depthLen - 1; i <= j; i++, j--)
      {
        char ch = depth[i];
        depth[i] = depth[j];
        depth[j] = ch;
      }

    for (int i = 0, j = (int) indexAtDepthLen - 1; i <= j; i++, j--)
      {
        char ch = indexAtDepth[i];
        indexAtDepth[i] = indexAtDepth[j];
        indexAtDepth[j] = ch;
      }

    /* Return the strings */
    char **depths = malloc(sizeof(char*) * 2);
    depths[0] = depth;
    depths[1] = indexAtDepth;

    return depths;
}


static char **_mats_input_get_header_line_dim(char *line)
{
    /*

        Get the dimensions of the matrix

    */

    char c;

    bool getDim1 = false;

    char *dim1 = NULL;
    size_t dim1Len = 0;

    char *dim2 = NULL;
    size_t dim2Len = 0;

    for (int i = (int) strlen(line) - 1; i > -1; i--)
      {
        /* Current character */
        c = line[i];

        /* Break at the delimiter */
        if (c == ',' && !getDim1)
            getDim1 = true;

        else if (c == ',' && getDim1)
          {
            /* Shorten the string */
            line[i] = '\0';
            break;
          }

        /* Only keep digits */
        if (isdigit(c))
          {
            if (getDim1)
              {
                dim1 = realloc(dim1, ++dim1Len + 1);
                dim1[dim1Len - 1] = c;
                dim1[dim1Len] = '\0';
              }

            else
              {
                dim2 = realloc(dim2, ++dim2Len + 1);
                dim2[dim2Len - 1] = c;
                dim2[dim2Len] = '\0';
              }
          }
      }

    /* Reverse the strings */
    for (int i = 0, j = (int) dim1Len - 1; i <= j; i++, j--)
      {
        char ch = dim1[i];
        dim1[i] = dim1[j];
        dim1[j] = ch;
      }

    for (int i = 0, j = (int) dim2Len - 1; i <= j; i++, j--)
      {
        char ch = dim2[i];
        dim2[i] = dim2[j];
        dim2[j] = ch;
      }

    /* Return the strings */
    char **dim = malloc(sizeof(char*) * 2);
    dim[0] = dim1;
    dim[1] = dim2;

    return dim;
}


static char *_mats_input_get_header_line_mtype(char *line)
{
    /*

        Get the type of the matrix

    */

    char c;

    char *type = NULL;
    size_t typeLen = 0;

    for (int i = (int) strlen(line) - 1; i > -1; i--)
      {
        /* Current character */
        c = line[i];

        /* Break at the delimiter */
        if (c == ',')
          {
            /* Shorten the string */
            line[i] = '\0';
            break;
          }

        /* Skip white space */
        if (isspace(c))
            continue;

        /* Set the type */
        type = realloc(type, ++typeLen + 1);
        type[typeLen - 1] = c;
        type[typeLen] = '\0';
      }

    /* Reverse the string */
    for (int i = 0, j = (int) typeLen - 1; i <= j; i++, j--)
      {
        char ch = type[i];
        type[i] = type[j];
        type[j] = ch;
      }

    return type;
}


static char *_mats_input_get_header_line_label(char *line)
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
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     MATRICES STRUCT     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   New / Free Struct   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *mats_new(size_t size)
{
    /*

        Create a new collection of matrices

    */

    /* Matrices struct */
    mats_t *mats = malloc(sizeof(mats_t));

    /* Set number of matrices */
    mats -> size = size;

    /* Matrix */
    if (size == 0)
        mats -> mat = NULL;

    else
      {
        mats -> mat = malloc(sizeof(mat_t*) * size);

        for (size_t i = 0; i < size; i++)
            mats -> mat[i] = NULL;
      }

    return mats;
}


/*  ------------------------------------------------------------------------------------------------------  */


mats_t *mats_free(mats_t *mats)
{
    /*

        Free mats

    */

    if (mats == NULL)
        return NULL;

    free(mats -> mat);
    free(mats);

    return NULL;
}


mats_t *mats_free_full(mats_t *mats)
{
    /*

        Free mats at its matrices

    */

    if (mats == NULL)
        return NULL;

    for (size_t i = 0; i < mats -> size; i++)
        mats -> mat[i] = mat_free(mats -> mat[i]);

    free(mats -> mat);
    free(mats);

    return NULL;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Setters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mats_set_mat(mats_t *mats, size_t loc, mat_t *mat)
{
    /*

        Insert a matrix into a collection of matrices at location loc.

    */

    /* Must reallocate memory to have location within bounds */
    if (loc >= mats -> size)
      {
        mats -> size = loc + 1;
        mats -> mat = realloc(mats -> mat, sizeof(mat_t*) * mats -> size);
      }

    mats -> mat[loc] = mat;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Getters   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mat_t *mats_get_mat(mats_t *mats, size_t loc)
{
    /*

        Get the matrix of mats at location loc

    */

    if (loc >= mats -> size)
        return NULL;

    return mats -> mat[loc];
}


mat_t *mats_get_mat_by_label(mats_t *mats, const char *label)
{
    /*

        Get the matrix of mats by label

    */

    for (size_t i = 0; i < mats -> size; i++)
      {
        mat_t *mat = mats_get_mat(mats, i);

        if (mat -> label != NULL && !strcmp(mat -> label, label))
            return mat;
      }

    return NULL;
}


size_t mats_get_size(mats_t *mats)
{
    /*

        Get the size of the collection of matrices

    */

    return mats -> size;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Concatenate Matrices   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


mats_t *mats_cat_mat(size_t num, ...)
{
    /*

        Concatenate matrices to a mats_t struct

    */

    va_list ap;

    size_t size = 0;
    mats_t *matsCat = mats_new(size);

    va_start(ap, num);

    for (size_t i = 0; i < num; i++)
      {
        mat_t *mat = va_arg(ap, mat_t*);

        if (mat == NULL)
            continue;

        mats_set_mat(matsCat, size++, mat);
      }

    va_end(ap);

    return matsCat;
}


mats_t *mats_cat_mats(size_t num, ...)
{
    /*

        Concatenate matrix collections (mats_t structs) to a single matrix collection

    */

    va_list ap;

    size_t size = 0;
    mats_t *matsCat = mats_new(size);

    va_start(ap, num);

    for (size_t i = 0; i < num; i++)
      {
        mats_t *mats = va_arg(ap, mats_t*);

        if (mats == NULL)
            continue;

        for (size_t j = 0; j < mats -> size; j++)
          {
            mats_set_mat(matsCat, size++, mats -> mat[j]);
          }
      }

    va_end(ap);

    return matsCat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ----------------------------------------   Append a Matrix   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int mats_append_empty(mats_t *mats, char *type, size_t *dim, bool block, char *label)
{
    /*

        Append an empty matrix to a collection of matrices

    */

    /* New matrix with label */
    mat_t *mat = mat_new(type, dim, block);
    mat_set_label(mat, label);

    /* Insert the new matrix at the end of mats */
    mats_set_mat(mats, mats -> size, mat);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

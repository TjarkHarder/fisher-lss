#ifndef MAT_H_INCLUDED
#define MAT_H_INCLUDED

#include "common.h"

#include "misc.h"



/*  ----------------------------------------------------  */
/*  ----------------  Macro Definitions  ---------------  */
/*  ----------------------------------------------------  */


#ifndef HEADER_LEN
#define HEADER_LEN 16384 // 2^14
#endif


#ifndef HEADER_LINE_LEN
#define HEADER_LINE_LEN 128 // 2^7
#endif



/*  ----------------------------------------------------  */
/*  -----------------  LAPACK Routines  ----------------  */
/*  ----------------------------------------------------  */


void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);

void dsytrf_(char *UPLO, int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
void dsytri_(char *UPLO, int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);

void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);




/*  ----------------------------------------------------  */
/*  --------------------  Structures  ------------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Type struct of a matrix

    */

    /* Identifier of matrix (e.g. full / f; diagonal / d, symmetric / s, upper-triangular / ut, lower-triangular / lt */
    const char *id;

    /* Category of the type (0: Removed entries are zeros; 1: Removed entries are repeated entries) */
    int category;


    // Must check return values!
      // 0: element exists
      // -1: element does not exist
      // 1: error

    /* Function transforming (i,j) position in matrix to index in matrix array */
    int (*index)(size_t*, size_t*, size_t*);

    /* Function transforming index in matrix array to (i,j) position in matrix */
    int (*loc)(size_t*, size_t, size_t*);

    /* Check type by comparing elements (i.e. if certain number of elements are zero, return false: e.g. off-diagonal is zero -> diagonal matrix) */
    bool (*comp)(void*, size_t*, double);

    /* Row and Column bounds, given a Column or Row Index, respectively */
    int (*rbounds)(size_t*, size_t, size_t*);
    int (*cbounds)(size_t*, size_t, size_t*);

} mat_type_t;


typedef struct
{
    /*

        Further matrix parameters in case the matrix is a block matrix

    */

    /* Function that gives the correct representation (as it would appear in a full matrix) of a block at some location (i,j) */
    void *(*rep)(void*, size_t*);

    /* Dimensions of blocks */
    size_t **bdim; // Shape (2, dim[i]) : dimensions of blocks in column / row
    size_t **bfdim; // Shape (2, dim[i]) : total dimensions of blocks in column / row

    size_t *fdim; // Shape (2) : total dimensions

} mat_block_t;


typedef struct
{
    /*

        Single location trafo function

    */

    const char *id;

    int (*loc)(void*, size_t*);
    int (*dim)(void*, size_t*);

} mat_ltrf_t;


typedef struct
{
    /*

        Keep track of location transformations (e.g. transpose)

    */

    size_t size;
    mat_ltrf_t **mltrf;

} mat_ltrfs_t;


typedef struct
{
    /*

        Matrix from file

    */

    char *fileName;

    char *type;
    size_t loc;
    bool binary;

} mat_file_t;


typedef struct
{
    /*

        Matrix struct for either one matrix or several matrices

        For several matrices matrix is pointer to **mat_t struct with *dim = number of matrices and type = collection / c

    */

    /* Label of the matrix */
    char *label;


    /* Location trasnformations */
    mat_ltrfs_t *mltrfs;

    /* Type of the matrix */
    mat_type_t *mtype;

    /* Block parameters of the matrix (if matrix is in block form) */
    mat_block_t *mblock;

    /* Matrix from file */
    mat_file_t *mfile;


    /* Pointer to either **mat_t struct or *double array */
    void *matrix;

    /* Dimensions of matrix */
    size_t *dim;

} mat_t;


typedef struct
{
    /*

        Reduce a matrix struct

    */

    /* Allow a small error when comparing an element with another or zero */
    double err;

} mat_reduce_t;



typedef struct
{
    /*

        Collection of matrices

    */

    size_t size;
    mat_t **mat;

} mats_t;




/*  ----------------------------------------------------  */
/*  ---------    Get the Matrix Type Struct   ----------  */
/*  ----------------------------------------------------  */


mat_type_t *mat_type_get_struct(
    const char *id);

mat_type_t *mat_type_get_struct_tp(
    const char *id);



/*  ----------------------------------------------------  */
/*  --------    Get the Location Trafo Struct   --------  */
/*  ----------------------------------------------------  */


mat_ltrf_t *mat_ltrf_get_struct(
    const char *id);





/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%    MATRIX STRUCT    %%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  -----------------    New / Free   ------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_new(
    const char *type,
    size_t *dim,
    bool block);

mat_t *mat_fnew(
    const char *type,
    size_t *dim,
    bool block);

mat_t *mat_new_id(
    mat_t *mat);

/*  ----------------------------------------------------  */

mat_t *mat_free(
    mat_t *mat);

mat_t *mat_ffree(
    mat_t *mat);



/*  ----------------------------------------------------  */
/*  ------------------    Compare    -------------------  */
/*  ----------------------------------------------------  */


bool mat_comp_matrix(
    mat_t *mat1,
    mat_t *mat2,
    double tol);

bool mat_comp_bstruct(
    mat_t *mat1,
    mat_t *mat2);

bool mat_comp_bstruct_cr(
    mat_t *mat1,
    mat_t *mat2);



/*  ----------------------------------------------------  */
/*  ------------------    Setters    -------------------  */
/*  ----------------------------------------------------  */


int mat_set_label(
    mat_t *mat,
    char *label);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int mat_set_value(
    mat_t *mat,
    size_t *loc,
    double value);

int mat_set_matrix(
    mat_t *mat1,
    mat_t *mat2);

int mat_set_array(
    mat_t *mat,
    double *array);


int mat_set_block(
    mat_t *mat,
    size_t *loc,
    mat_t *matBlock);

int mat_fset_block(
    mat_t *mat,
    size_t *loc,
    mat_t *matBlock);



/*  ----------------------------------------------------  */
/*  ------------------    Getters    -------------------  */
/*  ----------------------------------------------------  */


char *mat_get_label(
    mat_t *mat);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int mat_get_dim(
    mat_t *mat,
    size_t *dimMat);

int mat_mget_dim(
    mat_t *mat,
    size_t *dimMem);


int mat_get_bdim(
    mat_t *mat,
    size_t depth,
    size_t **loc,
    size_t *bdimMat);


int mat_get_bfdim(
    mat_t *mat,
    size_t *loc,
    size_t *bfdimMat);

int mat_mget_bfdim(
    mat_t *mat,
    size_t *loc,
    size_t *bfdimMem);


int mat_get_fdim(
    mat_t *mat,
    size_t *fdimMat);

int mat_mget_fdim(
    mat_t *mat,
    size_t *fdimMem);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int mat_get_rbounds(
    mat_t *mat,
    size_t col,
    size_t *rbounds);

int mat_get_cbounds(
    mat_t *mat,
    size_t row,
    size_t *cbounds);

/*  ----------------------------------------------------  */

int mat_get_loc(
    mat_t *mat,
    size_t indexMat,
    size_t *locMat);

int mat_mget_loc(
    mat_t *mat,
    size_t *locMat,
    size_t *locMem);


int mat_get_index(
    mat_t *mat,
    size_t *locMat,
    size_t *indexMat);

int mat_mget_index(
    mat_t *mat,
    size_t *locMat,
    size_t *indexMem);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double mat_get_value(
    mat_t *mat,
    size_t *loc);

double mat_get_ivalue(
    mat_t *mat,
    size_t *loc,
    size_t *index);

/*  ----------------------------------------------------  */

mat_t *mat_get_block(
    mat_t *mat,
    size_t *loc);

mat_t *mat_get_iblock(
    mat_t *mat,
    size_t *loc,
    size_t *index);


mat_t *mat_fget_block(
    mat_t *mat,
    size_t *loc);

mat_t *mat_fget_iblock(
    mat_t *mat,
    size_t *loc,
    size_t *index);


mat_t *mat_mget_block(
    mat_t *mat,
    size_t *loc);

mat_t *mat_mget_iblock(
    mat_t *mat,
    size_t *loc,
    size_t *index);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


size_t mat_get_depth(
    mat_t *mat);



/*  ----------------------------------------------------  */
/*  -----------------    Submatrix    ------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_sub(
    mat_t *mat,
    size_t *indices1,
    size_t size1,
    size_t *indices2,
    size_t size2,
    mat_reduce_t *reduce);



/*  ----------------------------------------------------  */
/*  ----------------    Copy Matrix    -----------------  */
/*  ----------------------------------------------------  */


mat_t *mat_cp(
    mat_t *mat);

mat_t *mat_fcp(
    mat_t *mat);


mat_t *mat_cp_struct(
    mat_t *mat,
    const char *id,
    int idDepth);

mat_t *mat_cp_structs(
    mats_t *mats,
    const char *id,
    int idDepth);

mat_t *mat_cp_struct_rc(
    mat_t *mat1,
    mat_t *mat2);


mat_t *mat_cp_struct_block(
    mat_t *mat,
    size_t depth,
    size_t **loc,
    const char *id,
    int idDepth);

mat_t *mat_cp_structs_block(
    mats_t *mats,
    size_t depth,
    size_t **loc,
    const char *id,
    int idDepth);

mat_t *mat_cp_struct_block_rc(
    mat_t *mat1,
    mat_t *mat2,
    size_t depth,
    size_t **loc);


mat_t *mat_cp_struct_sim(
    mat_t *mat1,
    mat_t *mat2,
    const char *id,
    int idDepth);



/*  ----------------------------------------------------  */
/*  ---------------    Swap Matrices    ----------------  */
/*  ----------------------------------------------------  */


int mat_swap(
    mat_t *mat1,
    mat_t *mat2);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


double mat_pnorm(
    mat_t *mat,
    size_t p);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_reduce(
    mat_t *mat,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


size_t mat_count(
    mat_t *mat);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_diag(
    mat_t *mat);

mat_t *mat_diag_ip(
    mat_t *mat);

/*  ----------------------------------------------------  */

mat_t *mat_diag_func(
    mat_t *mat,
    double (*func)(double));

double mat_diag_func_sqrt(
    double val);

double mat_diag_func_inv(
    double val);

double mat_diag_func_invsqrt(
    double val);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_trafo(
    mat_t *mat,
    const char *type,
    size_t *dim,
    bool block,
    mat_reduce_t *reduce);

mat_t *mat_trafo_ip(
    mat_t *mat,
    const char *type,
    size_t *dim,
    bool block,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int mat_print(
    mat_t *mat,
    int precision);



/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%    LINEAR ALGEBRA    %%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  -------    Basic Linear Algebra Functions    -------  */
/*  ----------------------------------------------------  */


mat_t *mat_tp(
    mat_t *mat,
    mat_reduce_t *reduce);

mat_t *mat_ftp(
    mat_t *mat,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_add_ip(
    char job,
    mat_t *mat1,
    mat_t *mat2,
    mat_t *matAdd,
    mat_reduce_t *reduce);

mat_t *mat_add(
    char job,
    mat_t *mat1,
    mat_t *mat2,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_mult(
    mat_t *mat1,
    mat_t *mat2,
    double value,
    mat_reduce_t *reduce);

mat_t *mat_mult_ip(
    char job,
    mat_t *mat1,
    mat_t *mat2,
    mat_t *matMult,
    double value,
    mat_reduce_t *reduce);

/*  ----------------------------------------------------  */

mat_t *mat_mult_mat(
    mat_t *mat1,
    mat_t *mat2,
    mat_reduce_t *reduce);

mat_t *mat_mult_mat_ip(
    char job,
    mat_t *mat1,
    mat_t *mat2,
    mat_t *matMult,
    mat_reduce_t *reduce);

/*  ----------------------------------------------------  */

mat_t *mat_mult_val(
    mat_t *mat,
    double value,
    mat_reduce_t *reduce);

mat_t *mat_mult_val_ip(
    mat_t *mat,
    double value,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_perm(
    char job,
    mat_t *mat,
    size_t *perm,
    mat_reduce_t *reduce);



/*  ----------------------------------------------------  */
/*  -------------    Correlation Matrix    -------------  */
/*  ----------------------------------------------------  */


mat_t *mat_corr(
    mat_t *mat,
    mat_reduce_t *reduce);

mat_t *mat_corr_ip(
    mat_t *mat,
    mat_reduce_t *reduce);



/*  ----------------------------------------------------  */
/*  -----------    Matrix Decompositions    ------------  */
/*  ----------------------------------------------------  */


mats_t *mat_qr(
    mat_t *mat,
    mat_reduce_t *reduce);

/*  ----------------------------------------------------  */

mats_t *mat_eig(
    mat_t *mat,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  --   Linear Equations Solver + Matrix Inversion   --  */
/*  ----------------------------------------------------  */


mat_t *mat_backsub(
    mat_t *matA,
    mat_t *matB,
    mat_reduce_t *reduce);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mat_t *mat_inv(
    mat_t *mat,
    mat_reduce_t *reduce);

mat_t *mat_inv_lapack(
    mat_t *mat,
    mat_reduce_t *reduce);

/*  ----------------------------------------------------  */

mat_t *mat_inv_cov(
    mat_t *inv(mat_t*, mat_reduce_t*),
    mat_t *mat,
    mat_reduce_t *reduce);



/*  ----------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%    MATRICES STRUCT    %%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ----------------------------------------------------  */


/*  ----------------------------------------------------  */
/*  -----------------    New / Free   ------------------  */
/*  ----------------------------------------------------  */


mats_t *mats_new(
    size_t size);

/*  ----------------------------------------------------  */

mats_t *mats_free(
    mats_t *mats);

mats_t *mats_free_full(
    mats_t *mats);



/*  ----------------------------------------------------  */
/*  ------------------    Setters    -------------------  */
/*  ----------------------------------------------------  */


int mats_set_mat(
    mats_t *mats,
    size_t loc,
    mat_t *mat);



/*  ----------------------------------------------------  */
/*  ------------------    Getters    -------------------  */
/*  ----------------------------------------------------  */


mat_t *mats_get_mat(
    mats_t *mats,
    size_t loc);

mat_t *mats_get_mat_by_label(
    mats_t *mats,
    const char *label);

/*  ----------------------------------------------------  */

size_t mats_get_size(
    mats_t *mats);



/*  ----------------------------------------------------  */
/*  ------------    Concatenate Matrices    ------------  */
/*  ----------------------------------------------------  */


mats_t *mats_cat_mat(
    size_t num,
    ...);

mats_t *mats_cat_mats(
    size_t num,
    ...);



/*  ----------------------------------------------------  */
/*  --------------    Append a Matrix    ---------------  */
/*  ----------------------------------------------------  */


int mats_append_empty(
    mats_t *mats,
    char *type,
    size_t *dim,
    bool block,
    char *label);



/*  ----------------------------------------------------  */
/*  ----------    Input / Output Matrices    -----------  */
/*  ----------------------------------------------------  */


int mats_output(
    char *fileName,
    mats_t *mats,
    int *precision);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


mats_t *mats_input(
    char *fileName,
    char *label,
    bool fill);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // MAT_H_INCLUDED

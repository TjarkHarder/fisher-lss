#ifndef SHAPE_H_INCLUDED
#define SHAPE_H_INCLUDED


#include "common.h"

#include "misc.h"



/*  ----------------------------------------------------  */
/*  -------------------   Structures   -----------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Shape struct (eg line, triangle, ...)

    */

    size_t dim;

    double *length;
    double *orientation;
    double *angle;

    double volume;
    double volumeErr;

    bool parity; // Parity invariance

} shape_t;



/*  ----------------------------------------------------  */
/*  -------------------   New / Free   -----------------  */
/*  ----------------------------------------------------  */


shape_t *shape_new(
    size_t dim);

shape_t *shape_free(
    shape_t *shape);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


shape_t *shape_cp(
    shape_t *shape);

shape_t *shape_cp_parity(
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  ----------------------   Swap   --------------------  */
/*  ----------------------------------------------------  */


int shape_swap(
    shape_t **shape1,
    shape_t **shape2);

int shape_swap_vertex(
    shape_t *shape,
    size_t index1,
    size_t index2);



/*  ----------------------------------------------------  */
/*  --------------------   Setters   -------------------  */
/*  ----------------------------------------------------  */


bool shape_set_vertex_length(
    shape_t *shape,
    size_t index,
    double length);

bool shape_set_vertex_orientation(
    shape_t *shape,
    size_t index,
    double orientation);

bool shape_set_vertex_angle(
    shape_t *shape,
    size_t index1,
    size_t index2,
    double angle);


bool shape_set_volume(
    shape_t *shape,
    double volume);

bool shape_set_volume_err(
    shape_t *shape,
    double volumeErr);


bool shape_set_parity(
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  --------------------   Getters   -------------------  */
/*  ----------------------------------------------------  */


double shape_get_vertex_length(
    shape_t *shape,
    size_t index);

double shape_qget_vertex_length(
    shape_t *shape,
    size_t index);

/*  ----------------------------------------------------  */

double shape_get_vertex_orientation(
    shape_t *shape,
    size_t index);

double shape_qget_vertex_orientation(
    shape_t *shape,
    size_t index);

/*  ----------------------------------------------------  */

double shape_get_vertex_angle(
    shape_t *shape,
    size_t index1,
    size_t index2);

double shape_qget_vertex_angle(
    shape_t *shape,
    size_t index1,
    size_t index2);

/*  ----------------------------------------------------  */

double shape_get_volume(
    shape_t *shape);

double shape_get_volume_err(
    shape_t *shape);


bool shape_get_parity(
    shape_t *shape);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


size_t shape_get_vertex_angle_index(
    size_t dim,
    size_t index1,
    size_t index2);



/*  ----------------------------------------------------  */
/*  ------------------   Check Shape   -----------------  */
/*  ----------------------------------------------------  */


bool shape_check(
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  ------------------   Check Shape   -----------------  */
/*  ----------------------------------------------------  */


int shape_parity(
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  ------------------   Sort Shape   ------------------  */
/*  ----------------------------------------------------  */


int shape_sort(
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  ----------------   Compare Shapes   ----------------  */
/*  ----------------------------------------------------  */


bool shape_ecomp(
    shape_t *shape1,
    shape_t *shape2);

int shape_comp(
    shape_t *shape1,
    shape_t *shape2);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


bool shape_comp_vertex_bin(
    shape_t *shape1,
    size_t index1,
    shape_t *shape2,
    size_t index2,
    double lengthTol,
    double orientationTol);

bool **shape_comp_vertexes_bin(
    shape_t *shape1,
    shape_t *shape2,
    double lengthTol,
    double orientationTol);

/*  ----------------------------------------------------  */

bool shape_comp_vertex(
    shape_t *shape1,
    size_t index1,
    shape_t *shape2,
    size_t index2);

bool **shape_comp_vertexes(
    shape_t *shape1,
    shape_t *shape2);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


bool shape_comp_vertex_parity_bin(
    shape_t *shape1,
    size_t index1,
    shape_t *shape2,
    size_t index2,
    double lengthTol,
    double orientationTol);

bool **shape_comp_vertexes_parity_bin(
    shape_t *shape1,
    shape_t *shape2,
    double lengthTol,
    double orientationTol);

/*  ----------------------------------------------------  */

bool shape_comp_vertex_parity(
    shape_t *shape1,
    size_t index1,
    shape_t *shape2,
    size_t index2);

bool **shape_comp_vertexes_parity(
    shape_t *shape1,
    shape_t *shape2);



/*  ----------------------------------------------------  */
/*  ---------------   Degenerate Edges   ---------------  */
/*  ----------------------------------------------------  */


size_t *shape_degen_vertexes(
    shape_t *shape);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif  // SHAPE_H_INCLUDED

#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include "common.h"

#include "integrate.h"



/*  ----------------------------------------------------  */
/*  ------------------  Time Structure  ----------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Time in __:__:__:___ format

    */

    unsigned int hr;
    unsigned int min;
    unsigned int sec;
    unsigned int msec;

} misc_time_t;



/*  ----------------------------------------------------  */
/*  ---------------   Extern Functions   ---------------  */
/*  ----------------------------------------------------  */


extern double integrand_growth(
    double *var,
    size_t dim,
    void *params);



/*  ----------------------------------------------------  */
/*  -------------   Skip Lines in Stream   -------------  */
/*  ----------------------------------------------------  */


int misc_stream_skip(
    FILE *stream,
    size_t lines);



/*  ----------------------------------------------------  */
/*  --------------   Concatenate Strings   -------------  */
/*  ----------------------------------------------------  */


char *misc_scat(
    size_t size, ...);



/*  ----------------------------------------------------  */
/*  -----------------   Copy String   ------------------  */
/*  ----------------------------------------------------  */

int misc_scp(
    char **str1,
    const char *str2);



/*  ----------------------------------------------------  */
/*  ----------   Count Characters in String   ----------  */
/*  ----------------------------------------------------  */


size_t misc_scnch(
    const char *str,
    char ch);



/*  ----------------------------------------------------  */
/*  --------   Search for String within String   -------  */
/*  ----------------------------------------------------  */


bool misc_sin_rm(
    char *str1,
    const char *str2);

/*  ----------------------------------------------------  */

bool misc_sin(
    const char *str1,
    const char *str2);

bool misc_sins(
    char **array,
    size_t size,
    const char *str,
    size_t *index);

/*  ----------------------------------------------------  */

bool misc_sinci(
    const char *str1,
    const char *str2);

/*  ----------------------------------------------------  */

char *misc_ssrm(
    char *str,
    const char *sub);

bool misc_srm(
    char **array,
    size_t *size,
    const char *str);



/*  ----------------------------------------------------  */
/*  ----------   Search for String in Array   ----------  */
/*  ----------------------------------------------------  */


bool misc_ssearch(
    char **array,
    size_t size,
    const char *str,
    size_t *index);

size_t misc_ssearch_index(
    char **array,
    size_t size,
    const char *str,
    bool *success);



/*  ----------------------------------------------------  */
/*  --------   Turn Double to String and Back   --------  */
/*  ----------------------------------------------------  */


char *misc_dtos(
    double num,
    unsigned int tol);

double misc_stod(
    const char *str,
    char **remstr,
    bool *success);



/*  ----------------------------------------------------  */
/*  ----------------   Binary Search   -----------------  */
/*  ----------------------------------------------------  */


bool misc_bsearch(
    double *array,
    size_t size,
    double value,
    double tol,
    size_t *index);

size_t misc_bsearch_index(
    double *array,
    size_t size,
    double value,
    double tol,
    bool *success);

/*  ----------------------------------------------------  */

int misc_insort(
    double **array,
    size_t *size,
    double value,
    double tol,
    bool *success);



/*  ----------------------------------------------------  */
/*  -----------------   Swap Pointers   ----------------  */
/*  ----------------------------------------------------  */


int misc_swap(
    void *a,
    void *b,
    const char *type);



/*  ----------------------------------------------------  */
/*  ------------------   Sort Array   ------------------  */
/*  ----------------------------------------------------  */


size_t *misc_sort(
    double *array,
    size_t size);



/*  ----------------------------------------------------  */
/*  ----------------   Random Numbers   ----------------  */
/*  ----------------------------------------------------  */


size_t *misc_random_sample(
    size_t start,
    size_t stop,
    size_t *size,
    bool repetition);



/*  ----------------------------------------------------  */
/*  ----------------   Combinatorics   -----------------  */
/*  ----------------------------------------------------  */

size_t misc_factorial(
    size_t n);

size_t misc_binom(
    size_t n,
    size_t m);

size_t misc_multinom(
    size_t n,
    size_t *m,
    size_t size);



/*  ----------------------------------------------------  */
/*  ---------------------   Time   ---------------------  */
/*  ----------------------------------------------------  */


misc_time_t *misc_time_new();

misc_time_t *misc_time_free(
    misc_time_t *time);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int misc_time_format(
    misc_time_t *time,
    double timeSec);

misc_time_t *misc_time_format_out(
    double timeSec);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // MISC_H_INCLUDED

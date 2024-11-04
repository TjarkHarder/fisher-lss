#ifndef PRINT_H_INCLUDED
#define PRINT_H_INCLUDED

#include "common.h"

#include "misc.h"
#include "kernels.h"


/*  ----------------------------------------------------  */
/*  -------------------  Structures   ------------------  */
/*  ----------------------------------------------------  */


typedef struct
{
    /*

        Print flags

    */

    bool main;
    bool sub;
    bool fid;

} print_flags_t;


typedef struct
{
    /*

        Parameters for the print functions

    */

    /* Id of the process */
    const char *id;


    /* Print flags */
    print_flags_t *flags;


    /* Start and end of the timer */
    double startTimer;
    double endTimer;

    /* Execution time and remaining time */
    double execTimeStep;
    double execTime;

    /* Formatted execution time and remaining time */
    misc_time_t *execTimeFormat;
    misc_time_t *remTimeFormat;


    /* Thread */
    int thread;

    /* Progress */
    double progressStep;
    double progress;

    /* State of the calculation */
    int subFinished;
    int stageFinished;

    /* Kern */
    kern_t *kern;


    /* Pointer to additional variables */
    void *var;

    /* Pointer to data sizes */
    size_t *sizes;

} print_t;



/*  ----------------------------------------------------  */
/*  ---------------  Print Flags Struct   --------------  */
/*  ----------------------------------------------------  */


print_flags_t *print_flags_new(void);

print_flags_t *print_flags_free(
    print_flags_t *flags);

/*  ----------------------------------------------------  */

print_flags_t *print_flags_cp(
    print_flags_t *flags);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int print_flags_set_main(
    print_flags_t *flags,
    bool main);

int print_flags_set_sub(
    print_flags_t *flags,
    bool sub);

int print_flags_set_fid(
    print_flags_t *flags,
    bool fid);



/*  ----------------------------------------------------  */
/*  ------------------  Print Struct   -----------------  */
/*  ----------------------------------------------------  */


print_t *print_new(void);

print_t *print_free(
    print_t *print);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int print_set_id(
    print_t *print,
    const char *id);

int print_set_flags(
    print_t *print,
    print_flags_t *flags);

/*  ----------------------------------------------------  */

int print_set_thread(
    print_t *print,
    int thread);

int print_set_kern(
    print_t *print,
    kern_t *kern);

/*  ----------------------------------------------------  */

int print_set_progress(
    print_t *print,
    double progress);

int print_set_subfinished(
    print_t *print,
    int subFinished);

int print_set_stagefinished(
    print_t *print,
    int stageFinished);

/*  ----------------------------------------------------  */

int print_set_var(
    print_t *print,
    void *var);

int print_set_sizes(
    print_t *print,
    size_t *sizes);


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


int print_update_progress(
    print_t *print,
    double progress);




/*  ----------------------------------------------------  */
/*  ----------   Extern Variables/Functions   ----------  */
/*  ----------------------------------------------------  */

/* -- */


/*  ----------------------------------------------------  */
/*  ------------------   Initialise   ------------------  */
/*  ----------------------------------------------------  */


int prnt_ini(void);
int prnt_free(void);



/*  ----------------------------------------------------  */
/*  ----------------   Print Functions   ---------------  */
/*  ----------------------------------------------------  */


int print_spec(
    print_t *print);

int print_cov(
    print_t *print);

int print_fish(
    print_t *print);



/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


#endif // PRINT_H_INCLUDED

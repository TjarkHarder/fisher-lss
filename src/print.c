/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PRINT.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "print.h"


// TODO: Improve this...

/* Size of print box */
const unsigned int boxSize = 150;





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     STRUCT FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Print Flags Struct   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


print_flags_t *print_flags_new(void)
{
    /*

        Create a new print_flags_t struct

    */

    print_flags_t *flags = malloc(sizeof(print_flags_t));

    flags -> main = false;
    flags -> sub = false;
    flags -> fid = false;

    return flags;
}


print_flags_t *print_flags_free(print_flags_t *flags)
{
    /*

        Free flags

    */

    if (flags == NULL)
        return NULL;

    free(flags);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */


print_flags_t *print_flags_cp(print_flags_t *flags)
{
    /*

        Copy flags

    */

    print_flags_t *flagsCp = print_flags_new();

    flagsCp -> main = flags -> main;
    flagsCp -> sub = flags -> sub;
    flagsCp -> fid = flags -> fid;

    return flagsCp;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int print_flags_set_main(print_flags_t *flags, bool main)
{
    /*

        Set the main flag for flags

    */

    flags -> main = main;

    return 0;
}


int print_flags_set_sub(print_flags_t *flags, bool sub)
{
    /*

        Set the sub flag for flags

    */

    flags -> sub = sub;

    return 0;
}


int print_flags_set_fid(print_flags_t *flags, bool fid)
{
    /*

        Set the fid flag for flags

    */

    flags -> fid = fid;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Print Struct   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


print_t *print_new(void)
{
    /*

        Create a new print_t struct

    */

    print_t *print = malloc(sizeof(print_t));

    print -> id = NULL;


    print -> flags = NULL;


    print -> startTimer = 0.;
    print -> endTimer = 0.;

    print -> execTimeStep = 0.;
    print -> execTime = 0.;

    print -> execTimeFormat = misc_time_new();
    print -> remTimeFormat = misc_time_new();


    print -> thread = 0;

    print -> progressStep = 0.;
    print -> progress = 0.;

    print -> subFinished = 0;
    print -> stageFinished = 0;

    print -> kern = NULL;

    print -> var = NULL;
    print -> sizes = NULL;

    return print;
}


print_t *print_free(print_t *print)
{
    /*

        Free an print_t struct

    */

    if (print == NULL)
        return print;

    print -> flags = print_flags_free(print -> flags);

    print -> execTimeFormat = misc_time_free(print -> execTimeFormat);
    print -> remTimeFormat = misc_time_free(print -> remTimeFormat);

    free(print);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int print_set_id(print_t *print, const char *id)
{
    /*

        Set print -> id

    */

    print -> id = id;

    return 0;
}


int print_set_flags(print_t *print, print_flags_t *flags)
{
    /*

        Set print -> flags

    */

    print -> flags = print_flags_free(print -> flags);
    print -> flags = print_flags_cp(flags);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int print_set_thread(print_t *print, int thread)
{
    /*

        Set the current thread's id

    */

    print -> thread = thread;

    return 0;
}


int print_set_kern(print_t *print, kern_t *kern)
{
    /*

        Set the pointer to the kern_t struct

    */

    print -> kern = kern;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int print_set_progress(print_t *print, double progress)
{
    /*

        Set the progress in print

    */

    print -> progressStep = print -> progress - progress;
    print -> progress = progress;

    return 0;
}


int print_set_subfinished(print_t *print, int subFinished)
{
    /*

        Set print -> subFinished

    */

    print -> subFinished = subFinished;

    return 0;
}


int print_set_stagefinished(print_t *print, int stageFinished)
{
    /*

        Set print -> stageFinished

    */

    print -> stageFinished = stageFinished;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int print_set_var(print_t *print, void *var)
{
    /*

        Set the var pointer in print

    */

    print -> var = var;

    return 0;
}


int print_set_sizes(print_t *print, size_t *sizes)
{
    /*

        Set the sizes pointer in print

    */

    print -> sizes = sizes;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int print_update_progress(print_t *print, double progress)
{
    /*

        Update the progress

    */

    print -> progressStep = progress;
    print -> progress += progress;

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LOCAL PARAMETERS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



static double _execTime;





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INITIALISE / FREE MODULE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


int prnt_ini(void)
{
    /*

        Initialise the print module.

    */

    /* Global execution time */
    _execTime = 0.;

    return 0;
}


int prnt_free(void)
{
    /*

        Free the local variables.

    */

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PRINT MODULES     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



static int _print_str(const char *str)
{
    /*

        Print a string of length "boxSize" by adding white space to the end

    */

    if (strlen(str) >= boxSize)
      {
        printf("%s\n", str);
      }

    else
      {
        printf("%s%*c\n", str, (int) (boxSize - strlen(str)), ' ');
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _print_progress(print_t *print);
static int _print_fiducials(print_t *print);

static int _print_spec_sub(print_t *print);
static int _print_cov_sub(print_t *print);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Main Print Functions   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int print_spec(print_t *print)
{
    /*

        Print information to screen for the spectra

    */

    /* Print nothing */
    if (!(print -> flags -> main))
        return 0;

    /* Allocate memory */
    char *str = malloc(boxSize);

    /* Started calculating the stage */
    if (print -> stageFinished == 0)
      {
        /* Start the timer */
        print -> startTimer = omp_get_wtime();

        /* Reset the execution time */
        print -> execTime = 0.;

        /* Print the fiducials */
        _print_fiducials(print);

        /* Print a header */
        _print_str("");
        _print_str("* * *");
        _print_str("*");

        if (print -> sizes != NULL)
          {
            snprintf(str, boxSize, "*  Calculating %s for (%ld, %ld, %ld) values of (z, k, mu) for a total of %ld values.",
                                    print -> id, print -> sizes[0], print -> sizes[1], print -> sizes[2], print -> sizes[3]);
            _print_str((const char*) str);

            _print_str("*");
            _print_str("*");
          }

        /* Stage is currently being calculated (no progress here) */
        print -> stageFinished = -1;
      }

    /* Currently calculating the stage */
    else if (print -> stageFinished == -1)
      {
        /* End the timer */
        print -> endTimer = omp_get_wtime();

        /* Execution time */
        print -> execTimeStep = print -> endTimer - print -> startTimer;
        print -> execTime += print -> execTimeStep;

        /* Reset the timer */
        print -> startTimer = print -> endTimer;

        /* Print what each thread does */
        _print_spec_sub(print);

        /* Print the current progress */
        _print_progress(print);
      }

    /* Finished calculating the stage */
    else if (print -> stageFinished == 1)
      {
        /* End the timer */
        print -> endTimer = omp_get_wtime();

        /* Execution time */
        print -> execTimeStep = print -> endTimer - print -> startTimer;
        print -> execTime += print -> execTimeStep;

        /* Global execution time */
        _execTime += print -> execTime;

        /* Formatted global execution time */
        misc_time_t *timeGlobal = misc_time_format_out(_execTime);

        /* Print the current progress */
        _print_progress(print);

        /* Print finish message */
        _print_str("*");
        _print_str("*");

        snprintf(str, boxSize, "*  Finished calculating %s.",
                                print -> id);
        _print_str((const char*) str);

        _print_str("*");

        /* Total execution time */
        snprintf(str, boxSize, "*  Total Execution Time: %02uhr %02umin %02us %03ums",
                                timeGlobal -> hr, timeGlobal -> min, timeGlobal -> sec, timeGlobal -> msec);
        _print_str((const char*) str);

        /* Print a footer */
        _print_str("*");

        _print_str("* * *");

        _print_str("");
        _print_str("");

        /* Free memory */
        free(timeGlobal);
      }

    /* Free memory */
    free(str);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int print_cov(print_t *print)
{
    /*

        Print information to screen for the covariance matrices

    */

    /* Print nothing */
    if (!(print -> flags -> main))
        return 0;

    /* Allocate memory */
    char *str = malloc(boxSize);

    /* Started calculating the stage */
    if (print -> stageFinished == 0)
      {
        /* Start the timer */
        print -> startTimer = omp_get_wtime();

        /* Reset the execution time */
        print -> execTime = 0.;

        /* Print the fiducials */
        _print_fiducials(print);

        /* Print a header */
        _print_str("");
        _print_str("* * *");
        _print_str("*");

        if (print -> sizes != NULL)
          {
            snprintf(str, boxSize, "*  Calculating %s for (%ld, %ld, %ld, %ld, %ld) values of (z, k1, mu1, k2, mu2) for a total of %ld values.",
                                    print -> id, print -> sizes[0], print -> sizes[1], print -> sizes[2], print -> sizes[3], print -> sizes[4], print -> sizes[5]);
            _print_str((const char*) str);

            _print_str("*");
            _print_str("*");
          }

        /* Stage is currently being calculated (no progress here) */
        print -> stageFinished = -1;
      }

    /* Currently calculating the stage */
    else if (print -> stageFinished == -1)
      {
        /* End the timer */
        print -> endTimer = omp_get_wtime();

        /* Execution time */
        print -> execTimeStep = print -> endTimer - print -> startTimer;
        print -> execTime += print -> execTimeStep;

        /* Reset the timer */
        print -> startTimer = print -> endTimer;

        /* Print what each thread does */
        _print_cov_sub(print);

        /* Print the current progress */
        _print_progress(print);
      }

    /* Finished calculating the stage */
    else if (print -> stageFinished == 1)
      {
        /* End the timer */
        print -> endTimer = omp_get_wtime();

        /* Execution time */
        print -> execTimeStep = print -> endTimer - print -> startTimer;
        print -> execTime += print -> execTimeStep;

        /* Global execution time */
        _execTime += print -> execTime;

        /* Formatted global execution time */
        misc_time_t *timeGlobal = misc_time_format_out(_execTime);

        /* Print the current progress */
        _print_progress(print);

        /* Print finish message */
        _print_str("*");
        _print_str("*");

        snprintf(str, boxSize, "*  Finished calculating %s.",
                                print -> id);
        _print_str((const char*) str);

        _print_str("*");

        /* Total execution time */
        snprintf(str, boxSize, "*  Total Execution Time: %02uhr %02umin %02us %03ums",
                                timeGlobal -> hr, timeGlobal -> min, timeGlobal -> sec, timeGlobal -> msec);
        _print_str((const char*) str);

        /* Print a footer */
        _print_str("*");

        _print_str("* * *");

        _print_str("");
        _print_str("");

        /* Free memory */
        free(timeGlobal);
      }

    /* Free memory */
    free(str);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int print_fish(print_t *print)
{
    /*

        Print information to screen for the fisher matrices

    */

    /* Print nothing */
    if (!(print -> flags -> main))
        return 0;

    /* Allocate memory */
    char *str = malloc(boxSize);

    /* Started calculating the stage */
    if (print -> stageFinished == 0)
      {
        /* Start the timer */
        print -> startTimer = omp_get_wtime();

        /* Reset the execution time */
        print -> execTime = 0.;

        /* Print the fiducials */
        _print_fiducials(print);

        /* Print a header */
        _print_str("");
        _print_str("* * *");
        _print_str("*");

        if (print -> sizes != NULL)
          {
            snprintf(str, boxSize, "*  Calculating %s for (%ld, %ld, %ld) values of (z, k, mu) for a total of %ld values.",
                                    print -> id, print -> sizes[0], print -> sizes[1], print -> sizes[2], print -> sizes[3]);
            _print_str((const char*) str);

            _print_str("*");
            _print_str("*");
          }

        /* Stage is currently being calculated (no progress here) */
        print -> stageFinished = -1;
      }

    /* Currently calculating the stage */
    else if (print -> stageFinished == -1)
      {
        /* End the timer */
        print -> endTimer = omp_get_wtime();

        /* Execution time */
        print -> execTimeStep = print -> endTimer - print -> startTimer;
        print -> execTime += print -> execTimeStep;

        /* Reset the timer */
        print -> startTimer = print -> endTimer;

        /* Print the current progress */
        _print_progress(print);
      }

    /* Finished calculating the stage */
    else if (print -> stageFinished == 1)
      {
        /* End the timer */
        print -> endTimer = omp_get_wtime();

        /* Execution time */
        print -> execTimeStep = print -> endTimer - print -> startTimer;
        print -> execTime += print -> execTimeStep;

        /* Global execution time */
        _execTime += print -> execTime;

        /* Formatted global execution time */
        misc_time_t *timeGlobal = misc_time_format_out(_execTime);

        /* Print the current progress */
        _print_progress(print);

        /* Print finish message */
        _print_str("*");
        _print_str("*");

        snprintf(str, boxSize, "*  Finished calculating %s.",
                                print -> id);
        _print_str((const char*) str);

        _print_str("*");

        /* Total execution time */
        snprintf(str, boxSize, "*  Total Execution Time: %02uhr %02umin %02us %03ums",
                                timeGlobal -> hr, timeGlobal -> min, timeGlobal -> sec, timeGlobal -> msec);
        _print_str((const char*) str);

        /* Print a footer */
        _print_str("*");

        _print_str("* * *");

        _print_str("");
        _print_str("");

        /* Free memory */
        free(timeGlobal);
      }

    /* Free memory */
    free(str);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Print Progress   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _print_progress(print_t *print)
{
    /*

        Print the progress of the current stage.

    */

    /* Format the execution time */
    misc_time_format(print -> execTimeFormat, print -> execTime);

    if (print -> stageFinished == -1)
      {
        if (print -> progress < 0.005)
          {
            printf("*  ###   Current progress: %.1f%%  |  Elapsed time: %02uhr %02umin %02us %03ums  |  Remaining time: --hr --min --s ---ms   ###\r",
                    print -> progress * 100, print -> execTimeFormat -> hr, print -> execTimeFormat -> min, print -> execTimeFormat -> sec, print -> execTimeFormat -> msec);
          }

        else
          {
            /* Remaining time */
            misc_time_format(print -> remTimeFormat, (1. / print -> progress - 1.) * print -> execTime);

            printf("*  ###   Current progress: %.1f%%  |  Elapsed time: %02uhr %02umin %02us %03ums  |  Remaining time: %02uhr %02umin %02us %03ums   ###\r",
                    print -> progress * 100, print -> execTimeFormat -> hr, print -> execTimeFormat -> min, print -> execTimeFormat -> sec, print -> execTimeFormat -> msec,
                                             print -> remTimeFormat -> hr, print -> remTimeFormat -> min, print -> remTimeFormat -> sec, print -> remTimeFormat -> msec);
          }

        fflush(stdout);
      }

    else
      {
        char *str = malloc(boxSize);

        snprintf(str, boxSize, "*  ###   Current progress: %.1f%%  |  Elapsed time: %02uhr %02umin %02us %03ums   ###",
                                print -> progress * 100, print -> execTimeFormat -> hr, print -> execTimeFormat -> min, print -> execTimeFormat -> sec, print -> execTimeFormat -> msec);
        _print_str(str);

        free(str);
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------   Print Fiducials   --------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _print_fiducials(print_t *print)
{
    /*

        Print the fiducials for spectra modules

    */

    /* Print nothing */
    if (!print -> flags -> fid)
        return 0;


    /* Allocate memory */
    char *str = malloc(boxSize);


    /* Print header */

    _print_str("");

    snprintf(str, boxSize, "***  Fiducials (at redshift %.2f)  ***", print -> kern -> z);
    _print_str((const char*) str);

    _print_str("");


    /* Print the Lambda CDM fiducials */

    _print_str("LCDM:");
    _print_str("");


    snprintf(str, boxSize, "Omega_M0 = %.4f", print -> kern -> lcdm -> omegaM0);
    _print_str((const char*) str);

    snprintf(str, boxSize, "ns = %.4f", print -> kern -> lcdm -> ns);
    _print_str((const char*) str);

    snprintf(str, boxSize, "Growth Index = %.4f", print -> kern -> lcdm -> growthIndex);
    _print_str((const char*) str);


    _print_str("");
    _print_str("");


    /* Print the growth factor */

    snprintf(str, boxSize, "Linear growth factor:");
    _print_str((const char*) str);

    _print_str("");


    snprintf(str, boxSize, "Growth = %.4f", print -> kern -> growth);
    _print_str((const char*) str);


    _print_str("");
    _print_str("");


    /* Print the fiducials bootstrap parameters */

    snprintf(str, boxSize, "Bootstrap parameters:");
    _print_str((const char*) str);

    _print_str("");


    snprintf(str, boxSize, "a2Ga = %.4f", print -> kern -> btst -> a2Ga);
    _print_str((const char*) str);

    snprintf(str, boxSize, "d2Ga = %.4f", print -> kern -> btst -> d2Ga);
    _print_str((const char*) str);

    snprintf(str, boxSize, "a3GaA = %.4f", print -> kern -> btst -> a3GaA);
    _print_str((const char*) str);

    snprintf(str, boxSize, "a3GaB = %.4f", print -> kern -> btst -> a3GaB);
    _print_str((const char*) str);

    snprintf(str, boxSize, "d3GaA = %.4f", print -> kern -> btst -> d3GaA);
    _print_str((const char*) str);

    snprintf(str, boxSize, "d3GaB = %.4f", print -> kern -> btst -> d3GaB);
    _print_str((const char*) str);

    snprintf(str, boxSize, "h = %.4f", print -> kern -> btst -> h);
    _print_str((const char*) str);


    _print_str("");
    _print_str("");


    /* Print the fiducial RSD parameters */

    snprintf(str, boxSize, "Bias and RSD parameters:");
    _print_str((const char*) str);

    _print_str("");


    snprintf(str, boxSize, "b1 = %.4f", print -> kern -> bias -> b1);
    _print_str((const char*) str);

    snprintf(str, boxSize, "b2 = %.4f", print -> kern -> bias -> b2);
    _print_str((const char*) str);

    snprintf(str, boxSize, "bG2 = %.4f / c2Ga = %.4f", print -> kern -> bias -> bG2, print -> kern -> bias -> c2Ga);
    _print_str((const char*) str);

    snprintf(str, boxSize, "bGam3 = %.4f", print -> kern -> bias -> bGam3);
    _print_str((const char*) str);

    snprintf(str, boxSize, "sigs = sig0s (1+z)/H(z) = %.4f Mpc/h", print -> kern -> rsd -> sigs);
    _print_str((const char*) str);

    snprintf(str, boxSize, "sigv = %.4f Mpc/h", print -> kern -> rsd -> sigv);
    _print_str((const char*) str);


    _print_str("");
    _print_str("");


    /* Print the fiducial counter term parameters */

    _print_str("Counter term parameters:");

    _print_str("");


    snprintf(str, boxSize, "c0 = %.2f (Mpc/h)^2", print -> kern -> ctr -> c0);
    _print_str((const char*) str);

    snprintf(str, boxSize, "c2 = %.2f (Mpc/h)^2", print -> kern -> ctr -> c2);
    _print_str((const char*) str);

    snprintf(str, boxSize, "c4 = %.2f (Mpc/h)^4", print -> kern -> ctr -> c4);
    _print_str((const char*) str);

    _print_str("");
    _print_str("");


    /* Print the fiducial survey parameters */

    snprintf(str, boxSize, "Survey parameters:");
    _print_str((const char*) str);

    _print_str("");


    snprintf(str, boxSize, "V = %.2f (Gpc/h)^3", print -> kern -> surv -> V);
    _print_str((const char*) str);

    snprintf(str, boxSize, "n = %.2e (h/Mpc)^3", print -> kern -> surv -> n);
    _print_str((const char*) str);


    _print_str("");
    _print_str("");


    /* Free memory */
    free(str);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Print Threads   ---------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _print_spec_sub(print_t *print)
{
    /*

        Print what a thread is doing

    */

    /* Print nothing */
    if (!print -> flags -> sub)
        return 0;

    char *str = malloc(boxSize);

    char *strK = misc_scat(1, "(");
    char *strMu = misc_scat(1, "(");

    for (size_t i = 0; i < print -> kern -> specOrder; i++)
      {
        /* Copy the previous strings */
        char *strKCp = misc_scat(1, strK);
        char *strMuCp = misc_scat(1, strMu);
        free(strK);
        free(strMu);

        /* Double to string */
        char *strKVal = misc_dtos(print -> kern -> k[i], __XSBUFFER__);
        char *strMuVal = misc_dtos(print -> kern -> mu[i], __XSBUFFER__);

        /* Concatenate the strings */
        if (i < print -> kern -> specOrder - 1)
          {
            strK = misc_scat(3, strK, strKVal, ", ");
            strMu = misc_scat(3, strMu, strMuVal, ", ");
          }

        else
          {
            strK = misc_scat(3, strK, strKVal, ")");
            strMu = misc_scat(3, strMu, strMuVal, ")");
          }

        /* Free memory */
        free(strKCp);
        free(strMuCp);
        free(strKVal);
        free(strMuVal);
      }

    /* Thread started calculating something */
    if (print -> subFinished == 0)
      {
        snprintf(str, boxSize, "*  Thread %d is calculating '%s' for z = %.2f, k = %s and mu = %s!",
                                print -> thread, print -> id, print -> kern -> z, strK, strMu);
        _print_str((const char*) str);

        _print_str("*");
      }

    /* Thread finished calculating something */
    else
      {
        snprintf(str, boxSize, "*  Thread %d finished calculating '%s' for z = %.2f, k = %s and mu = %s!",
                                print -> thread, print -> id, print -> kern -> z, strK, strMu);
        _print_str((const char*) str);

        _print_str("*");
      }

    free(str);

    free(strK);
    free(strMu);

    return 0;
}


static int _print_cov_sub(print_t *print)
{
    /*

        Print what a thread is doing

    */

    /* Print nothing */
    if (!print -> flags -> sub)
        return 0;

    char *str = malloc(boxSize);

    char *strK = misc_scat(1, "(");
    char *strMu = misc_scat(1, "(");

    for (size_t i = 0; i < print -> kern -> specOrder; i++)
      {
        /* Copy the previous strings */
        char *strKCp = misc_scat(1, strK);
        char *strMuCp = misc_scat(1, strMu);
        free(strK);
        free(strMu);

        /* Double to string */
        char *strKVal = misc_dtos(print -> kern -> k[i], __XSBUFFER__);
        char *strMuVal = misc_dtos(print -> kern -> mu[i], __XSBUFFER__);

        /* Concatenate the strings */
        if (i < print -> kern -> specOrder - 1)
          {
            strK = misc_scat(3, strK, strKVal, ", ");
            strMu = misc_scat(3, strMu, strMuVal, ", ");
          }

        else
          {
            strK = misc_scat(3, strK, strKVal, ")");
            strMu = misc_scat(3, strMu, strMuVal, ")");
          }

        /* Free memory */
        free(strKCp);
        free(strMuCp);
        free(strKVal);
        free(strMuVal);
      }

    /* Thread started calculating something */
    if (print -> subFinished == 0)
      {
        snprintf(str, boxSize, "*  Thread %d is calculating '%s' for z = %.2f, k = %s and mu = %s!",
                                print -> thread, print -> id, print -> kern -> z, strK, strMu);
        _print_str((const char*) str);

        _print_str("*");
      }

    /* Thread finished calculating something */
    else
      {
        snprintf(str, boxSize, "*  Thread %d finished calculating '%s' for z = %.2f, k = %s and mu = %s!",
                                print -> thread, print -> id, print -> kern -> z, strK, strMu);
        _print_str((const char*) str);

        _print_str("*");
      }

    free(str);

    free(strK);
    free(strMu);

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     MISC.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "misc.h"




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     GENERAL     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


int misc_stream_skip(FILE *stream, size_t lines)
{
    /*

        Skip 'lines' number of lines in a stream

    */

    char c;
    size_t skippedLines = 0;

    do
      {
        c = (char) fgetc(stream);

        if (c == '\n')
            skippedLines++;
      }

    while (skippedLines < lines && c != EOF);

    if (c == EOF)
      {
        printf("Reached the end of the file before being able to skip %ld lines.\n", lines);
        exit(1);

        return 1;
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Concatenate Strings   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *misc_scat(size_t size, ...)
{
    /*

        Concatenate several strings

    */

    char *strCat = malloc(1);
    *strCat = '\0';

    va_list ap;

    va_start(ap, size);

    for (size_t i = 0; i < size; i++)
      {
        char *str = va_arg(ap, char*);

        /* Skip strings that have no memory */
        if (str == NULL)
            continue;

        size_t strLen = strlen(strCat) + strlen(str) + 1;

        strCat = realloc(strCat, strLen);
        strncat(strCat, str, strLen);
      }

    va_end(ap);

    /* String did not contain any memory */
    if (strCat[0] == '\0')
      {
        free(strCat);
        strCat = NULL;
      }

    return strCat;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -------------------------------------------   Copy Strings   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int misc_scp(char **str1, const char *str2)
{
    /*

        Copy str2 to str1

    */

    if (str2 == NULL || strlen(str2) == 0)
      {
        free(*str1);
        *str1 = NULL;

        return 0;
      }

    *str1 = realloc(*str1, strlen(str2) + 1);
    snprintf(*str1, strlen(str2) + 1, "%s", str2);

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------   Count Characters in String   ----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t misc_scnch(const char *str, char ch)
{
    /*

        Count the number of times ch appears in str

    */

    if (str == NULL || strlen(str) == 0)
        return 0;

    size_t count = 0;

    for (size_t i = 0; i < strlen(str); i++)
      {
        if (str[i] == ch)
            count++;
      }

    return count;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   Search for String within String   --------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool misc_sin_rm(char *str1, const char *str2)
{
    /*

        Check whether str2 is in str1 and return str1 stripped off of the first appearance of str2

    */

    /* Initalise sucess */
    bool success = false;

    /* Lengths of the strings */
    size_t lenStr1 = strlen(str1);
    size_t lenStr2 = strlen(str2);

    /* Cannot have str2 longer than str1 */
    if (lenStr1 < lenStr2)
      {
        return success;
      }

    /* Check if str2 is in str1 and get the starting and ending position inside str1 */
    char c1;
    char c2;

    size_t i = 0;
    size_t j = 0;

    size_t start = 0;
    size_t end = 0;

    do
      {
        c1 = str1[i];
        c2 = str2[j];

        /* Characters match */
        if (c1 == c2)
          {
            if (j == 0)
                start = i;

            j++;
          }

        else
            j = 0;

        /* Found str2 in str1 */
        if (j == lenStr2)
          {
            end = i;

            break;
          }

        i++;
      }

    while (c1 != '\0' && c2 != '\0');

    /* Remove str2 from str1 (in strCp) */
    if (end != 0)
      {
        success = true;

        for (i = start; i < lenStr1 - (end + 1) + start + 1; i++)
          {
            str1[i] = str1[i + (end + 1) - start];
          }
      }

    return success;
}


/*  ------------------------------------------------------------------------------------------------------  */


bool misc_sin(const char *str1, const char *str2)
{
    /*

        Check whether str2 is in str1 (difference to misc_sin_rm: does not change str1)

    */

    size_t lenStr1 = strlen(str1);
    size_t lenStr2 = strlen(str2);

    if (lenStr1 < lenStr2)
        return false;

    char c1;
    char c2;

    size_t i = 0;
    size_t j = 0;

    do
      {
        c1 = str1[i];
        c2 = str2[j];

        if (c1 == c2)
            j++;

        else
            j = 0;

        i++;
      }

    while (c1 != '\0' && c2 != '\0' && j != lenStr2);

    if (j == lenStr2)
        return true;

    return false;
}


bool misc_sins(char **array, size_t size, const char *str, size_t *index)
{
    /*

        Check if the string 'str' is part of a string of an array of strings and overwrite index with the earliest hit

    */

    for (size_t i = 0; i < size; i++)
      {
        if (misc_sin((const char*) array[i], str))
          {
            if (index != NULL)
                *index = i;

            return true;
          }
      }

    return false;
}


/*  ------------------------------------------------------------------------------------------------------  */


bool misc_sinci(const char *str1, const char *str2)
{
    /*

        Check whether str2 is in str1 irrespective of upper or lower case - S(tring) IN C(ase) I(ndependent)

    */

    size_t lenStr1 = strlen(str1);
    size_t lenStr2 = strlen(str2);

    if (lenStr1 < lenStr2)
        return false;

    char c1;
    char c2;

    size_t i = 0;
    size_t j = 0;

    do
      {
        c1 = str1[i];
        c2 = str2[j];

        if (tolower(c1) == tolower(c2))
            j++;

        else
            j = 0;

        i++;
      }

    while (c1 != '\0' && c2 != '\0' && j != lenStr2);

    if (j == lenStr2)
        return true;

    return false;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *misc_ssrm(char *str, const char *sub)
{
    /*

        Remove a sub string from a string

        Credit: "https://stackoverflow.com/questions/47116974/remove-a-substring-from-a-string-in-c"

    */

    size_t len = strlen(sub);

    if (len > 0)
      {
        char *p = str;

        while ((p = strstr(p, sub)) != NULL)
          {
            memmove(p, p + len, strlen(p + len) + 1);
          }
      }

    return str;
}


bool misc_srm(char **array, size_t *size, const char *str)
{
    /*

        Remove a string from an array of strings. Returns false if the string does not exist in the array.

    */

    bool success;

    /* Get the index of the string in the array */
    size_t index = misc_ssearch_index(array, *size, str, &success);

    /* Not in the array */
    if (!success)
        return false;

    /* Decrease the size of the array */
    *size -= 1;

    /* Push the string to the end of the array */
    for (size_t i = index; i < *size; i++)
      {
        misc_swap(&array[i], &array[i+1], "s");
      }

    array = realloc(array, sizeof(char*) * *size);

    return true;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Search for String in Array   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


bool misc_ssearch(char **array, size_t size, const char *str, size_t *index)
{
    /*

        Search an array of strings for the string str and return the result of the search

    */

    /* No strings in the array */
    if (array == NULL)
      {
        if (index != NULL) *index = size;

        return false;
      }

    /* Loop through the elements and break if the string has been found */
    for (size_t i = 0; i < size; i++)
      {
        if (!strcmp(array[i], str))
          {
            if (index != NULL) *index = i;

            return true;
          }
      }

    return false;
}


size_t misc_ssearch_index(char **array, size_t size, const char *str, bool *success)
{
    /*

        Search an array of strings for the earliest encounter of the string str. If str is not in the array set success to false
        and return the length of the array

    */

    /* Initialise success in case the string does not exist in the array */
    if (success != NULL) *success = false;

    /* No strings in the array */
    if (array == NULL)
        return 0;

    /* Initialise the index in case the string does not exist in the array */
    size_t index = size;

    /* Loop through the elements and break if the string has been found */
    for (size_t i = 0; i < size; i++)
      {
        if (!strcmp(array[i], str))
          {
            if (success != NULL) *success = true;

            index = i;
            break;
          }
      }

    return index;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------   Turn Double to String And Back   ---------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


char *misc_dtos(double num, unsigned int tol)
{
    /*

        Turn a double into a string in scientifc notation of maximum tolerance tol

    */

    /* Final precision */
    int prec = (int) tol;

    /* Get the full string */
    char *str = malloc((size_t) tol + __SBUFFER__);
    snprintf(str, (size_t) tol + __SBUFFER__, "%.*g", prec, num);
    size_t len = strlen(str);

    /* Reduce the precision if there are zeros at the end of the string */
    bool reduce = false;

    for (size_t i = 0; i < len; i++)
      {
        if (reduce)
          {
            /* Decrease precision */
            if (str[len - i - 1] == '0')
                prec--;

            /* Break the for loop */
            else
                break;
          }

        else
          {
            /* Can start reducing the precision if an 'e' or 'E' has been hit */
            if (str[len - i - 1] == 'e' || str[len - i - 1] == 'E')
                reduce = true;
          }
      }

    /* Get the string with the final precision and reallocate memory */
    snprintf(str, (size_t) tol + __SBUFFER__, "%.*g", prec, num);
    str = realloc(str, strlen(str) + 1);

    return str;
}


double misc_stod(const char *str, char **remstr, bool *success)
{
    /*

        Extract the first number from str and return it. If there is no number return 0. and set success to false

    */

    char c;
    size_t ind = 0;

    char *strF = malloc(strlen(str) + 1);
    size_t indF = 0;

    bool startF = false;

    do
      {
        c = str[ind];
        ind++;

        /* Char is an int from 0 to 9 */
        if (c - '0' >= 0 && c - '0' <= 9)
            startF = true;

        /* Store the string of the double in strF */
        if (startF)
          {
            strF[indF] = c;
            indF++;
          }
      }

    while (c != '\0');

    /* Null terminate the string */
    strF[indF] = '\0';

    /* Get the double value (and remainder of the string) */
    double res;
    char *rem;
    res = strtod(strF, &rem);

    if (remstr != NULL)
      {
        misc_scp(remstr, rem);
      }

    free(strF);

    /* Set the success */
    if (indF == 0 && success != NULL)
        *success = false;

    else if (indF != 0 && success != NULL)
        *success = true;

    return res;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Binary Search   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _bsearch_bounds(double *array, int *low, int *high, double value, double tol);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


bool misc_bsearch(double *array, size_t size, double value, double tol, size_t *index)
{
    /*

        Perform a binary search for an element (value) in an array (sorted small to large) of given size while
        allowing an absolute tolerance on value.

        This function returns whether the element was found or not. If it was found, the index is stored in (index),
        otherwise the index at which the element would have been expected given its value.

    */

    /* No memory allocated */
    if (array == NULL)
      {
        if (index != NULL) *index = 0;

        return false;
      }

    /* Bounds of binary search */
    int low = 0;
    int high = (int) size - 1;

    /* Success of binary search */
    bool success;

    do
      {
        /* Check the new bounds */
        success = _bsearch_bounds(array, &low, &high, value, tol);
      }

    while (!success && low <= high);

    /* Index of element is low */
    if (index != NULL) *index = (size_t) low;

    return success;
}


/*  ------------------------------------------------------------------------------------------------------  */


size_t misc_bsearch_index(double *array, size_t size, double value, double tol, bool *success)
{
    /*

        Perform a binary search for an element (value) in an array (sorted small to large) of given size
        while allowing an absolute tolerance on the value.

        This function returns the index of the found element and sets success to true, or if not found
        returns the index at which the element would have been expected given its value with success being
        set to false.

    */

    /* No memory allocated */
    if (array == NULL)
      {
        if (success != NULL) *success = false;

        return 0;
      }

    /* Index of element in array */
    size_t index;

    /* Bounds of binary search */
    int low = 0;
    int high = (int) size - 1;

    /* Success of binary search */
    bool success_;

    do
      {
        /* Check the new bounds */
        success_ = _bsearch_bounds(array, &low, &high, value, tol);
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


static bool _bsearch_bounds(double *array, int *low, int *high, double value, double tol)
{
    /*

        Look up the mid point in an array between positions low and high and check if value is equal, higher or lower
        and change low and high accordingly.

    */

    /* Index of the mid point */
    int mid = (*high - *low) / 2 + *low;

    /* Difference between mid point and value */
    double diff = array[mid] - value;

    /* Success on finding the value */
    bool success = (fabs(diff) <= tol) ? true : false;

    /* Update high and low */
    *high = (fabs(diff) <= tol) ? mid : (diff > 0.) ? mid - 1 : *high;
    *low = (fabs(diff) <= tol) ? mid : (diff < 0.) ? mid + 1 : *low;

    return success;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int misc_insort(double **array, size_t *size, double value, double tol, bool *success)
{
    /*

        Given an ordered array of some size, insert an element (value) such that the array will still be ordered
        while allowing an absolute tolerance on the value.

        Set success to true if the value was inserted or false if it already exists in the array

    */

    /* Success variable */
    bool success_;

    /* Get the index of the value */
    size_t index = misc_bsearch_index(*array, *size, value, tol, &success_);

    /* Value is already in the array */
    if (success_)
      {
        if (success != NULL) *success = false;

        return 0;
      }


    /* Must insert the value */

    /* Reallocate memory */
    *array = realloc(*array, sizeof(double) * (*size + 1));

    for (int i = (int) *size; i > (int) index; i--)
      {
        /* Push the i'th element forwards */
        (*array)[i] = (*array)[i - 1];
      }

    /* Insert value */
    (*array)[index] = value;

    /* Increment the size */
    *size += 1;

    /* Set success to true */
    if (success != NULL) *success = true;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------   Swap Pointers   -----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int misc_swap(void *a, void *b, const char *type)
{
    /*

        Swap two elements

    */


    if (!strcmp(type, "f"))
      {
        double *x = (double*) a;
        double *y = (double*) b;

        double buf = *x;

        *x = *y;
        *y = buf;
      }

    else if (!strcmp(type, "*f"))
      {
        double **x = (double**) a;
        double **y = (double**) b;

        double *buf = *x;

        *x = *y;
        *y = buf;
      }

    else if (!strcmp(type, "ld"))
      {
        size_t *x = (size_t*) a;
        size_t *y = (size_t*) b;

        size_t buf = *x;

        *x = *y;
        *y = buf;
      }

    else if (!strcmp(type, "*ld"))
      {
        size_t **x = (size_t**) a;
        size_t **y = (size_t**) b;

        size_t *buf = *x;

        *x = *y;
        *y = buf;
      }

    else if (!strcmp(type, "b"))
      {
        bool *x = (bool*) a;
        bool *y = (bool*) b;

        bool buf = *x;

        *x = *y;
        *y = buf;
      }

    else if (!strcmp(type, "s"))
      {
        char **x = (char**) a;
        char **y = (char**) b;

        char *buf = *x;

        *x = *y;
        *y = buf;
      }

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------------   Sort Array   ------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static bool _sort_comp(double a, double b);

static int _sort_kern(double *array, size_t *ind, int low, int high, bool (*comp)(double, double));
static int _sort_partition(double *array, size_t *ind, int low, int high, bool (*comp)(double, double));

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t *misc_sort(double *array, size_t size)
{
    /*

        Sort an array from small to large

    */

    /* Keep track of the swapped indices */
    size_t *ind = malloc(sizeof(size_t) * size);

    for (size_t i = 0; i < size; i++)
        ind[i] = i;

    /* Sort the array */
    _sort_kern(array, ind, 0, (int) size - 1, _sort_comp);

    return ind;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static bool _sort_comp(double a, double b)
{
    /*

        Compare two elements and return the sign of a - b (Small To Large)

    */

    return (a - b) <= 0.;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _sort_partition(double *array, size_t *ind, int low, int high, bool (*comp)(double, double))
{
    /*

        Partition the array

    */

    /* Push i'th value to the right */
    int i = low - 1;

    for (int j = low; j < high; j++)
      {
        /* Push j'th value to the left if comp is true */
        if (comp(array[j], array[high]))
          {
            misc_swap(&array[i + 1], &array[j], "f");
            misc_swap(&ind[i + 1], &ind[j], "ld");

            i++;
          }
      }

    /* Swap i'th value with pivot's value */
    misc_swap(&array[i + 1], &array[high], "f");
    misc_swap(&ind[i + 1], &ind[high], "ld");

    return i + 1;
}


static int _sort_kern(double *array, size_t *ind, int low, int high, bool (*comp)(double, double))
{
    /*

        Kernel of the quicksort function

    */

    /* Sort as long as pivot is larger than low */
    if (low < high)
      {
        /* Place last pivot, and get new pivot */
        int pivot = _sort_partition(array, ind, low, high, comp);

        /* Sort to the left of the previous pivot */
        _sort_kern(array, ind, low, pivot - 1, comp);

        /* Sort to the right of the last pivot */
        _sort_kern(array, ind, pivot + 1, high, comp);
      }

    return 0;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     RANDOM NUMBERS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */



size_t *misc_random_sample(size_t start, size_t stop, size_t *size, bool repetition)
{
    /*

        Draw a random uniform integer sample from [start, stop] 'size' number of times. If repetition is false,
        only include unique values. If size < stop - start + 1, size is updated and all integer values
        in [start, stop] are returned.

    */

    size_t *sample = NULL;

    /* With repetition */
    if (repetition)
      {
        /* Allocate memory */
        sample = malloc(sizeof(size_t) * *size);

        for (size_t i = 0; i < *size; i++)
          {
            sample[i] = gsl_rng_uniform_int(_randomNumberGen_, stop - start + 1) + start;
          }
      }

    /* Without repetition */
    else
      {
        /* Maximal possible sample size for unique numbers in sample */
        size_t sampleSizeMax = stop - start + 1;

        /* Requested size is too large -> store only unique numbers */
        if (*size >= sampleSizeMax)
          {
            /* Update size */
            *size = sampleSizeMax;

            /* Allocate memory */
            sample = malloc(sizeof(size_t) * *size);

            for (size_t i = start, index = 0; i < stop + 1; i++)
                sample[index++] = i;
          }

        else
          {
            /* Allocate memory */
            sample = malloc(sizeof(size_t) * *size);

            /* Only need to sample max(size, sampleSizeMax - size); In he latter case we sample rejected number in [start, stop] */
            size_t sampleSize = (*size > sampleSizeMax - *size) ? sampleSizeMax : *size;

            /* Keep track 0f visited numbers */
            size_t visitedNumbersSize = 0;
            double *visitedNumbers = NULL;
            bool successInsert;

            for (size_t i = 0; i < sampleSize; i++)
              {
                do
                  {
                    /* Random number */
                    sample[i] = gsl_rng_uniform_int(_randomNumberGen_, stop - start + 1) + start;

                    /* Search for and insert the index i (as double) in visitedNumbers */
                    misc_insort(&visitedNumbers, &visitedNumbersSize, (double) sample[i], __ABSTOL__, &successInsert);
                  }

                while (!successInsert); // Sample random numbers until unique was found
              }

            /* Get the accepted numbers in the case where rejected numbers were sampled */
            if (*size != sampleSize)
              {
                for (size_t i = start, index = 0; i < stop + 1; i++)
                  {
                    /* Skip rejected numbers */
                    if ((size_t) visitedNumbers[i] == i) continue;

                    /* Store accepted numbers */
                    sample[index++] = i;
                  }
              }

            /* Free memory */
            free(visitedNumbers);
          }
      }

    return sample;
}





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COMBINATORICS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


size_t misc_factorial(size_t n)
{
    /*

        Calculate n! = 1 * ... * n

    */

    /* Default case */
    if (n == 0)
        return 1;

    /* Calculate the factorial of n */
    size_t factorial = 1;

    for (size_t i = 1; i < n; i++)
      {
        factorial *= (i + 1);
      }

    return factorial;
}


size_t misc_binom(size_t n, size_t m)
{
    /*

        Caluclate the binomial coefficient (n choose m)

    */

    /* Error */
    if (n < m)
      {
        printf("Cannot compute (%ld choose %ld).\n", n, m);
        exit(1);

        return 1;
      }

    /* Must only perform min(m, n-m) multiplications */
    size_t k = (m > n - m) ? n - m : m;

    /* Default case */
    if (k == 0)
        return 1;

    /* Calculate the binomial coefficient */
    size_t binom = 1;

    for (size_t i = 0; i < k; i++)
      {
        binom *= n - i;
      }

    binom /= misc_factorial(k);

    return binom;
}


size_t misc_multinom(size_t n, size_t *m, size_t size)
{
    /*

        Calculate the multinomial coefficient

            (n choose m1, ..., mk) = n! / (m1! * ... * mk!) = exp( log Γ(n + 1) - log Γ(m1 + 1) - ... - log Γ(mk + 1) )

    */

    double binom = gamma((double) (n + 1));

    for (size_t i = 0; i < size; i++)
      {
        binom -= gamma((double) (m[i] + 1));
      }

    return (size_t) round(exp(binom));
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     TIME     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


misc_time_t *misc_time_new()
{
    /*

        New misc_time_t struct

    */

    misc_time_t *time = malloc(sizeof(misc_time_t));

    time -> msec = 0.;
    time -> sec = 0.;
    time -> min = 0.;
    time -> hr = 0.;

    return time;
}


misc_time_t *misc_time_free(misc_time_t *time)
{
    /*

        Free a misc_time_t struct

    */

    /* Check for NULL */
    if (time == NULL)
        return NULL;

    /* Free time */
    free(time);

    return NULL;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int misc_time_format(misc_time_t *time, double timeSec)
{
    /*

        Given a time in seconds transform it into hours, minutes, seconds and milliseconds

    */

    unsigned int timeSecFloor = (unsigned int) floor(timeSec);

    time -> msec = (unsigned int) floor((timeSec - timeSecFloor) * 1000);

    time -> hr = (timeSecFloor - timeSecFloor % 3600) / 3600;
    timeSecFloor = timeSecFloor % 3600;

    time -> min = (timeSecFloor - timeSecFloor % 60) / 60;
    timeSecFloor = timeSecFloor % 60;

    time -> sec = timeSecFloor;

    return 0;
}


misc_time_t *misc_time_format_out(double timeSec)
{
    /*

        Given a time in seconds transform it into hours, minutes, seconds and milliseconds

    */

    misc_time_t *time = malloc(sizeof(misc_time_t));

    unsigned int timeSecFloor = (unsigned int) floor(timeSec);

    time -> msec = (unsigned int) floor((timeSec - timeSecFloor) * 1000);

    time -> hr = (timeSecFloor - timeSecFloor % 3600) / 3600;
    timeSecFloor = timeSecFloor % 3600;

    time -> min = (timeSecFloor - timeSecFloor % 60) / 60;
    timeSecFloor = timeSecFloor % 60;

    time -> sec = timeSecFloor;

    return time;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

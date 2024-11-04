/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     INPUT.C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


#include "kernels.h"





/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     KERNEL STRUCT     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------    (Internal) Kern Struct    ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static kern_t *_kernels_work_new(size_t kernOrder)
{
    /*

        Create a new kern_t struct used by the kernels_... functions that require all combinations of the ki_ at given kernel order:

        Have (n = kernOrder)

                2^n - 1 = sum_k=1^n (n choose k)

        values for k, mu and

                (3^n - (2^(n+1) - 1)) / 2 = 1/2 sum_k1,k2=1^n (n choose k1,k2)

        values for nu.

        Note: Given k1_, ..., kn_ we require all sums ki1_ + ... + kij_ which result in the equations above.

    */

    kern_t *kern = malloc(sizeof(kern_t));

    kern -> kernOrder = kernOrder;

    kern -> k = calloc((size_t) ((1 << kernOrder) - 1), sizeof(double));
    kern -> mu = calloc((size_t) ((1 << kernOrder) - 1), sizeof(double));
    kern -> nu = calloc(((size_t) pow(3., (double) kernOrder) - (size_t) ((1 << (kernOrder + 1)) - 1)) / 2, sizeof(double));

    return kern;
}


static kern_t *_kernels_work_free(kern_t *kern)
{
    /*

        Free kern used by the kernels_... functions

        Note: Can never free fiducials !!! Only consider them as pointers to the original fiducial struct

    */

    /* Check for NULL */
    if (kern == NULL)
        return NULL;

    free(kern -> k);
    free(kern -> nu);
    free(kern -> mu);

    /* Free kern itself */
    free(kern);

    /* Set kern to NULL so calling free again is no problem */
    kern = NULL;

    return kern;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ####################################   Function Declarations   #######################################  */

static int _kernels_work_var_k(kern_t *kernWork, size_t **indices, size_t *size, size_t **offset, size_t *depth);
static int _kernels_work_var_nu(kern_t *kernWork, size_t **indices, size_t *size, size_t **offset, size_t *depth);

static int _kernels_work_var_k_cp(kern_t *kernWorkIn, kern_t *kernWorkOut, size_t **indices, size_t *size, size_t **offset, size_t *depth);
static int _kernels_work_var_nu_cp(kern_t *kernWorkIn, kern_t *kernWorkOut, size_t **indices, size_t *size, size_t **offset, size_t *depth);

static size_t _kernels_work_index_k(size_t kernOrder, size_t *indices, size_t size);
static size_t _kernels_work_index_nu(size_t kernOrder, size_t **indices, size_t *size);

/*  ######################################################################################################  */
/*  ------------------------------------------------------------------------------------------------------  */


static int _kernels_work_var(kern_t *kernVar, size_t kernOrder)
{
    /*

        Get the working kern_t structs up to kernOrder for the variables in kernVar

        Assumes that kernOrder <= kernVar -> kernOrder!

    */


    /* Working kern_t structs */
    kern_t **kernWork = (kern_t**) kernVar -> kernWork;

    /* Recursion depth */
    size_t depth = 0;

    /* Indices */
    size_t **indices = malloc(sizeof(size_t) * 2);
    indices[0] = malloc(sizeof(size_t) * kernVar -> kernOrder);
    indices[1] = malloc(sizeof(size_t) * kernVar -> kernOrder);

    /* Offset */
    size_t **offset = malloc(sizeof(size_t*) * 2);

    offset[0] = malloc(sizeof(size_t) * 3);
    offset[1] = malloc(sizeof(size_t) * 3);

    /* Fill working kern_t structs */
    for (int n = (int) (kernOrder - 1); n > -1; n--)
      {
        /* Must-have fiducials */
        kernWork[n] -> btst = kernVar -> btst;
        kernWork[n] -> bias = kernVar -> bias;
        kernWork[n] -> rsd = kernVar -> rsd;

        /* Initialise offset */
        offset[0][0] = 0;
        offset[1][0] = 0;

        /* Calculate all variables for the highest order kern_t struct */
        if ((size_t) n == kernOrder - 1)
          {
            /* ki, mui, nuij */
            for (size_t i = 0; i < kernWork[n] -> kernOrder; i++)
              {
                /* ki + mui */
                kernWork[n] -> k[offset[0][0]] = kernels_qget_k(kernVar, i);
                kernWork[n] -> mu[offset[0][0]] = kernels_qget_mu(kernVar, i);

                /* Incremenet offset[0][0] */
                offset[0][0]++;

                for (size_t j = i + 1; j < kernWork[n] -> kernOrder; j++)
                  {
                    /* nuij */
                    kernWork[n] -> nu[offset[1][0]] = kernels_qget_nu(kernVar, i, j);

                    /* Incremenet offset[1][0] */
                    offset[1][0]++;
                  }
              }

            /* ki1...im_ = ki1_ + ... + kim1_ for m1 > 1 */
            for (size_t m1 = 2; m1 <= kernWork[n] -> kernOrder; m1++)
              {
                /* First k..._ vector has m1 added wavevectors */
                size_t size[2] = {m1, 1};

                /* Calculate all possible k, mu, and nu */
                _kernels_work_var_k(kernWork[n], indices, size, offset, &depth);
              }
          }

        /* Copy all variables from the highest order kern_t struct */
        else
          {
            /* ki1...im_ = ki1_ + ... + kim1_ */
            for (size_t m1 = 1; m1 <= kernWork[n] -> kernOrder; m1++)
              {
                /* First k..._ vector has m1 added wavevectors */
                size_t size[2] = {m1, 1};

                /* Copy all possible k, mu, and nu */
                _kernels_work_var_k_cp(kernWork[kernOrder - 1], kernWork[n], indices, size, offset, &depth);
              }
          }
      }

    /* No longer need to compute work kern_t structs */
    kernVar -> computeWork = false;

    /* Free memory */
    free(indices[0]);
    free(indices[1]);
    free(offset[0]);
    free(offset[1]);

    free(indices);
    free(offset);


    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _kernels_work_var_k(kern_t *kernWork, size_t **indices, size_t *size, size_t **offset, size_t *depth)
{
    /*

        Calculate all possible

            k = |ki1_ + ... + kim1_|

            mu = (ki1_.e_ + ... + kim1_.e_) / k = (ki1 * mui1 + ... + kim1 * muim1) / k

            nu = (ki1_ + ... + kim1_).(kj1_ + ... + kjm2_) / (|ki1_ + ... + kim1_| * |kj1_ + ... + kjm2_|).


        The array "indices" must be initialised as

            size_t **indices = malloc(sizeof(size_t*) * 2);

            indices[0] = malloc(sizeof(size_t) * kernWork -> kernOrder);
            indices[1] = malloc(sizeof(size_t) * kernWork -> kernOrder);

        The array "size" must be initialised as

            size_t *size = malloc(sizeof(size_t) * 2);

            size[0] = m1;

        for m1 as given above.

        The array "offset" must be initialised as

            size_t **offset = malloc(sizeof(size_t*) * 2);

            offset[0] = calloc(3, sizeof(size_t));
            offset[1] = calloc(3, sizeof(size_t));

        The pointer "depth" must be initialised as

            size_t *depth = calloc(1, sizeof(size_t));

        or passed to _kernels_work_var_k as

            size_t depth = 0;

            _kernels_work_var_k(kernWork, indices, size, offset, &depth);

    */

    /* Nothing to do */
    if (*size == 1)
        return 0;

    /* Loop through all possible i1, ..., im1 */
    if (*depth < *size)
      {
        /* Initial index for i ; Must start at 0 for i1, and at in-1 + 1 for in with n = *depth */
        size_t ini = (*depth == 0) ? 0 : indices[0][*depth - 1] + 1;

        for (size_t i = ini; i < kernWork -> kernOrder; i++)
          {
            /* Index at location *depth in indices array is i */
            indices[0][*depth] = i;

            /* Get the next index at *depth + 1 */
            *depth += 1;
            _kernels_work_var_k(kernWork, indices, size, offset, depth);

            /* Reset *depth to adjust the index in indices array */
            *depth -= 1;
          }
      }

    /* Calculate ki1_ + ... + kim_ */
    else
      {
        /**

            Calculate k_ = k1_ + k2_, where k1_ = ki1_ + ... + kim1-1_ and k2_ = kim1_

        **/


        /**  Index of k1_ and k2_  **/

        /* Reduce size[0] for k1_ */
        *size -= 1;

        /* Index of k1_ */
        offset[0][1] = _kernels_work_index_k(kernWork -> kernOrder, *indices, *size);

        /* Index of k2_ */
        offset[0][2] = indices[0][*size];


        /**  Index of nu12  **/

        /* Reduce size[1] for k2_ = kim1_ */
        size_t bufSize = size[1];
        size[1] = 1;

        /* Let indices[1] point to index of k2 */
        size_t *buf = indices[1];
        indices[1] = &offset[0][2];

        /* Index of nu12 */
        size_t indexNu = _kernels_work_index_nu(kernWork -> kernOrder, indices, size);

        /* Reset size[1] and indices[1] */
        size[1] = bufSize;
        indices[1] = buf;

        /* Reset size[0] */
        *size += 1;


        /**  Calculate k = |k_1 + k2_| = sqrt( k1^2 + k2^2 + 2 * k1 * k2 * nu12 ) and mu = (mu1 * k1 + mu2 * k2) / k  **/

        /* Read k1, k2, nu12 */
        double k1 = kernWork -> k[offset[0][1]];
        double k2 = kernWork -> k[offset[0][2]];
        double nu12 = kernWork -> nu[indexNu];

        /* Calculate k and mu ; If k = 0, mu is undefined -> mu = 0 to avoid issues involving multiplications of mu with finite quantities */
        kernWork -> k[offset[0][0]] = (fabs(k1 - k2) > __ABSTOL__ || fabs(nu12 + 1) > __ABSTOL__) ? sqrt( fabs(k1*k1 + k2*k2 + 2. * k1 * k2 * nu12) ) : __ABSTOL__;
        kernWork -> mu[offset[0][0]] = (kernWork -> k[offset[0][0]] > __ABSTOL__) ? (k1 * kernWork -> mu[offset[0][1]] + k2 * kernWork -> mu[offset[0][2]]) / kernWork -> k[offset[0][0]] : 0.;


        /**

            Calculate all possible (ki1_ + ... + kim1_).(kj1_ + ... + kjm2_) / (|ki1_ + ... + kim1_| * |kj1_ + ... + kjm2_|)

        **/

        /* Set depth to 0 for nu functions */
        *depth = 0;

        /* Restrict m2 to m2 <= m1 */
        for (size[1] = 1; size[1] <= size[0]; size[1]++)
          {
            _kernels_work_var_nu(kernWork, indices, size, offset, depth);
          }

        /* Reset depth */
        *depth = *size;

        /* Incrmeent offset[0][0] */
        offset[0][0] += 1;
      }

    return 0;
}


static int _kernels_work_var_nu(kern_t *kernWork, size_t **indices, size_t *size, size_t **offset, size_t *depth)
{
    /*

        Calculate all possible

            nu = (ki1_ + ... + kim1_).(kj1_ + ... + kjm2_) / (|ki1_ + ... + kim1_| * |kj1_ + ... + kjm2_|)

    */

    /* Loop through all possible j1, ..., jm1 */
    if (*depth < size[1])
      {
        /* Keep track of the i1, ..., im1 ; indices[1][*depth] holds the value of i from the previous recursion */
        size_t i = (*depth == 0) ? ( (size[0] == size[1]) ? 1 : 0) : indices[1][*depth];

        /* Initial index for j */
        size_t ini = (*depth == 0) ? ( (size[0] == size[1]) ? indices[0][0] + 1 : 0) : indices[1][*depth - 1] + 1;

        for (size_t j = ini; j < kernWork -> kernOrder; j++)
          {
            /* Cannot have indices that are occupied by i1, ..., im1 */
            if (i < size[0] && j == indices[0][i])
              {
                i++;
                continue;
              }

            /* Index at location *depth in indices array is j */
            indices[1][*depth] = j;

            /* Store the current value of i in indices[1][*depth + 1] to initialise i in the next recursion */
            if (*depth < size[1])
                indices[1][*depth + 1] = i;

            /* Get the next index at *depth + 1 */
            *depth += 1;
            _kernels_work_var_nu(kernWork, indices, size, offset, depth);

            /* Reset *depth */
            *depth -= 1;
          }
      }

    else
      {
        /**
            Calculate

                nu = (k1_ + k2_).k3_ / (|k1_ + k2_| k3) = (nu13 * k1 + nu23 * k2) / |k1_ + k3_|

            where k1_ = ki1_ + ... + kim1-1_ and k2_ = kim1_, as well as k3_ = kj1_ + ... + kjm2_ if m1 > m2. If m1 = m2, we must calculate nu as

                nu = k1_.(k2_ + k3_) / (k1 |k2_ + k3_|) = (nu12 * k2 + nu13 * k3) / |k2_ + k3_|,

            where k1_ = ki1_ + ... + kim1_, as well as k2_ = kj1_ + ... + kjm2-1_ and k3_ = kjm2_

        **/

        /* If k = |ki1_ + ... + kim1_| = 0, we get issues for nu (see below) -> nu = 0 to avoid issues involving multiplications of nu with finite quantities */
        if (kernWork -> k[offset[0][0]] <= __ABSTOL__)
          {
            kernWork -> nu[offset[1][0]] = 0.;
            offset[1][0] += 1;

            return 0;
          }

        /* m1 > m2 */
        if (size[0] > size[1])
          {
            /**  Index of nu13  **/

            /* Reduce size[0] for k1_ */
            *size -= 1;

            /* Index of nu13 */
            offset[1][1] = _kernels_work_index_nu(kernWork -> kernOrder, indices, size);

            /* Reset size[0] */
            *size += 1;


            /**  Index of nu23  **/

            /* Reduce size[0] for k2_ = kim1_ */
            size_t bufSize = size[0];
            size[0] = 1;

            /* Let indices[0] point to index of k2 */
            size_t *bufIndices = indices[0];
            indices[0] = &offset[0][2];

            /* Index of nu23 */
            offset[1][2] = _kernels_work_index_nu(kernWork -> kernOrder, indices, size);

            /* Reset indices[0] and size[0] */
            indices[0] = bufIndices;
            size[0] = bufSize;


            /**  Calculate nu = (k1_ + k2_).k3_ / (|k1_ + k2_| k3) = (nu13 * k1 + nu23 * k2) / |k1_ + k2_| */

            /* Calculate nu */
            kernWork -> nu[offset[1][0]] = (kernWork -> nu[offset[1][1]] * kernWork -> k[offset[0][1]] + kernWork -> nu[offset[1][2]] * kernWork -> k[offset[0][2]]) / kernWork -> k[offset[0][0]];
          }

        /* m1 = m2 */
        else
          {
            /**  Indices of k2 and nu12  **/

            /* Reduce size[1] for k2_ = kj1_ + ... + kjm2-1_ */
            size[1] -= 1;

            /* Index of k2 */
            offset[0][1] = _kernels_work_index_k(kernWork -> kernOrder, indices[1], size[1]);

            /* Index of nu12 */
            offset[1][1] = _kernels_work_index_nu(kernWork -> kernOrder, indices, size);

            /* Reset size[1] */
            size[1] += 1;


            /**  Indices of k3 and nu13, and value of nu23  **/

            /* Reduce size[1] for k3_ = kjm2_ */
            size_t bufSize = size[1] - 1;
            size[1] = 1;

            /* Let indices[1] point to index of k3 */
            indices[1] += bufSize;

            /* Index of k3 */
            offset[0][2] = indices[1][0];

            /* Index of nu13 */
            offset[1][2] = _kernels_work_index_nu(kernWork -> kernOrder, indices, size);


            /* Set size[0] for k2_ = kj1_ + ... + kjm2-1_ */
            size_t bufSize_ = size[0];
            size[0] = bufSize;

            /* Let indices[0] point to index of kj1_ */
            size_t *bufIndices = indices[0];
            indices[0] = indices[1] - bufSize;

            /* Value of nu23 */
            double nu23 = kernWork -> nu[_kernels_work_index_nu(kernWork -> kernOrder, indices, size)];

            /* Reset indices and size */
            indices[0] = bufIndices;
            indices[1] -= bufSize;
            size[0] = bufSize_;
            size[1] = bufSize + 1;


            /**  Calculate nu = k1_.(k2_ + k3_) / (k1 |k2_ + k3_|) = (nu12 * k2 + nu13 * k3) / sqrt(k2^2 + k3^2 + 2 * k2 * k3 * nu23)  -- TODO: This way k23 is calculated twice... Must first finish calculation of all kj1...jm2 before starting the m1=m2 nu calculations (only used for kernOrder >= 4, i.e., not important for one-loop Power Spectrum).  **/

            /* Calculate k23 = |k2_ + k3_| = sqrt(k2^2 + k3^2 + 2 * k2 * k3 * nu23) */
            double k2 = kernWork -> k[offset[0][1]];
            double k3 = kernWork -> k[offset[0][2]];
            double k23 = sqrt(k2*k2 + k3*k3 + 2. * k2*k3*nu23);

            /* Calculate nu */
            kernWork -> nu[offset[1][0]] = (kernWork -> nu[offset[1][1]] * kernWork -> k[offset[0][1]] + kernWork -> nu[offset[1][2]] * kernWork -> k[offset[0][2]]) / k23;
          }

        /* Increment offset[1][0] */
        offset[1][0] += 1;
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


static int _kernels_work_var_k_cp(kern_t *kernWorkIn, kern_t *kernWorkOut, size_t **indices, size_t *size, size_t **offset, size_t *depth)
{
    /*

        Copy all possible

            k_ = ki1_ + ... + kim1_

        from kernWorkIn to kernWorkOut. Must have kernWork2 -> kernOrder >= kernWork1 -> kernOrder.

    */

    /* Loop through all possible i1, ..., im1 */
    if (*depth < *size)
      {
        /* Initial index for i ; Must start at 0 for i1, and at in-1 + 1 for in with n = *depth */
        size_t ini = (*depth == 0) ? 0 : indices[0][*depth - 1] + 1;

        for (size_t i = ini; i < kernWorkOut -> kernOrder; i++)
          {
            /* Index at location *depth in indices array is i */
            indices[0][*depth] = i;

            /* Get the next index at *depth + 1 */
            *depth += 1;
            _kernels_work_var_k_cp(kernWorkIn, kernWorkOut, indices, size, offset, depth);

            /* Reset *depth to adjust the index in indices array */
            *depth -= 1;
          }
      }

    /* Copy ki1_ + ... + kim_ */
    else
      {
        /**

            Copy k_ = k1_ + k2_, where k1_ = ki1_ + ... + kim1-1_ and k2_ = kim1_

        **/

        /* Copy k and mu */
        size_t index = _kernels_work_index_k(kernWorkIn -> kernOrder, *indices, *size);

        kernWorkOut -> k[offset[0][0]] = kernWorkIn -> k[index];
        kernWorkOut -> mu[offset[0][0]] = kernWorkIn -> mu[index];


        /**

            Copy all possible (ki1_ + ... + kim1_).(kj1_ + ... + kjm2_) / (|ki1_ + ... + kim1_| * |kj1_ + ... + kjm2_|)

        **/

        /* Set depth to 0 for nu functions */
        *depth = 0;

        /* Restrict m2 to m2 <= m1 */
        for (size[1] = 1; size[1] <= size[0]; size[1]++)
          {
            _kernels_work_var_nu_cp(kernWorkIn, kernWorkOut, indices, size, offset, depth);
          }

        /* Reset depth */
        *depth = *size;

        /* Incrmeent offset[0][0] */
        offset[0][0] += 1;
      }

    return 0;
}


static int _kernels_work_var_nu_cp(kern_t *kernWorkIn, kern_t *kernWorkOut, size_t **indices, size_t *size, size_t **offset, size_t *depth)
{
    /*

        Copy all possible

            nu = (ki1_ + ... + kim1_).(kj1_ + ... + kjm2_) / (|ki1_ + ... + kim1_| * |kj1_ + ... + kjm2_|)

        from kernWorkIn to kernWorkOut. Must have kernWork2 -> kernOrder >= kernWork1 -> kernOrder.

    */

    /* Loop through all possible j1, ..., jm1 */
    if (*depth < size[1])
      {
        /* Keep track of the i1, ..., im1 ; indices[1][*depth] holds the value of i from the previous recursion */
        size_t i = (*depth == 0) ? ( (size[0] == size[1]) ? 1 : 0) : indices[1][*depth];

        /* Initial index for j */
        size_t ini = (*depth == 0) ? ( (size[0] == size[1]) ? indices[0][0] + 1 : 0) : indices[1][*depth - 1] + 1;

        for (size_t j = ini; j < kernWorkOut -> kernOrder; j++)
          {
            /* Cannot have indices that are occupied by i1, ..., im1 */
            if (i < size[0] && j == indices[0][i])
              {
                i++;
                continue;
              }

            /* Index at location *depth in indices array is j */
            indices[1][*depth] = j;

            /* Store the current value of i in indices[1][*depth + 1] to initialise i in the next recursion */
            if (*depth < size[1])
                indices[1][*depth + 1] = i;

            /* Get the next index at *depth + 1 */
            *depth += 1;
            _kernels_work_var_nu_cp(kernWorkIn, kernWorkOut, indices, size, offset, depth);

            /* Reset *depth */
            *depth -= 1;
          }
      }

    else
      {
        /* Copy nu */
        kernWorkOut -> nu[offset[1][0]] = kernWorkIn -> nu[_kernels_work_index_nu(kernWorkIn -> kernOrder, indices, size)];

        /* Increment offset[1][0] */
        offset[1][0] += 1;
      }

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


static size_t _kernels_work_index_k(size_t kernOrder, size_t *indices, size_t size)
{
    /*

        Get the index in the kernWork -> k and kernWork -> mu arrays

        The array indices must contain the indices of the summed variables

            k_ = ki1_ + ... + kim_

        such that

            indices = {i1, ..., im}

        with

            i1 < ... < im <= kernOrder.

        Here, m = size.

    */

    /* Default case */
    if (size == 1)
        return indices[0];

    /* Index */
    size_t index = 0;

    /* Offset all highest blocks */
    for (size_t i = 1; i < size; i++)
      {
        index += misc_binom(kernOrder, i);
      }

    /* Offset to deepest block */
    for (size_t n = 0; n < size; n++)
      {
        size_t ini = (n == 0) ? 0 : indices[n - 1] + 1;

        for (size_t i = ini; i < indices[n]; i++)
          {
            index += misc_binom(kernOrder - i - 1, size - n - 1);
          }
      }

    return index;
}


static size_t _kernels_work_index_nu(size_t kernOrder, size_t **indices, size_t *size)
{
    /*

        Get the index in the kernWork -> nu array

        The array indices must contain the indices of the summed variables

            nu(ki1_ + ... + kim1_, kj1_ + ... + kjm2_)

        such that

            indices = {{i1, ..., im1}, {j1, ..., jm2}}

        with

            i1 < ... < im1 <= *kernOrder
            j1 < ... < jm2 <= *kernOrder

        and

            m1 <= m2.

        If m1 == m2, i1 < j1, and if m1 != m2, j1 may be less than i1, but does not have to be.

        Here, size = {m1, m2}.

    */

    /* Must make sure that i1 < j1 for m1 == m2 */
    bool swap = false;

    if (size[0] == size[1] && indices[0][0] > indices[1][0])
      {
        size_t *buf = indices[0];
        indices[0] = indices[1];
        indices[1] = buf;

        swap = true;
      }

    /* Index */
    size_t index = 0;

    /* Offset to correct m1 in {i1, ..., im1; j1, ..., jm2} block */
    if (size[0] != 1)
      {
        /* Indices to be passed to misc_multinom() */
        size_t multinomialIndices[3];

        for (size_t i1 = 1; i1 < size[0]; i1++)
          {
            /* First multinomial index is i1 */
            multinomialIndices[0] = i1;

            for (size_t i2 = 1; i2 <= i1; i2++)
              {
                /* Second multinomial index is at i2 */
                multinomialIndices[1] = i2;

                /* Final multinomial index is at n - i1 - i2 */
                multinomialIndices[2] = kernOrder - i1 - i2;

                /* Offset is number of ways to choose i1 and i2 out of n without choosing duplicates; If i1 == i2 we must therefore divide by 2 to correct for overcounting */
                index += misc_multinom(kernOrder, multinomialIndices, 3) / ((i1 == i2) ? 2 : 1);
              }
          }
      }

    /* Offset to correct i1, ..., im1 */
    for (size_t n = 0; n < size[0]; n++)
      {
        /* Initial index of i ; i1' must start at 0, while all higher i2', ..., im1' must start at i1 + , ..., im1-1 + 1 (i.e., the increment of the previous final index), respectively */
        size_t ini = (n == 0) ? 0 : indices[0][n - 1] + 1;

        for (size_t i = ini; i < indices[0][n]; i++)
          {
            /* Number of possible ways to draw in+1', ..., im1' for given in' */
            size_t binom1 = misc_binom(kernOrder - i - 1, size[0] - n - 1);
            size_t binom2 = 0;

            for (size_t m = 1; m <= size[0]; m++)
              {
                /* Adjust the number of values j1', ..., jm2'' can be drawn from ; n - m1 - i1 or n - m1 - i1' if m1 = m2', depending on whether n > 1 or n = 0, respectively, or n - m1 if m1 = m2' */
                size_t kernOrderAdjusted = kernOrder - size[0] - ( (size[0] != m) ? 0 : (n == 0) ? i : indices[0][0] );

                /* Not enough values to choose from */
                if (m > kernOrderAdjusted)
                    continue;

                /* Number of possible ways to draw j1', ..., jm2'' for given i1', ..., im1' */
                binom2 += misc_binom(kernOrderAdjusted, m);
              }

            index += binom1 * binom2;
          }
      }

    /* Offset to correct m2 in {i1, ..., im1; j1, ..., jm2} block */
    for (size_t m = 1; m < size[1]; m++)
      {
        size_t offset = (size[0] != m) ? 0 : indices[0][0];

        index += misc_binom(kernOrder - size[0] - offset, m);
      }

    /* Offset */

    /* Index for {i1, ..., im1} array (starts at 1 if m1 == m2, since j loop starts at i1 + 1 for m1 == m2) */
    size_t i = (size[0] == size[1]) ? 1 : 0;

    /* Index of {j1, ..., jm} array when correcting for by jn -> jn - j1 and removing the holes due to {i1, ..., im} array */
    size_t k = 0;

    /* Adjust the number of values j1', ..., jm2' can be drawn from ; n - m1 - i1 if m1 = m2, or n - m1 if m1 = m2 */
    size_t kernOrderAdjusted = kernOrder - size[0] - ( (size[0] != size[1]) ? 0 : indices[0][0] );

    for (size_t n = 0; n < size[1]; n++)
      {
        /* Initial index of j ; j1' must start at i1 + 1 is m1 == m2, and else at 0, while all higher j2', ..., jm' must start at j1 + , ..., jm1-1 + 1, respectively */
        size_t ini = (n == 0) ? ( (size[0] == size[1]) ? indices[0][0] + 1 : 0 ) : indices[1][n - 1] + 1;

        for (size_t j = ini; j < indices[1][n]; j++)
          {
            /* Cannot count indices already contained in {i1, ..., im1} */
            if (i < size[0] && j == indices[0][i])
              {
                i++;
                continue;
              }

            /* Number of possible ways to draw jn+1', ..., jm2' for given jn' */
            index += misc_binom(kernOrderAdjusted - ++k, size[1] - n - 1);
          }
      }

    /* Undo a possible swap */
    if (swap)
      {
        size_t *buf = indices[0];
        indices[0] = indices[1];
        indices[1] = buf;
      }

    return index;
}


#if false

static size_t _kernels_work_index_nu_(size_t kernOrder, size_t **indices, size_t *size)
{
    /*

        !!! OLD VERSION !!!

        Get the index in the kernWork -> nu array

        The array indices must contain the indices of the summed variables

            nu(ki1_ + ... + kim1_, kj1_ + ... + kjm2_)

        such that

            indices = {{i1, ..., im1}, {j1, ..., jm2}}

        with

            i1 < ... < im1 <= *kernOrder
            j1 < ... < jm2 <= *kernOrder

        and

            m1 <= m2.

        If m1 == m2, i1 < j1, and if m1 != m2, j1 may be less than i1, but does not have to be.

        Here, size = {m1, m2}.

    */

    /* Index */
    size_t index = 0;

    /* Offset all higher blocks */
    if (size[0] != 1)
      {
        /* Indices to be passed to misc_multinom() */
        size_t multinomialIndices[3];

        /* Final index of i1 ; If m2 == 1, the last block has i1 = m1-1, but if m2 != 1, the last block has i1 = m1 */
        size_t fin1 = (size[1] == 1) ? size[0] - 1 : size[0];

        for (size_t i1 = 1; i1 <= fin1; i1++)
          {
            /* First multinomial index is i1 */
            multinomialIndices[0] = i1;

            /* Final index of i2 ; If m2 == 1 (max is always i1, since i1 >= i2) or i1 < fin1 (i.e., not yet at correct block), the last block has i2 = i1, and else i2 = m2 - 1 (one before last block) */
            size_t fin2 = (i1 < fin1 || size[1] == 1) ? i1 : size[1] - 1;

            for (size_t i2 = 1; i2 <= fin2; i2++)
              {
                /* Second multinomial index is at i2 */
                multinomialIndices[1] = i2;

                /* Final multinomial index is at n - i1 - i2 */
                multinomialIndices[2] = kernOrder - i1 - i2;

                /* Offset is number of ways to choose i1 and i2 out of n without choosing duplicates; If i1 == i2 we must therefore divide by 2 to correct for overcounting */
                index += misc_multinom(kernOrder, multinomialIndices, 3) / ((i1 == i2) ? 2 : 1);
              }
          }
      }

    /* Offset to correct i1, ..., im1 */
    for (size_t n = 0; n < size[0]; n++)
      {
        /* Initial index of i ; i1' must start at 0, while all higher i2', ..., im1' must start at i1 + , ..., im1-1 + 1 (i.e., the increment of the previous final index), respectively */
        size_t ini = (n == 0) ? 0 : indices[0][n - 1] + 1;

        for (size_t i = ini; i < indices[0][n]; i++)
          {
            /* Number of indices to choose for j1, ..., jm2 must be reduced by i1' if n = 0, or else i1, to ensure only i1' < j1' are counted for m1 == m2 */
            size_t offset = (size[0] != size[1]) ? 0 : (n == 0) ? i : indices[0][0];

            /* For fixed in1, can only choose from n - i' - 1 (must have -1 since i1' starts at 0) values for in+1', ..., im1', and n - m1 - offset for j1', ..., jm2' */
            index += misc_binom(kernOrder - i - 1, size[0] - n - 1) * misc_binom(kernOrder - size[0] - offset, size[1]);
          }
      }

    /* Offset to correct j1, ..., jm2 */

    /* Index for {i1, ..., im1} array (starts at 1 if m1 == m2, since j loop starts at i1 + 1 for m1 == m2) */
    size_t i = (size[0] == size[1]) ? 1 : 0;

    /* Index of {j1, ..., jm} array when correcting for by jn -> jn - j1 and removing the holes due to {i1, ..., im} array */
    size_t k = 0;

    for (size_t n = 0; n < size[1]; n++)
      {
        /* Initial index of j ; j1' must start at i1 + 1 is m1 == m2, and else at 0, while all higher j2', ..., jm' must start at j1 + , ..., jm1-1 + 1, respectively */
        size_t ini = (n == 0) ? ( (size[0] == size[1]) ? indices[0][0] + 1 : 0 ) : indices[1][n - 1] + 1;

        for (size_t j = ini; j < indices[1][n]; j++)
          {
            /* Cannot count indices already contained in {i1, ..., im1} */
            if (j == indices[0][i])
              {
                i++;
                continue;
              }

            /* Number of indices to choose for j1, ..., jm2 must be reduced by i1 to ensure only i1 < j1' are counted for m1 == m2 */
            size_t offset = (size[0] != size[1]) ? 0 : indices[0][0];

            /* For fixed jn', can only choose from n - m1 - offset - k - 1 values for jn+1', ..., jm2' */
            index += misc_binom(kernOrder - size[0] - offset - k++ - 1, size[1] - n - 1);
          }
      }

    return index;
}

#endif



/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------    (External) Kern Struct    ------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


kern_t *kernels_new(size_t specOrder, size_t loopOrder)
{
    /*

        Create a new kern_t struct at a given kernel order (see master's thesis):

            kernOrder = specOrder + 2 * loopOrder - 1

        Have (n == max(specOrder, kernOrder))

                n = (n choose 1)

        values for k, mu and

                n (n - 1) / 2 = (n choose 2)

        values for nu.

    */

    kern_t *kern = malloc(sizeof(kern_t));

    /* Kern and Spec order */
    kern -> specOrder = specOrder;
    kern -> kernOrder = specOrder + 2 * loopOrder - 1;

    kern -> maxOrder = kern -> specOrder > kern -> kernOrder ? kern -> specOrder : kern -> kernOrder;

    /* Variables */
    kern -> z = 0.;

    kern -> k = calloc(kern -> maxOrder, sizeof(double));
    kern -> mu = calloc(kern -> maxOrder, sizeof(double));
    kern -> nu = calloc(kern -> maxOrder * (kern -> maxOrder - 1) / 2, sizeof(double));

    /* Fiducials */
    kern -> growth = 0.;

    kern -> lcdm = NULL;
    kern -> btst = NULL;
    kern -> bias = NULL;
    kern -> rsd = NULL;
    kern -> ctr = NULL;
    kern -> surv = NULL;

    /* Working kern_t structs */
    kern_t **kernWork = malloc(sizeof(kern_t*) * kern -> kernOrder);

    for (size_t i = 0; i < kern -> kernOrder; i++)
        kernWork[i] = _kernels_work_new(i + 1);

    kern -> kernWork = kernWork;

    /* Need to compute variables for working kern_t structs */
    kern -> computeWork = true;

    return kern;
}


kern_t *kernels_new_order(size_t specOrder, size_t kernOrder)
{
    /*

        Create a new kern_t struct at a given kernel order.

        Have (n == max(specOrder, kernOrder))

                n = (n choose 1)

        values for k, mu and

                n (n - 1) / 2 = (n choose 2)

        values for nu.

    */

    kern_t *kern = malloc(sizeof(kern_t));

    /* Kern and Spec order */
    kern -> specOrder = specOrder;
    kern -> kernOrder = kernOrder;

    kern -> maxOrder = kern -> specOrder > kern -> kernOrder ? kern -> specOrder : kern -> kernOrder;

    /* Variables */
    kern -> z = 0.;

    kern -> k = calloc(kern -> maxOrder, sizeof(double));
    kern -> mu = calloc(kern -> maxOrder, sizeof(double));
    kern -> nu = calloc(kern -> maxOrder * (kern -> maxOrder - 1) / 2, sizeof(double));

    /* Fiducials */
    kern -> growth = 0.;

    kern -> lcdm = NULL;
    kern -> btst = NULL;
    kern -> bias = NULL;
    kern -> rsd = NULL;
    kern -> ctr = NULL;
    kern -> surv = NULL;

    /* Working kern_t structs */
    kern_t **kernWork = malloc(sizeof(kern_t*) * kern -> kernOrder);

    for (size_t i = 0; i < kern -> kernOrder; i++)
        kernWork[i] = _kernels_work_new(i + 1);

    kern -> kernWork = kernWork;

    /* Need to compute variables for working kern_t structs */
    kern -> computeWork = true;

    return kern;
}


kern_t *kernels_free(kern_t *kern)
{
    /*

        Free kern

    */

    /* Check for NULL */
    if (kern == NULL)
        return NULL;

    /* Variables */
    free(kern -> k);
    free(kern -> nu);
    free(kern -> mu);

    /* Fiducials */
    kern -> lcdm = fid_lcdm_free(kern -> lcdm);
    kern -> btst = fid_btst_free(kern -> btst);
    kern -> bias = fid_bias_free(kern -> bias);
    kern -> rsd = fid_rsd_free(kern -> rsd);
    kern -> ctr = fid_ctr_free(kern -> ctr);
    kern -> surv = fid_surv_free(kern -> surv);

    /* Working kern_t structs */
    kern_t **kernWork = (kern_t**) kern -> kernWork;

    for (size_t i = 0; i < kern -> kernOrder; i++)
        kernWork[i] = _kernels_work_free(kernWork[i]);

    free(kernWork);

    /* Free kern itself */
    free(kern);

    /* Set kern to NULL so calling free again is no problem */
    kern = NULL;

    return kern;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int kernels_set_z(kern_t *kern, double z)
{
    /*

        Set the redshift z of kern

    */

    /* Set the redshift */
    kern -> z = z;

    /* Set the fiducials */

    /* Growth */
    kern -> growth = _fidGrowth_(&z, _fidParamsGrowth_);

    /* LCDM */
    kern -> lcdm = fid_lcdm_free(kern -> lcdm);
    kern -> lcdm = _fidLCDM_(&z, _fidParamsLCDM_);

    /* BTST */
    kern -> btst = fid_btst_free(kern -> btst);
    kern -> btst = _fidBTST_(&z, _fidParamsBTST_);

    /* Bias */
    kern -> bias = fid_bias_free(kern -> bias);
    kern -> bias = _fidBias_(&z, _fidParamsBias_);

    /* RSD */
    kern -> rsd = fid_rsd_free(kern -> rsd);
    kern -> rsd = _fidRSD_(&z, _fidParamsRSD_);

    /* Ctr */
    kern -> ctr = fid_ctr_free(kern -> ctr);
    kern -> ctr = _fidCtr_(&z, _fidParamsCtr_);

    /* Surv */
    kern -> surv = fid_surv_free(kern -> surv);
    kern -> surv = _fidSurv_(&z, _fidParamsSurv_);

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int kernels_set_k(kern_t *kern, size_t index, double k)
{
    /*

        Set the index'th mode k of kern

    */

    if (index >= kern -> maxOrder)
      {
        printf("Cannot set the %ld'th k from a 'kern_t' struct if only %ld are expected.\n", index, kern -> maxOrder);
        exit(1);

        return 1;
      }

    /* Set the mode */
    kern -> k[index] = k;

    return 0;
}


int kernels_qset_k(kern_t *kern, size_t index, double k)
{
    /*

        Set the index'th mode k of kern (does not test for valid index)

    */

    /* Set the mode */
    kern -> k[index] = k;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int kernels_set_mu(kern_t *kern, size_t index, double mu)
{
    /*

        Set the index'th angle mu for kern

    */

    if (index >= kern -> maxOrder)
      {
        printf("Cannot set the %ld'th mu from a 'kern_t' struct if only %ld are expected.\n", index, kern -> maxOrder);
        exit(1);

        return 1;
      }

    /* Set the angle */
    kern -> mu[index] = mu;

    return 0;
}


int kernels_qset_mu(kern_t *kern, size_t index, double mu)
{
    /*

        Set the index'th angle mu for kern (does not test for valid index)

    */

    /* Set the angle */
    kern -> mu[index] = mu;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int kernels_set_nu(kern_t *kern, size_t index1, size_t index2, double nu)
{
    /*

        Set the angle nu(ki_, kj_) between different k modes for kern.

    */

    if (index1 == index2)
      {
        printf("Cannot set the angle between the same k vectors from a 'kern_t' struct (angle is always 1. and thus is ignored).\n");
        exit(1);

        return 1;
      }

    if (index1 >= kern -> maxOrder || index2 >= kern -> maxOrder)
      {
        printf("Cannot set the angle between the %ld'th and %ld'th k vectors from a 'kern_t' struct if only %ld k vectors are expected.\n", index1, index2, kern -> maxOrder);
        exit(1);

        return 1;
      }

    /* Set nu */
    kern -> nu[kernels_get_nu_index(kern -> maxOrder, index1, index2)] = nu;

    return 0;
}


int kernels_qset_nu(kern_t *kern, size_t index1, size_t index2, double nu)
{
    /*

        Set the angle nu(ki_, kj_) between different k modes for kern. (does not test for valid index1 and index2)

    */

    /* Set nu */
    kern -> nu[kernels_get_nu_index(kern -> maxOrder, index1, index2)] = nu;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_get_z(kern_t *kern)
{
    /*

        Get the redshift z from kern

    */

    return kern -> z;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_get_k(kern_t *kern, size_t index)
{
    /*

        Get the index'th k from kern

    */

    if (index >= kern -> maxOrder)
      {
        printf("Cannot get the %ld'th k from a 'kern_t' struct if only %ld are expected.\n", index, kern -> maxOrder);
        exit(1);

        return NAN;
      }

    return kern -> k[index];
}


double kernels_qget_k(kern_t *kern, size_t index)
{
    /*

        Get the index'th k from kern (does not test for valid index)

    */

    return kern -> k[index];
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_get_mu(kern_t *kern, size_t index)
{
    /*

        Get the index'th mu from kern

    */

    if (index >= kern -> maxOrder)
      {
        printf("Cannot get the %ld'th mu from a 'kern_t' struct if only %ld are expected.\n", index, kern -> maxOrder);
        exit(1);

        return NAN;
      }

    return kern -> mu[index];
}


double kernels_qget_mu(kern_t *kern, size_t index)
{
    /*

        Get the index'th mu from kern (does not test for valid index)

    */

    return kern -> mu[index];
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_get_nu(kern_t *kern, size_t index1, size_t index2)
{
    /*

        Get the nu angle between the index1'th and index2'th k vectors from kern

    */

    if (index1 == index2)
      {
        printf("Cannot get the angle between the same k vectors from a 'kern_t' struct (angle is always 1. and thus is ignored).\n");
        exit(1);

        return NAN;
      }

    if (index1 >= kern -> maxOrder || index2 >= kern -> maxOrder)
      {
        printf("Cannot get the angle between the %ld'th and %ld'th k vectors from a 'kern_t' struct only %ld k vectors are expected.\n", index1, index2, kern -> maxOrder);
        exit(1);

        return NAN;
      }

    return kern -> nu[kernels_get_nu_index(kern -> maxOrder, index1, index2)];
}


double kernels_qget_nu(kern_t *kern, size_t index1, size_t index2)
{
    /*

        Get the nu angle between the index1'th and index2'th k vectors from kern (does not test for valid index1 and index2)

    */

    return kern -> nu[kernels_get_nu_index(kern -> maxOrder, index1, index2)];
}


/*  ------------------------------------------------------------------------------------------------------  */


size_t kernels_get_nu_index(size_t order, size_t index1, size_t index2)
{
    /*

        Get the index of the nuij variable in the nu array

    */

    /* Let i be the smaller index and j the larger */
    size_t i = index1;
    size_t j = index2;

    if (i > j)
      {
        i = index2;
        j = index1;
      }

    return i * (order - 1) - i * (i + 1) / 2 + j - 1;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     BOOTSTRAP KERNELS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Kernel Components   ----------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_alpha(kern_t *kernVar)
{
    /*

        Antisymmetric alpha kernel:

            α^(a)(k1_, k2_) = α^(a)(k1, k2, nu12)
                            = 1/2 nu12 (k2/k1 - k1/k2)

    */

    double alpha = kernVar -> nu[0] / 2. * (kernVar -> k[1] / kernVar -> k[0] - kernVar -> k[0] / kernVar -> k[1]);

    return alpha;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dalpha_k(kern_t *kernVar)
{
    /*

        Derivative of antisymmetric alpha kernel with respect to k1 (first variable):

            dα^(a)/dk1(k1, k2, nu12) = -1/2 nu12 (k2/k1 + k1/k2) / k1

        To get the derivative w.r.t. k2 (second variable) simlpy exchange k1 <-> k2 in the k array and multiply the result by -1:

            dα^(a)/dk2(k1, k2, nu12) = -dα^(a)/dk2(k2, k1, nu12)

    */

    double dalpha = - kernVar -> nu[0] / 2. * (kernVar -> k[1] / kernVar -> k[0] + kernVar -> k[0] / kernVar -> k[1]) / kernVar -> k[0];

    return dalpha;
}


double kernels_dalpha_nu(kern_t *kernVar)
{
    /*

        Derivative of antisymmetric alpha kernel with respect to nu12:

            dα^(a)/dnu12(k1, k2, nu12) = 1/2 (k2/k1 - k1/k2)

    */

    double dalpha = 0.5 * (kernVar -> k[1] / kernVar -> k[0] - kernVar -> k[0] / kernVar -> k[1]);

    return dalpha;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_beta(kern_t *kernVar)
{
    /*

        Beta kernel:

            β(k1_, k2_) = β(k1, k2, nu12)
                        = (k1_ + k2_)^2 k1_.k2_ / (2 k1^2 k2^2) = nu12 (nu12 + 1/2 (k1/k2 + k2/k1))

    */

    double beta = kernVar -> nu[0] * ( kernVar -> nu[0] + (kernVar -> k[0] / kernVar -> k[1] + kernVar -> k[1] / kernVar -> k[0]) / 2. );

    return beta;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dbeta_k(kern_t *kernVar)
{
    /*

        Derivative of the beta kernel with respect to k1 (first variable):

            dβ/dk1(k1, k2, nu12) = nu12/2 (k1/k2 - k2/k1) / k1

        To get the derivative w.r.t. k2 (second variable) simlpy exchange k1 <-> k2 in the k array:

            dβ/dk2(k1, k2, nu12) = dβ/dk2(k2, k1, nu12)

    */

    double dbeta = kernVar -> nu[0] / 2. * (kernVar -> k[0] / kernVar -> k[1] - kernVar -> k[1] / kernVar -> k[0]) / kernVar -> k[0];

    return dbeta;
}


double kernels_dbeta_nu(kern_t *kernVar)
{
    /*

        Derivative of the beta function with respect to nu12:

            dβ/dnu12(k1, k2, nu12) = 2 nu12 + 1/2 (k1/k2 + k2/k1)

    */

    double dbeta = 2. * kernVar -> nu[0] + (kernVar -> k[0] / kernVar -> k[1] + kernVar -> k[1] / kernVar -> k[0]) / 2.;

    return dbeta;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_gamma(kern_t *kernVar)
{
    /*

        Gamma kernel:

            γ(k1_, k2_) = γ(nu12)
                        = 1 - nu12^2

    */

    double gamma = 1.0 - pow(kernVar -> nu[0], 2);

    return gamma;
}


double kernels_dgamma_nu(kern_t *kernVar)
{
    /*

        Derivative of the gamma kernel with respect to nu12:

            dγ/dnu12(nu12) = - 2 nu12

    */

    double dgamma = - 2. * kernVar -> nu[0];

    return dgamma;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Second order Kernels   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_h2(kern_t *kernVar, double *h2Kernels)
{
    /*

        Bootstrap second order density + velocity kernels:

            F2(k1_, k2_) = F2(k1, k2, nu12)
                         = β(k1, k2, nu12) + 1/2 a^(2)_γ γ(nu12)

            G2(k1_, k2_) = G2(k1, k2, nu12)
                         = β(k1, k2, nu12) + 1/2 d^(2)_γ γ(nu12)

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* β(k1_, k2_) */
    double beta = kernels_beta(kernVarH2);

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarH2);

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    h2Kernels[0] = beta + 0.5 * kernVar -> btst -> a2Ga * gamma;
    h2Kernels[1] = beta + 0.5 * kernVar -> btst -> d2Ga * gamma;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_dh2_k(kern_t *kernVar, double *dh2Kernels)
{
    /*

        Bootstrap second order density + velocity kernels derivative with respect to k1 (first variable):

            dF2/dk1(k1, k2, nu12) = dβ/dk1(k1, k2, nu12)

            dG2/dk1(k1, k2, nu12) = dβ/dk1(k1, k2, nu12)

        To get the derivative w.r.t. k2 (second variable) simlpy exchange k1 <-> k2 in the k array:

            dF2/dk2(k1, k2, nu12) = dF2/dk2(k2, k1, nu12)

            dG2/dk2(k1, k2, nu12) = dG2/dk2(k2, k1, nu12)

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* dβ/dk1(k1_, k2_) */
    double dbeta = kernels_dbeta_k(kernVarH2);

    /* dF2/dk1(k1_, k2_) + dG2/dk1(k1_, k2_) */
    dh2Kernels[0] = dbeta;
    dh2Kernels[1] = dbeta;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh2_nu(kern_t *kernVar, double *dh2Kernels)
{
    /*

        Bootstrap second order density + velocity kernels derivative with respect to nu12:

            dF2/dnu12(k1, k2, nu12) = dβ/dnu12(k1, k2, nu12) + 1/2 a^(2)_γ dγ/dnu12(nu12)

            dG2/dnu12(k1, k2, nu12) = dβ/dnu12(k1, k2, nu12) + 1/2 d^(2)_γ dγ/dnu12(nu12)

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* dβ/dnu12(k1_, k2_) */
    double dbeta = kernels_dbeta_nu(kernVarH2);

    /* dγ/dnu12(k1_, k2_) */
    double dgamma = kernels_dgamma_nu(kernVarH2);

    /* dF2/dnu12(k1_, k2_) + dG2/dnu12(k1_, k2_) */
    dh2Kernels[0] = dbeta + 0.5 * kernVar -> btst -> a2Ga * dgamma;
    dh2Kernels[1] = dbeta + 0.5 * kernVar -> btst -> d2Ga * dgamma;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_dh2_a2ga(kern_t *kernVar, double *dh2Kernels)
{
    /*

        Bootstrap second order density + velocity kernels derivative with respect to a^(2)_γ:

            dF2/d(a^(2)_γ)(k1, k2, nu12) = dF2/d(a^(2)_γ)(nu12) = 1/2 γ(nu12)

            dG2/d(a^(2)_γ)(k1, k2, nu12) = 0

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarH2);

    /* dF2/da^(2)_γ(k1_, k2_) + dG2/da^(2)_γ(k1_, k2_) */
    dh2Kernels[0] = 0.5 * gamma;
    dh2Kernels[1] = 0.;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh2_d2ga(kern_t *kernVar, double *dh2Kernels)
{
    /*

        Bootstrap second order density + velocity kernels derivative with respect to d^(2)_γ:

            dF2/d(d^(2)_γ)(k1, k2, nu12) = 0

            dG2/d(d^(2)_γ)(k1, k2, nu12) = dG2/d(d^(2)_γ)(nu12) = 1/2 γ(nu12)

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarH2);

    /* dF2/dd^(2)_γ(k1_, k2_) + dG2/dd^(2)_γ(k1_, k2_) */
    dh2Kernels[0] = 0.;
    dh2Kernels[1] = 0.5 * gamma;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Third order Kernels   ---------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_h3(kern_t *kernVar, double *h3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels:

            F3(k1_, k2_, k3_) = F3(k1, k2, k3, nu12, nu23, nu31)
                              = 1/3 ( β(k1_ + k2_, k3_) β(k1_, k2_) + β(k1_ + k3_, k2_) β(k1_, k3_) + β(k2_ + k3_, k1_) β(k2_, k3_) )
                                 + 1/6 (1/2 a^(3)_γa + a^(3)_γb) ( γ(k1_ + k2_, k3_) γ(k1_, k2_) + γ(k1_ + k3_, k2_) γ(k1_, k3_) + γ(k2_ + k3_, k1_) γ(k2_, k3_))
                                 + 1/6 (1/2 a^(3)_γa - a^(3)_γb) ( α^(a)(k1_ + k2_, k3_) γ(k1_, k2_) + α^(a)(k1_ + k3_, k2_) γ(k1_, k3_) + α^(a)(k2_ + k3_, k1_) γ(k2_, k3_) )
                                 - 1/6 (1/2 a^(3)_γa - a^(3)_γb - 2h) ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                 + 1/3 (1/2 a^(3)_γa - a^(3)_γb + a^(2)_γ - h) ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

            G3(k1_, k2_, k3_) = G3(k1, k2, k3, nu12, nu23, nu31)
                              = 1/3 ( β(k1_ + k2_, k3_) β(k1_, k2_) + β(k1_ + k3_, k2_) β(k1_, k3_) + β(k2_ + k3_, k1_) β(k2_, k3_) )
                                 + 1/6 (1/2 d^(3)_γa + d^(3)_γb) ( γ(k1_ + k2_, k3_) γ(k1_, k2_) + γ(k1_ + k3_, k2_) γ(k1_, k3_) + γ(k2_ + k3_, k1_) γ(k2_, k3_))
                                 + 1/6 (1/2 d^(3)_γa - d^(3)_γb) ( α^(a)(k1_ + k2_, k3_) γ(k1_, k2_) + α^(a)(k1_ + k3_, k2_) γ(k1_, k3_) + α^(a)(k2_ + k3_, k1_) γ(k2_, k3_) )
                                 - 1/6 (1/2 d^(3)_γa - d^(3)_γb - 2h) ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                 + 1/3 (1/2 d^(3)_γa - d^(3)_γb + d^(2)_γ - h) ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3, k12, k13, k23}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Bootstrap parameters */
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    h3Kernels[0] = 0.;
    h3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSin;

    double betaSum;
    double gammaSum;
    double alphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            alphaSum = kernels_alpha(kernVarH2);
            gammaSum = kernels_gamma(kernVarH2);

            /* F3(k1_, k2_, k3_) */
            h3Kernels[0] += 1./3. * (betaSum * betaSin)
                            + (btst -> a3GaA / 2. + btst -> a3GaB) / 6. * (gammaSum * gammaSin)
                            + (btst -> a3GaA / 2. - btst -> a3GaB) / 6. * (alphaSum * gammaSin)
                            - (btst -> a3GaA / 2. - btst -> a3GaB - 2. * btst -> h) / 6. * (betaSum * gammaSin)
                            + (btst -> a3GaA / 2. - btst -> a3GaB + btst -> a2Ga - btst -> h) / 3. * (gammaSum * betaSin);

            /* G3(k1_, k2_, k3_) */
            h3Kernels[1] += 1./3. * (betaSum * betaSin)
                            + (btst -> d3GaA / 2. + btst -> d3GaB) / 6. * (gammaSum * gammaSin)
                            + (btst -> d3GaA / 2. - btst -> d3GaB) / 6. * (alphaSum * gammaSin)
                            - (btst -> d3GaA / 2. - btst -> d3GaB - 2. * btst -> h) / 6. * (betaSum * gammaSin)
                            + (btst -> d3GaA / 2. - btst -> d3GaB + btst -> d2Ga - btst -> h) / 3. * (gammaSum * betaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_dh3_k(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernel derivative with respect to k1 (first variable)s:

            dF3/dk1(k1_, k2_, k3_) = dF3/dk1(k1, k2, k3, nu12, nu23, nu31)
                                   = 1/3 ( d(β(k1_ + k2_, k3_) β(k1_, k2_))/dk1 + d(β(k1_ + k3_, k2_) β(k1_, k3_))/dk1 + dβ/dk1(k2_ + k3_, k1_) β(k2_, k3_) )
                                      + 1/6 (1/2 a^(3)_γa + a^(3)_γb) ( dγ/dk1(k1_ + k2_, k3_) γ(k1_, k2_) + dγ/dk1(k1_ + k3_, k2_) γ(k1_, k3_))
                                      + 1/6 (1/2 a^(3)_γa - a^(3)_γb) ( dα^(a)/dk1(k1_ + k2_, k3_) γ(k1_, k2_) + dα^(a)/dk1(k1_ + k3_, k2_) γ(k1_, k3_) + dα^(a)/dk1(k2_ + k3_, k1_) γ(k2_, k3_) )
                                      - 1/6 (1/2 a^(3)_γa - a^(3)_γb - 2h) ( dβ/dk1(k1_ + k2_, k3_) γ(k1_, k2_) + dβ/dk1(k1_ + k3_, k2_) γ(k1_, k3_) + dβ/dk1(k2_ + k3_, k1_) γ(k2_, k3_) )
                                      + 1/3 (1/2 a^(3)_γa - a^(3)_γb + a^(2)_γ - h) ( d(γ(k1_ + k2_, k3_) β(k1_, k2_))/dk1 + d(γ(k1_ + k3_, k2_) β(k1_, k3_))/dk1 )

            dG3/dk1(k1_, k2_, k3_) = dG3/dk1(k1, k2, k3, nu12, nu23, nu31)
                                   = 1/3 ( d(β(k1_ + k2_, k3_) β(k1_, k2_))/dk1 + d(β(k1_ + k3_, k2_) β(k1_, k3_))/dk1 + dβ/dk1(k2_ + k3_, k1_) β(k2_, k3_) )
                                      + 1/6 (1/2 d^(3)_γa + d^(3)_γb) ( dγ/dk1(k1_ + k2_, k3_) γ(k1_, k2_) + dγ/dk1(k1_ + k3_, k2_) γ(k1_, k3_))
                                      + 1/6 (1/2 d^(3)_γa - d^(3)_γb) ( dα^(a)/dk1(k1_ + k2_, k3_) γ(k1_, k2_) + dα^(a)/dk1(k1_ + k3_, k2_) γ(k1_, k3_) + dα^(a)/dk1(k2_ + k3_, k1_) γ(k2_, k3_) )
                                      - 1/6 (1/2 d^(3)_γa - d^(3)_γb - 2h) ( dβ/dk1(k1_ + k2_, k3_) γ(k1_, k2_) + dβ/dk1(k1_ + k3_, k2_) γ(k1_, k3_) + dβ/dk1(k2_ + k3_, k1_) γ(k2_, k3_) )
                                      + 1/3 (1/2 d^(3)_γa - d^(3)_γb + d^(2)_γ - h) ( d(γ(k1_ + k2_, k3_) β(k1_, k2_))/dk1 + d(γ(k1_ + k3_, k2_) β(k1_, k3_))/dk1 )

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Bootstrap parameters */
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double dbetaSin;
    double gammaSin;

    double betaSum;
    double dbetaSum;
    double gammaSum;
    double dgammaSum;
    double dalphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            dbetaSin = (i == 0) ? kernels_dbeta_k(kernVarH2) : 0.;

            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* d(k12)/dk1 and d(k13)/dk1 */
            double dk = (i == 0) ? ( kernVarH2 -> k[0] + kernVarH2 -> k[1] * kernVarH2 -> nu[0] ) / kernVarH3 -> k[3 + index] : 0.;

            /* Variables ; Must swap locations k23 <-> k1 for derivatives (only acting on first) */
            kernVarH2 -> k[0] = (i == 0) ? kernVarH3 -> k[3 + index] : kernVarH3 -> k[2 - index];
            kernVarH2 -> k[1] = (i == 0) ? kernVarH3 -> k[2 - index] : kernVarH3 -> k[3 + index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* d(nu12,3)/dk1 and d(nu13,2)/dk1 */
            double dnu = (i == 0) ? ( kernVarH3 -> nu[1 - index] - dk * kernVarH2 -> nu[0] ) / kernVarH3 -> k[3 + index] : 0.;

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            dbetaSum = (i == 0) ? kernels_dbeta_k(kernVarH2) * dk + kernels_dbeta_nu(kernVarH2) * dnu : kernels_dbeta_k(kernVarH2);

            dalphaSum = (i == 0) ? kernels_dalpha_k(kernVarH2) * dk + kernels_dalpha_nu(kernVarH2) * dnu : -kernels_dalpha_k(kernVarH2); // Have additional - sign due to antisymmetry for k23 <-> k1

            gammaSum = kernels_gamma(kernVarH2);
            dgammaSum = (i == 0) ? kernels_dgamma_nu(kernVarH2) * dnu : 0.;

            /* dF3/dk1(k1_, k2_, k3_) */
            dh3Kernels[0] += 1./3. * (dbetaSum * betaSin + betaSum * dbetaSin)
                            + (btst -> a3GaA / 2. + btst -> a3GaB) / 6. * (dgammaSum * gammaSin)
                            + (btst -> a3GaA / 2. - btst -> a3GaB) / 6. * (dalphaSum * gammaSin)
                            - (btst -> a3GaA / 2. - btst -> a3GaB - 2. * btst -> h) / 6. * (dbetaSum * gammaSin)
                            + (btst -> a3GaA / 2. - btst -> a3GaB + btst -> a2Ga - btst -> h) / 3. * (dgammaSum * betaSin + gammaSum * dbetaSin);

            /* dG3/dk1(k1_, k2_, k3_) */
            dh3Kernels[1] += 1./3. * (dbetaSum * betaSin + betaSum * dbetaSin)
                            + (btst -> d3GaA / 2. + btst -> d3GaB) / 6. * (dgammaSum * gammaSin)
                            + (btst -> d3GaA / 2. - btst -> d3GaB) / 6. * (dalphaSum * gammaSin)
                            - (btst -> d3GaA / 2. - btst -> d3GaB - 2. * btst -> h) / 6. * (dbetaSum * gammaSin)
                            + (btst -> d3GaA / 2. - btst -> d3GaB + btst -> d2Ga - btst -> h) / 3. * (dgammaSum * betaSin + gammaSum * dbetaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh3_nu(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to nu12:

            dF3/dnu12(k1_, k2_, k3_) = dF3/dnu12(k1, k2, k3, nu12, nu23, nu31)
                                     = 1/3 ( d(β(k1_ + k2_, k3_) β(k1_, k2_))/dnu12 + dβ/dnu12(k1_ + k3_, k2_) β(k1_, k3_) + dβ/dnu12(k2_ + k3_, k1_) β(k2_, k3_) )
                                        + 1/6 (1/2 a^(3)_γa + a^(3)_γb) ( d(γ(k1_ + k2_, k3_) γ(k1_, k2_))/dnu12 + dγ/dnu12(k1_ + k3_, k2_) γ(k1_, k3_) + dγ/dnu12(k2_ + k3_, k1_) γ(k2_, k3_))
                                        + 1/6 (1/2 a^(3)_γa - a^(3)_γb) ( d(α^(a)(k1_ + k2_, k3_) γ(k1_, k2_))/dnu12 + dα^(a)/dnu12(k1_ + k3_, k2_) γ(k1_, k3_) + dα^(a)/dnu12(k2_ + k3_, k1_) γ(k2_, k3_) )
                                        - 1/6 (1/2 a^(3)_γa - a^(3)_γb - 2h) ( d(β(k1_ + k2_, k3_) γ(k1_, k2_))/dnu12 + dβ/dnu12(k1_ + k3_, k2_) γ(k1_, k3_) + dβ/dnu12(k2_ + k3_, k1_) γ(k2_, k3_) )
                                        + 1/3 (1/2 a^(3)_γa - a^(3)_γb + a^(2)_γ - h) ( d(γ(k1_ + k2_, k3_) β(k1_, k2_))/dnu12 + dγ/dnu12(k1_ + k3_, k2_) β(k1_, k3_) + dγ/dun12(k2_ + k3_, k1_) β(k2_, k3_) )

            dG3/dnu12(k1_, k2_, k3_) = dG3/dnu12(k1, k2, k3, nu12, nu23, nu31)
                                     = 1/3 ( d(β(k1_ + k2_, k3_) β(k1_, k2_))/dnu12 + dβ/dnu12(k1_ + k3_, k2_) β(k1_, k3_) + dβ/dnu12(k2_ + k3_, k1_) β(k2_, k3_) )
                                        + 1/6 (1/2 d^(3)_γa + d^(3)_γb) ( d(γ(k1_ + k2_, k3_) γ(k1_, k2_))/dnu12 + dγ/dnu12(k1_ + k3_, k2_) γ(k1_, k3_) + dγ/dnu12(k2_ + k3_, k1_) γ(k2_, k3_))
                                        + 1/6 (1/2 d^(3)_γa - d^(3)_γb) ( d(α^(a)(k1_ + k2_, k3_) γ(k1_, k2_))/dnu12 + dα^(a)/dnu12(k1_ + k3_, k2_) γ(k1_, k3_) + dα^(a)/dnu12(k2_ + k3_, k1_) γ(k2_, k3_) )
                                        - 1/6 (1/2 d^(3)_γa - d^(3)_γb - 2h) ( d(β(k1_ + k2_, k3_) γ(k1_, k2_))/dnu12 + dβ/dnu12(k1_ + k3_, k2_) γ(k1_, k3_) + dβ/dnu12(k2_ + k3_, k1_) γ(k2_, k3_) )
                                        + 1/3 (1/2 d^(3)_γa - d^(3)_γb + d^(2)_γ - h) ( d(γ(k1_ + k2_, k3_) β(k1_, k2_))/dnu12 + dγ/dnu12(k1_ + k3_, k2_) β(k1_, k3_) + dγ/dun12(k2_ + k3_, k1_) β(k2_, k3_) )


        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Bootstrap parameters */
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double dbetaSin;
    double gammaSin;
    double dgammaSin;

    double betaSum;
    double dbetaSum;
    double gammaSum;
    double dgammaSum;
    double alphaSum;
    double dalphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            dbetaSin = (j == 0) ? kernels_dbeta_nu(kernVarH2) : 0.;

            gammaSin = kernels_gamma(kernVarH2);
            dgammaSin = (j == 0) ? kernels_dgamma_nu(kernVarH2) : 0.;


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* d(k12)/dnu12 */
            double dk = (j == 0) ? ( kernVarH2 -> k[0] * kernVarH2 -> k[1] ) / kernVarH3 -> k[3 + index] : 0.;

            /* Variables ; Must swap locations k23 <-> k1 for derivatives (only acting on first) */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* d(nu12,3)/dnu12, d(nu13,2)/dnu12 and d(nu23,1)/dnu12 */
            double dnu = (j == 0) ? -kernVarH2 -> nu[0] * dk / kernVarH2 -> k[0] : kernVarH3 -> k[index - 1] / kernVarH2 -> k[0];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            dbetaSum = (j == 0) ? kernels_dbeta_k(kernVarH2) * dk + kernels_dbeta_nu(kernVarH2) * dnu : kernels_dbeta_nu(kernVarH2) * dnu;

            alphaSum = kernels_alpha(kernVarH2);
            dalphaSum = (j == 0) ? kernels_dalpha_k(kernVarH2) * dk + kernels_dalpha_nu(kernVarH2) * dnu : kernels_dalpha_nu(kernVarH2) * dnu;

            gammaSum = kernels_gamma(kernVarH2);
            dgammaSum = kernels_dgamma_nu(kernVarH2) * dnu;

            /* dF3/dnu12(k1_, k2_, k3_) */
            dh3Kernels[0] += 1./3. * (dbetaSum * betaSin + betaSum * dbetaSin)
                            + (btst -> a3GaA / 2. + btst -> a3GaB) / 6. * (dgammaSum * gammaSin + gammaSum * dgammaSin)
                            + (btst -> a3GaA / 2. - btst -> a3GaB) / 6. * (dalphaSum * gammaSin + alphaSum * dgammaSin)
                            - (btst -> a3GaA / 2. - btst -> a3GaB - 2. * btst -> h) / 6. * (dbetaSum * gammaSin + betaSum * dgammaSin)
                            + (btst -> a3GaA / 2. - btst -> a3GaB + btst -> a2Ga - btst -> h) / 3. * (dgammaSum * betaSin + gammaSum * dbetaSin);

            /* dG3/dnu12(k1_, k2_, k3_) */
            dh3Kernels[1] += 1./3. * (dbetaSum * betaSin + betaSum * dbetaSin)
                            + (btst -> d3GaA / 2. + btst -> d3GaB) / 6. * (dgammaSum * gammaSin + gammaSum * dgammaSin)
                            + (btst -> d3GaA / 2. - btst -> d3GaB) / 6. * (dalphaSum * gammaSin + alphaSum * dgammaSin)
                            - (btst -> d3GaA / 2. - btst -> d3GaB - 2. * btst -> h) / 6. * (dbetaSum * gammaSin + betaSum * dgammaSin)
                            + (btst -> d3GaA / 2. - btst -> d3GaB + btst -> d2Ga - btst -> h) / 3. * (dgammaSum * betaSin + gammaSum * dbetaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_dh3_a2ga(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to a^(2)_γ:

            dF3/da^(2)_γ(k1_, k2_, k3_) = dF3/da^(2)_γ(k1, k2, k3, nu12, nu23, nu31)
                                        = 1/3 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

            dG3/da^(2)_γ(k1_, k2_, k3_) = dG3/da^(2)_γ(k1, k2, k3, nu12, nu23, nu31)
                                        = 0

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Bootstrap parameters (not needed) */
//    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/da^(2)_γ(k1_, k2_, k3_) */
            dh3Kernels[0] += 1. / 3. * (gammaSum * betaSin);

            /* dG3/da^(2)_γ(k1_, k2_, k3_) */
            dh3Kernels[1] += 0.;

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh3_d2ga(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to d^(2)_γ:

            dF3/dd^(2)_γ(k1_, k2_, k3_) = dF3/dd^(2)_γ(k1, k2, k3, nu12, nu23, nu31)
                                        = 0

            dG3/dd^(2)_γ(k1_, k2_, k3_) = dG3/dd^(2)_γ(k1, k2, k3, nu12, nu23, nu31)
                                        = 1/3 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Bootstrap parameters (not needed) */
//    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/da^(2)_γ(k1_, k2_, k3_) */
            dh3Kernels[0] += 0.;

            /* dG3/da^(2)_γ(k1_, k2_, k3_) */
            dh3Kernels[1] += 1. / 3. * (gammaSum * betaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh3_h(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to h:

            dF3/dh(k1_, k2_, k3_) = dF3/dh(k1, k2, k3, nu12, nu23, nu31)
                                  = 1/3 ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                    - 1/3 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

            dG3/dh(k1_, k2_, k3_) = dG3/dh(k1, k2, k3, nu12, nu23, nu31)
                                  = 1/3 ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                    - 1/3 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Bootstrap parameters (not needed) */
//    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSin;

    double betaSum;
    double gammaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/dh(k1_, k2_, k3_) */
            dh3Kernels[0] += 1. / 3. * (betaSum * gammaSin)
                            - 1. / 3. * (gammaSum * betaSin);

            /* dG3/dh(k1_, k2_, k3_) */
            dh3Kernels[1] += 1. / 3. * (betaSum * gammaSin)
                            - 1. / 3. * (gammaSum * betaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


/*  ------------------------------------------------------------------------------------------------------  */


int kernels_btst_dh3_a3gaa(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to a^(3)_γa:

            dF3/da^(3)_γa(k1_, k2_, k3_) = dF3/da^(3)_γa(k1, k2, k3, nu12, nu23, nu31)
                                         = 1/12 ( γ(k1_ + k2_, k3_) γ(k1_, k2_) + γ(k1_ + k3_, k2_) γ(k1_, k3_) + γ(k2_ + k3_, k1_) γ(k2_, k3_))
                                          + 1/12 ( α^(a)(k1_ + k2_, k3_) γ(k1_, k2_) + α^(a)(k1_ + k3_, k2_) γ(k1_, k3_) + α^(a)(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          - 1/12 ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          + 1/6 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

            dG3/da^(3)_γa(k1_, k2_, k3_) = dG3/da^(3)_γa(k1, k2, k3, nu12, nu23, nu31)
                                         = 0

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSin;

    double betaSum;
    double gammaSum;
    double alphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            alphaSum = kernels_alpha(kernVarH2);
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/da^(3)_γa(k1_, k2_, k3_) */
            dh3Kernels[0] += 1. / 12. * (gammaSum * gammaSin)
                            + 1. / 12. * (alphaSum * gammaSin)
                            - 1. / 12. * (betaSum * gammaSin)
                            + 1. / 6. * (gammaSum * betaSin);

            /* dG3/da^(3)_γa(k1_, k2_, k3_) */
            dh3Kernels[1] += 0.;

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh3_a3gab(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to a^(3)_γb:

            dF3/da^(3)_γb(k1_, k2_, k3_) = dF3/da^(3)_γb(k1, k2, k3, nu12, nu23, nu31)
                                         = 1/6 ( γ(k1_ + k2_, k3_) γ(k1_, k2_) + γ(k1_ + k3_, k2_) γ(k1_, k3_) + γ(k2_ + k3_, k1_) γ(k2_, k3_))
                                          - 1/6 ( α^(a)(k1_ + k2_, k3_) γ(k1_, k2_) + α^(a)(k1_ + k3_, k2_) γ(k1_, k3_) + α^(a)(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          + 1/6 ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          - 1/3 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

            dG3/da^(3)_γb(k1_, k2_, k3_) = dG3/da^(3)_γb(k1, k2, k3, nu12, nu23, nu31)
                                         = 0

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSin;

    double betaSum;
    double gammaSum;
    double alphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            alphaSum = kernels_alpha(kernVarH2);
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/da^(3)_γb(k1_, k2_, k3_) */
            dh3Kernels[0] += 1. / 6. * (gammaSum * gammaSin)
                            - 1. / 6. * (alphaSum * gammaSin)
                            + 1. / 6. * (betaSum * gammaSin)
                            - 1. / 3. * (gammaSum * betaSin);

            /* dG3/da^(3)_γb(k1_, k2_, k3_) */
            dh3Kernels[1] += 0.;

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh3_d3gaa(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to d^(3)_γa:

            dF3/dd^(3)_γa(k1_, k2_, k3_) = dF3/dd^(3)_γa(k1, k2, k3, nu12, nu23, nu31)
                                         = 0

            dG3/dd^(3)_γa(k1_, k2_, k3_) = dG3/dd^(3)_γa(k1, k2, k3, nu12, nu23, nu31)
                                         = 1/12 ( γ(k1_ + k2_, k3_) γ(k1_, k2_) + γ(k1_ + k3_, k2_) γ(k1_, k3_) + γ(k2_ + k3_, k1_) γ(k2_, k3_))
                                          + 1/12 ( α^(a)(k1_ + k2_, k3_) γ(k1_, k2_) + α^(a)(k1_ + k3_, k2_) γ(k1_, k3_) + α^(a)(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          - 1/12 ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          + 1/6 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSin;

    double betaSum;
    double gammaSum;
    double alphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            alphaSum = kernels_alpha(kernVarH2);
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/dd^(3)_γa(k1_, k2_, k3_) */
            dh3Kernels[0] += 0.;

            /* dG3/dd^(3)_γa(k1_, k2_, k3_) */
            dh3Kernels[1] += 1. / 12. * (gammaSum * gammaSin)
                            + 1. / 12. * (alphaSum * gammaSin)
                            - 1. / 12. * (betaSum * gammaSin)
                            + 1. / 6. * (gammaSum * betaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}


int kernels_btst_dh3_d3gab(kern_t *kernVar, double *dh3Kernels)
{
    /*

        Bootstrap third order density + velocity kernels derivative with respect to d^(3)_γb:

            dF3/dd^(3)_γb(k1_, k2_, k3_) = dF3/dd^(3)_γb(k1, k2, k3, nu12, nu23, nu31)
                                         = 0

            dG3/dd^(3)_γb(k1_, k2_, k3_) = dG3/dd^(3)_γb(k1, k2, k3, nu12, nu23, nu31)
                                         = 1/6 ( γ(k1_ + k2_, k3_) γ(k1_, k2_) + γ(k1_ + k3_, k2_) γ(k1_, k3_) + γ(k2_ + k3_, k1_) γ(k2_, k3_))
                                          - 1/6 ( α^(a)(k1_ + k2_, k3_) γ(k1_, k2_) + α^(a)(k1_ + k3_, k2_) γ(k1_, k3_) + α^(a)(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          + 1/6 ( β(k1_ + k2_, k3_) γ(k1_, k2_) + β(k1_ + k3_, k2_) γ(k1_, k3_) + β(k2_ + k3_, k1_) γ(k2_, k3_) )
                                          - 1/3 ( γ(k1_ + k2_, k3_) β(k1_, k2_) + γ(k1_ + k3_, k2_) β(k1_, k3_) + γ(k2_ + k3_, k1_) β(k2_, k3_) )

        Note: kernVar must contain all scales kij = |ki_ + kj_| and all angles nuij := nu(ki_, kj_) = ki_.kj_ / (ki kj), nuijk := nu(kij_, kk_) = (ki_ + kj_).kk_ / (kij kk) in the
              following order:

                    k[] = {k1, k2, k3}
                    nu[] = {nu12, nu13, nu23, nu123, nu132, nu231}

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* Second order kernel struct */
    kern_t *kernVarH2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarH3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Initialise kernels */
    dh3Kernels[0] = 0.;
    dh3Kernels[1] = 0.;

    /* Declare variables for basis kernels */
    double betaSin;
    double gammaSin;

    double betaSum;
    double gammaSum;
    double alphaSum;


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarH2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarH2 -> kernOrder; j++)
          {
            /* Ignore apprent singularits */
            if (kernVarH3 -> k[3 + index] <= __ABSTOL__)
              {
                index++;
                continue;
              }

            /*

                    β(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[i];
            kernVarH2 -> k[1] = kernVarH3 -> k[j + 1];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[index];

            /* Kernels */
            betaSin = kernels_beta(kernVarH2);
            gammaSin = kernels_gamma(kernVarH2);


            /*

                    β(k12_, k3_) + γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables */
            kernVarH2 -> k[0] = kernVarH3 -> k[3 + index];
            kernVarH2 -> k[1] = kernVarH3 -> k[2 - index];
            kernVarH2 -> nu[0] = kernVarH3 -> nu[3 + index];

            /* Kernels */
            betaSum = kernels_beta(kernVarH2);
            alphaSum = kernels_alpha(kernVarH2);
            gammaSum = kernels_gamma(kernVarH2);

            /* dF3/dd^(3)_γb(k1_, k2_, k3_) */
            dh3Kernels[0] += 0.;

            /* dG3/dd^(3)_γb(k1_, k2_, k3_) */
            dh3Kernels[1] += 1. / 6. * (gammaSum * gammaSin)
                            - 1. / 6. * (alphaSum * gammaSin)
                            + 1. / 6. * (betaSum * gammaSin)
                            - 1. / 3. * (gammaSum * betaSin);

            index++;
          }
      }

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return 0;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     BIASED KERNELS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */


/*  ------------------------------------------------------------------------------------------------------  */
/*  -----------------------------------   Phenomenological Smoothing   -----------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_smooth(kern_t *kernVar)
{
    /*

        This function calculates the phenomenological smoothing in the form:

                exp(-0.5 (ki_.s_)^2 sigv^2) * exp(-0.5 (ki_.s_)^2 (sigs0 (1 + z)/H(z))^2)
                \___________  ____________/   \___________________  ____________________/
                            \/                                    \/
                            FoG                             Spectra Errors

        where

                k_ = k1_ + ... + kn_  =>  (ki_.s_)^2 = (k1 mu1)^2 + ... + (kn mun)^2

        for the n'th order spectrum (i.e. n = 2 -> power spectrum, n = 3 -> bispectrum, ...)

        and

                sigs0 (1 + z)/H(z) := sigs(z)


    */

    double result = 0.;

    for (size_t i = 0; i < kernVar -> specOrder; i++)
      {
        result += (kernVar -> k[i] == 0.) ? 0. : pow(kernVar -> k[i] * kernVar -> mu[i], 2.);
      }

    result = exp(- 0.5 * result * pow(kernVar -> rsd -> sigv, 2.)) * exp(- 0.5 * result * pow(kernVar -> rsd -> sigs , 2.));

    return result;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dsmooth_sigv(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the phenomenological smoothing with respect to sigv:

                -(ki mui)^2 sigv * exp(-0.5 (ki mui)^2 sigv^2) * exp(-0.5 (ki mui)^2 (sigs0 (1 + z)/H(z))^2)
                                   \___________  ____________/   \___________________  ____________________/
                                               \/                                    \/
                                              FoG                              Spectra Errors


    */

    double result = 0.;

    for (size_t i = 0; i < kernVar -> specOrder; i++)
        result += (kernVar -> k[i] == 0.) ? 0. : pow(kernVar -> k[i] * kernVar -> mu[i], 2.);

    result = - result * kernVar -> rsd -> sigv * exp(- 0.5 * result * pow(kernVar -> rsd -> sigv, 2.)) * exp(- 0.5 * result * pow(kernVar -> rsd -> sigs , 2.));

    return result;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   First order Kernels   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_z1(kern_t *kernVar)
{
    /*

        This function calculates the kernel Z1:

                Z1(k_) = Z1(mu)
                       = b1 + f mu^2

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 1) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];


    /* Calculate kernels */

    /* Z1(k1_) */
    double z1Kernel = kernVar -> bias -> b1 + kernVar -> rsd -> f * (kernVarZ1 -> mu[0] * kernVarZ1 -> mu[0]);

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return z1Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz1_mu(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the kernel Z1 with respect to mu:

                dZ1/dmu(mu) = 2 * f * mu

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 1) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];


    /* Calculate kernels */

    /* dZ1/dmu1(k1_) */
    double dz1Kernel = 2. * kernVar -> rsd -> f * kernVarZ1 -> mu[0];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz1Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz1_b1(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the kernel Z1 with respect to b1:

                dZ1/db1(mu) = 1

    */

    (void) kernVar;

    return 1.;
}


double kernels_dz1_f(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the kernel Z1 with respect to f (growth rate):

                dZ1/df(mu) = mu^2

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 1) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];


    /* Calculate kernels */

    /* dZ1/df(k1_) */
    double dz1Kernel = kernVarZ1 -> mu[0]*kernVarZ1 -> mu[0];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz1Kernel;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  --------------------------------------   Second order Kernels   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_z2(kern_t *kernVar)
{
    /*

        This function calculates the second order biased density kernel:

                Z2(k1_, k2_) = Z2(k1, k2, nu12, mu1, mu2)
                             = b1 F2(k1, k2, nu12) + mu12^2 f G2(k1, k2, nu12)
                                + f k12 mu12 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2))
                                - ((b1 a^(2)_γ - c^(2)_γ) / 2) γ(nu12)
                                + b2/2


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    double h2Kernels[2];
    kernels_btst_h2(kernVar, h2Kernels);

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarZ2);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[0];
    double z1Kernel1 = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[1];
    double z1Kernel2 = kernels_z1(kernVar);

    /* Z2(k1_, k2_) */
    double z2Kernel = bias -> b1 * h2Kernels[0] + rsd -> f * (kernVarZ2 -> mu[2] * kernVarZ2 -> mu[2]) * h2Kernels[1]
                       + rsd -> f * kernVarZ2 -> mu[2] * kernVarZ2 -> k[2] / 2.
                          * (kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernel1 + kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernel2)
                       - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * gamma
                       + bias -> b2 / 2.;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return z2Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz2_k(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. k1 (first variable):

                dZ2/dk1(k1_, k2_) = dZ2/dk1(k1, k2, nu12, mu1, mu2)
                                  =  b1 dF2/dk1(k1, k2, nu12) + d(mu12^2 f G2(k1, k2, nu12)/dk1
                                     + f d(k12 mu12 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2)))/dk1
                                  =  b1 dF2/dk1(k1, k2, nu12) + 2 dmu12/dk1 mu12 f G2(k1, k2, nu12) + (mu12)^2 f dG2/dk1(k1, k2, nu12)
                                     + f d(k12 mu12)/dk1 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2))
                                     - f k12 mu12 / 2 * (mu1/k1^2 Z1(mu2))


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

                   =>  d(k12)/dk1 = (k1 + k2 nu12) / k12
                   =>  d(mu12)/dk1 = mu1 / k12 - mu12 / k12 d(k12)/dk1
                   =>  d(k12 mu12)/dk1 = mu1

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Additional variables */
    double dk12 = (kernVarZ2 -> k[0] + kernVarZ2 -> k[1] * kernVarZ2 -> nu[0]) / kernVarZ2 -> k[2];
    double dmu12 = (kernVarZ2 -> mu[0] - kernVarZ2 -> mu[2] * dk12) / kernVarZ2 -> k[2];


    /* Calculate kernels */

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    double h2Kernels[2];
    kernels_btst_h2(kernVar, h2Kernels);

    /* dF2/dk1(k1_, k2_) + dG2/dk1(k1_, k2_) */
    double dh2Kernels[2];
    kernels_btst_dh2_k(kernVar, dh2Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[0];
    double z1Kernel1 = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[1];
    double z1Kernel2 = kernels_z1(kernVar);

    /* dZ2/dk1(k1_, k2_) */
    double dz2Kernel = bias -> b1 * dh2Kernels[0] + rsd -> f * kernVarZ2 -> mu[2] * (2. * dmu12 * h2Kernels[1] + kernVarZ2 -> mu[2] * dh2Kernels[1])
                       + rsd -> f * kernVarZ2 -> mu[0] / 2.
                          * ( kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernel1
                                + (kernVarZ2 -> mu[0] - kernVarZ2 -> mu[2] * kernVarZ2 -> k[2] / kernVarZ2 -> k[0]) / kernVarZ2 -> k[0] * z1Kernel2 );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


double kernels_dz2_nu(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. nu12:

                dZ2/dnu12(k1_, k2_) = dZ2/dnu12(k1, k2, nu12, mu1, mu2)
                                    =  b1 dF2/dnu12(k1, k2, nu12) + d(mu12^2 f G2(k1, k2, nu12)/dnu12
                                        + f d(k12 mu12 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2)))/dnu12
                                        - ((b1 a^(2)_γ - c^(2)_γ) / 2) dγ/dnu12(nu12)
                                    =  b1 dF2/dnu12(k1, k2, nu12) + 2 dmu12/dnu12 mu12 f G2(k1, k2, nu12) + (mu12)^2 f dG2/dnu12(k1, k2, nu12)
                                        - ((b1 a^(2)_γ - c^(2)_γ) / 2) dγ/dnu12(nu12)


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

                   =>  d(k12)/dnu12 = k1 k2 / k12
                   =>  d(mu12)/dnu12 = - mu12 / k12 d(k12)/dk1
                   =>  d(k12 mu12)/dnu12 = 0

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct (not needed) */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Additional variables */
    double dk12 = kernVarZ2 -> k[0] * kernVarZ2 -> k[1] / kernVarZ2 -> k[2];
    double dmu12 = - kernVarZ2 -> mu[2] * dk12 / kernVarZ2 -> k[2];


    /* Calculate kernels */

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    double h2Kernels[2];
    kernels_btst_h2(kernVar, h2Kernels);

    /* dF2/dknu12(k1_, k2_) + dG2/dnu12(k1_, k2_) */
    double dh2Kernels[2];
    kernels_btst_dh2_nu(kernVar, dh2Kernels);

    /* dγ/dnu12(k1_, k2_) */
    double dgamma = kernels_dgamma_nu(kernVarZ2);

    /* dZ2/dnu12(k1_, k2_) */
    double dz2Kernel = bias -> b1 * dh2Kernels[0] + rsd -> f * kernVarZ2 -> mu[2] * (2. * dmu12 * h2Kernels[1] + kernVarZ2 -> mu[2] * dh2Kernels[1])
                        - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * dgamma;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


double kernels_dz2_mu(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. mu1 (first variable):

                dZ2/dmu1(k1_, k2_) = dZ2/dmu1(k1, k2, nu12, mu1, mu2)
                                   =  d(mu12^2)/dmu1 f G2(k1, k2, nu12)
                                       + f d(k12 mu12 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2)))/dmu1
                                   =  2 d(mu12)/dmu1 mu12 f G2(k1, k2, nu12)
                                       + f d(k12 mu12)/dmu1 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2))
                                       + f k12 mu12 / 2 * (mu2/k2 dZ1/dmu1(mu1) + 1/k1 Z1(mu2))


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

                   =>  d(k12)/dmu1 = 0
                   =>  d(mu12)/dmu1 = k1 / k12
                   =>  d(k12 mu12)/dk1 = k1

    */

    /* Fiducials */
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Additional variables */
    double dmu12 = kernVarZ2 -> k[0] / kernVarZ2 -> k[2];


    /* Calculate kernels */

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    double h2Kernels[2];
    kernels_btst_h2(kernVar, h2Kernels);

    /* Z1(k1_) + dZ1/dmu1(k1#_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[0];
    double z1Kernel1 = kernels_z1(kernVar);
    double dz1Kernel1 = kernels_dz1_mu(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[1];
    double z1Kernel2 = kernels_z1(kernVar);

    /* dZ2/dnu12(k1_, k2_) */
    double dz2Kernel = 2. * rsd -> f * kernVarZ2 -> mu[2] * dmu12 * h2Kernels[1]
                        + rsd -> f * kernVarZ2 -> k[0] / 2.
                          * ( kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernel1 + kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernel2 )
                        + rsd -> f * kernVarZ2 -> k[2] * kernVarZ2 -> mu[2] / 2.
                          * ( kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * dz1Kernel1 + 1. / kernVarZ2 -> k[0] * z1Kernel2 );


    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz2_a2ga(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. a^(2)_γ:

                dZ2/da^(2)_γ(k1_, k2_) = dZ2/da^(2)_γ(k1, k2, nu12, mu1, mu2)
                                       = b1 dF2/da^(2)_γ(k1, k2, nu12) + f (mu12)^2 dG2/da^(2)_γ(k1, k2, nu12)
                                          - b1 / 2 γ(nu12)


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct (not needed) */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* dF2/da^(2)_γ(k1_, k2_) + dG2/da^(2)_γ(k1_, k2_) */
    double dh2Kernels[2];
    kernels_btst_dh2_a2ga(kernVar, dh2Kernels);

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarZ2);

    /* dZ2/da^(2)_γ(k1_, k2_) */
    double dz2Kernel = bias -> b1 * dh2Kernels[0] + rsd -> f * (kernVarZ2 -> mu[2] * kernVarZ2 -> mu[2]) * dh2Kernels[1]
                       - bias -> b1 / 2. * gamma;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


double kernels_dz2_d2ga(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. d^(2)_γ:

                dZ2/dd^(2)_γ(k1_, k2_) = dZ2/dd^(2)_γ(k1, k2, nu12, mu1, mu2)
                                       = b1 dF2/dd^(2)_γ(k1, k2, nu12) + f (mu12)^2 dG2/dd^(2)_γ(k1, k2, nu12)

        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct (not needed) */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* dF2/dd^(2)_γ(k1_, k2_) + dG2/dd^(2)_γ(k1_, k2_) */
    double dh2Kernels[2];
    kernels_btst_dh2_d2ga(kernVar, dh2Kernels);

    /* dZ2/da^(2)_γ(k1_, k2_) */
    double dz2Kernel = bias -> b1 * dh2Kernels[0] + rsd -> f * (kernVarZ2 -> mu[2] * kernVarZ2 -> mu[2]) * dh2Kernels[1];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz2_b1(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. b1:

                dZ2/db1(k1_, k2_) = dZ2/db1(k1, k2, nu12, mu1, mu2)
                                  = F2(k1, k2, nu12)
                                     + f k12 mu12 / 2 * (mu2/k2 dZ1/db1(mu1) + mu1/k1 dZ1/db1(mu2))
                                     - a^(2)_γ / 2 γ(nu12)


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

    */

    /* Fiducials */
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    double h2Kernels[2];
    kernels_btst_h2(kernVar, h2Kernels);

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarZ2);

    /* dZ1/db1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[0];
    double dz1Kernel1 = kernels_dz1_b1(kernVar);

    /* dZ1/db1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[1];
    double dz1Kernel2 = kernels_dz1_b1(kernVar);

    /* dZ2/db1(k1_, k2_) */
    double dz2Kernel = h2Kernels[0]
                        + rsd -> f * kernVarZ2 -> mu[2] * kernVarZ2 -> k[2] / 2.
                          * (kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * dz1Kernel1 + kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * dz1Kernel2)
                        - btst -> a2Ga / 2. * gamma;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


double kernels_dz2_b2(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. b2:

                dZ2/db2(k1_, k2_) = dZ2/db2(k1, k2, nu12, mu1, mu2)
                                  = 1/2

    */

    (void) kernVar;

    return 0.5;
}


double kernels_dz2_f(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the second order biased density kernel w.r.t. f (growth rate):

                dZ2/df(k1_, k2_) = dZ2/df(k1, k2, nu12, mu1, mu2)
                                 = mu12^2 G2(k1, k2, nu12)
                                    + k12 mu12 / 2 * (mu2/k2 Z1(mu1) + mu1/k1 Z1(mu2))
                                    + f k12 mu12 / 2 * (mu2/k2 dZ1/df(mu1) + mu1/k1 dZ1/df(mu2))


        where

                k12_ = k1_ + k2_

                   =>  k12^2 = k1^2 + k2^2 + 2 k1 k2 nu12
                   =>  mu12 = (k1 mu1 + k2 mu2) / k12

    */

    /* Fiducials */
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* F2(k1_, k2_) + G2(k1_, k2_) */
    double h2Kernels[2];
    kernels_btst_h2(kernVar, h2Kernels);

    /* Z1(k1_) + dZ1/df(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[0];
    double z1Kernel1 = kernels_z1(kernVar);
    double dz1Kernel1 = kernels_dz1_f(kernVar);

    /* Z1(k2_) + dZ1/df(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ2 -> mu[1];
    double z1Kernel2 = kernels_z1(kernVar);
    double dz1Kernel2 = kernels_dz1_f(kernVar);

    /* dZ2/df(k1_, k2_) */
    double dz2Kernel = (kernVarZ2 -> mu[2] * kernVarZ2 -> mu[2]) * h2Kernels[1]
                        + kernVarZ2 -> mu[2] * kernVarZ2 -> k[2] / 2.
                          * (kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernel1 + kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernel2)
                        + rsd -> f * kernVarZ2 -> mu[2] * kernVarZ2 -> k[2] / 2.
                          * (kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * dz1Kernel1 + kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * dz1Kernel2);

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}


double kernels_dz2_c2ga(kern_t *kernVar)
{
    /*

        This function calculates the second order biased density kernel:

                dZ2/dc^(2)_γ(k1_, k2_) = dZ2/dc^(2)_γ(k1, k2, nu12, mu1, mu2)
                                       = γ(nu12) / 2


    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 2) : false;

    /* First order kernel struct (not needed) */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];


    /* Calculate kernels */

    /* γ(k1_, k2_) */
    double gamma = kernels_gamma(kernVarZ2);

    /* dZ2/dc^(2)_γ(k1_, k2_) */
    double dz2Kernel = gamma / 2.;

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz2Kernel;
}



/*  ------------------------------------------------------------------------------------------------------  */
/*  ---------------------------------------   Third order Kernels   --------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_z3(kern_t *kernVar)
{
    /*

        This function calculates the third order biased density kernel:

                Z3(k1_, k2_, k3_) = Z3(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                  = b1 F3(k1_, k2_, k3_) + mu123^2 f G3(k1_, k2_, k3_)
                                     + f mu123 k123 / 3 * (mu1/k1 Z2'(k2_, k3_) + mu2/k2 Z2'(k1_, k3_) + mu3/k3 Z2'(k1_, k2_))
                                     + f mu123 k123 / 3 * (mu12/k12 G2(k1_, k2_) Z1(k3_) + mu23/k23 G2(k2_, k3_) Z1(k1_) + mu31/k31 G2(k3_, k1_) Z1(k2_))
                                     + b2 / 3 * (F2(k1_, k2_) + F2(k2_, k3_) + F2(k3_, k1_))
                                     - (b1 a2Ga - c2Ga) / 3 * (γ(k12_, k3_) F2(k1_, k2_) + γ(k23_, k1_) F2(k2_, k3_) + γ(k31_, k2_) F2(k3_, k1_))
                                     - 2 bGam3 / 3 * ( γ(k12_, k3_) (F2(k1_, k2_) - G2(k1_, k2_)) + γ(k23_, k1_) (F2(k2_, k3_) - G2(k2_, k3_)) + γ(k31_, k2_) (F2(k3_, k1_) - G2(k3_, k1_)) )
                                     + b3 / 6

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2/2

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h2Kernels[3][2];

    double gammaSin[3];
    double gammaSum[3];

    double z1Kernels[3];

    double z2Kernels[3];


    /* Calculate kernels */

    /* F3(k1_, k2_, k3_) + G3(k1_, k2_, k3_) */
    double h3Kernels[2];
    kernels_btst_h3(kernVar, h3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            z2Kernels[index] = bias -> b1 * h2Kernels[index][0] + rsd -> f * (kernVarZ3 -> mu[index + 3] * kernVarZ3 -> mu[index + 3]) * h2Kernels[index][1]
                            + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1]
                                                                                       + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i] )
                            - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * gammaSin[index]
                            + bias -> b2 / 2.;


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);


            index++;
          }
      }

    /* Z3(k1_, k2_, k3_) */
    double z3Kernel = bias -> b1 * h3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * h3Kernels[1]

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * z2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * z2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * z2Kernels[2] )

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * h2Kernels[0][1] * z1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * h2Kernels[1][1] * z1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * z1Kernels[0] )

                       + bias -> b2 / 3. * ( h2Kernels[0][0] + h2Kernels[1][0] + h2Kernels[2][0] )

                       - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 3. * ( gammaSum[0] * h2Kernels[0][0]
                                                                         + gammaSum[1] * h2Kernels[1][0]
                                                                         + gammaSum[2] * h2Kernels[2][0] )

                       - 2. * bias -> bGam3 / 3. * ( gammaSum[0] * (h2Kernels[0][0] - h2Kernels[0][1])
                                                  + gammaSum[1] * (h2Kernels[1][0] - h2Kernels[1][1])
                                                  + gammaSum[2] * (h2Kernels[2][0] - h2Kernels[2][1]) );


    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return z3Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */
/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz3_k(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to k1 (first variable):

                dZ3/dk1(k1_, k2_, k3_) = dZ3/dk1(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                       = b1 dF3/dk1(k1_, k2_, k3_) + d(mu123^2 f G3(k1_, k2_, k3_))/dk1
                                          + f d(mu123 k123)/dk1 / 3 * (mu1/k1 Z2'(k2_, k3_) + mu2/k2 Z2'(k1_, k3_) + mu3/k3 Z2'(k1_, k2_))
                                              + f mu123 k123 / 3 * (-mu1/k1^2 Z2'(k2_, k3_) + mu2/k2 dZ2'/dk1(k1_, k3_) + mu3/k3 dZ2'/dk1(k1_, k2_))
                                          + f d(mu123 k123)/dk1 / 3 * (mu12/k12 G2(k1_, k2_) Z1(k3_) + mu23/k23 G2(k2_, k3_) Z1(k1_) + mu31/k31 G2(k3_, k1_) Z1(k2_))
                                              + f mu123 k123 / 3 * (d(mu12/k12 G2(k1_, k2_))/dk1 Z1(k3_) + d(mu31/k31 G2(k3_, k1_))/dk1 Z1(k2_))
                                          + b2 / 3 * (dF2(k1_, k2_)/dk1 + dF2(k1_, k3_)/dk1)
                                          - (b1 a2Ga - c2Ga) / 3 * (d(γ(k12_, k3_) F2(k1_, k2_))/dk1 + d(γ(k31_, k2_) F2(k3_, k1_))/dk1)
                                          - 2 bGam3 / 3 * ( d(γ(k12_, k3_) (F2(k1_, k2_) - G2(k1_, k2_)))/dk1 + d(γ(k31_, k2_) (F2(k3_, k1_) - G2(k3_, k1_)))/dk1 )

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2/2
                dZ2'/dk1(k1_, k2_) = b1 dF2/dk1(k1_, k2_) + d((mu12)^2 f G2(k1_, k2_))/dk1
                                        + f d(mu123 k123)/dk1 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - f mu123 k123 / 4 * mu1/k1^2 Z1(k2_)

        and

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

                    =>  dk1i/k1 = (k1 + ki nu1i) / k1i
                    =>  dmu1i/dk1 = mu1 / k1i - mu1i/k1i dk1i/dk1
                    =>  d(mu1i k1i) = mu1

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

                    =>  dk123/dk1 = (k1 + k2 nu12 + k3 nu13) / k123
                    =>  dmu123/dk1 = mu1 / k123 - mu123/k123 dk123/dk1
                    =>  d(mu123 k123)/dk1 = mu1

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h3Kernels[2];
    double dh3Kernels[2];

    double h2Kernels[3][2];
    double dh2Kernels[3][2];

    double gammaSin[3];
    double gammaSum[3];
    double dgammaSum[3];

    double z1Kernels[3];

    double z2Kernels[3];
    double dz2Kernels[3];

    /* Additional variables */
    double dk[3];
    double dmu[3];
    double dnu[3];

    double dk123 = (kernVarZ3 -> k[0] + kernVarZ3 -> k[1] * kernVarZ3 -> nu[0] + kernVarZ3 -> k[2] * kernVarZ3 -> nu[1]) / kernVarZ3 -> k[6];
    double dmu123 = (kernVarZ3 -> mu[0] - kernVarZ3 -> mu[6] * dk123) / kernVarZ3 -> k[6];


    /* Calculate kernels */

    /* F3(k1_, k2_, k3_) + G3(k1_, k2_, k3_) */
    kernels_btst_h3(kernVar, h3Kernels);

    /* dF3/dk1(k1_, k2_, k3_) + dG3/dk1(k1_, k2_, k3_) */
    kernels_btst_dh3_k(kernVar, dh3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            dk[index] = (i == 0) ? (kernVarZ3 -> k[i] + kernVarZ3 -> k[j + 1] * kernVarZ3 -> nu[j]) / kernVarZ3 -> k[index + 3] : 0.;
            dmu[index] = (i == 0) ? (kernVarZ3 -> mu[i] - kernVarZ3 -> mu[index + 3] * dk[index]) / kernVarZ3 -> k[index + 3] : 0.;
            dnu[index] = (i == 0) ? (kernVarZ3 -> nu[1 - j] - kernVarZ3 -> nu[index + 3] * dk[index]) / kernVarZ3 -> k[index + 3] : 0.;

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            if (i == 0) kernels_btst_dh2_k(kernVar, dh2Kernels[index]);
            else {dh2Kernels[index][0] = 0; dh2Kernels[index][1] = 0;}

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            z2Kernels[index] = bias -> b1 * h2Kernels[index][0] + rsd -> f * (kernVarZ3 -> mu[3 + index]*kernVarZ3 -> mu[3 + index]) * h2Kernels[index][1]
                                + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * (kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1]
                                                                                                + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i])
                                - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * gammaSin[index]
                                + bias -> b2 / 2.;

            dz2Kernels[index] = rsd -> f * kernVarZ3 -> mu[0] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1] + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i] );
            dz2Kernels[index] += (i == 0) ? bias -> b1 * dh2Kernels[index][0]
                                                + rsd -> f * kernVarZ3 -> mu[index + 3] * (2. * dmu[index] * h2Kernels[index][1] + kernVarZ3 -> mu[index + 3] * dh2Kernels[index][1])
                                                - rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * kernVarZ2 -> mu[0] / (kernVarZ2 -> k[0]*kernVarZ2 -> k[0]) * z1Kernels[j + 1]
                                          : 0.;


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);
            dgammaSum[index] = kernels_dgamma_nu(kernVarZ2) * dnu[index];

            index++;
          }
      }

    /* dZ3/dk1(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * kernVarZ3 -> mu[6] * (2. * dmu123 * h3Kernels[1] + kernVarZ3 -> mu[6] * dh3Kernels[1])

               + rsd -> f * kernVarZ3 -> mu[0] / 3.
                  * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * z2Kernels[0]
                    + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * z2Kernels[1]
                    + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * z2Kernels[2] )
               + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                  * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                    + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                    + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] - kernVarZ3 -> mu[0] / pow(kernVarZ3 -> k[0], 2) * z2Kernels[2] )

               + rsd -> f * kernVarZ3 -> mu[0] / 3.
                  * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * h2Kernels[0][1] * z1Kernels[2]
                    + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * h2Kernels[1][1] * z1Kernels[1]
                    + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * z1Kernels[0] )
               + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                  * ( ((dmu[0] - kernVarZ3 -> mu[3] * dk[0] / kernVarZ3 -> k[3]) / kernVarZ3 -> k[3] * h2Kernels[0][1] + kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * dh2Kernels[0][1]) * z1Kernels[2]
                    + ((dmu[1] - kernVarZ3 -> mu[4] * dk[1] / kernVarZ3 -> k[4]) / kernVarZ3 -> k[4] * h2Kernels[1][1] + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * dh2Kernels[1][1]) * z1Kernels[1] )

               + bias -> b2 / 3. * ( dh2Kernels[0][0] + dh2Kernels[1][0] )

               - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 3. * ( dgammaSum[0] * h2Kernels[0][0] + gammaSum[0] * dh2Kernels[0][0]
                                                                 + dgammaSum[1] * h2Kernels[1][0] + gammaSum[1] * dh2Kernels[1][0])

               - 2. * bias -> bGam3 / 3. * ( dgammaSum[0] * (h2Kernels[0][0] - h2Kernels[0][1]) + gammaSum[0] * (dh2Kernels[0][0] - dh2Kernels[0][1])
                                          + dgammaSum[1] * (h2Kernels[1][0] - h2Kernels[1][1]) + gammaSum[1] * (dh2Kernels[1][0] - dh2Kernels[1][1]) );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_nu(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to nu12 (first variable):

                dZ3/dnu12(k1_, k2_, k3_) = dZ3/dnu12(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                         = b1 dF3/dnu12(k1_, k2_, k3_) + d(mu123^2 f G3(k1_, k2_, k3_))/dnu12
                                            + f mu123 k123 / 3 * mu3/k3 dZ2'/dnu12(k1_, k2_)
                                            + f mu123 k123 / 3 * d(mu12/k12 G2(k1_, k2_))/dnu12 Z1(k3_)
                                            + b2 / 3 * (dF2(k1_, k2_)/dnu12)
                                            - (b1 a2Ga - c2Ga) / 3 * ( d(γ(k12_, k3_) F2(k1_, k2_))/dnu12 + dγ/nu12(k13_, k2_) F2(k1_, k3_) + dγ/nu12(k23_, k1_) F2(k2_, k3_) )
                                            - 2 bGam3 / 3 * ( d(γ(k12_, k3_) (F2(k1_, k2_) - G2(k1_, k2_)))/dnu12 + dγ/nu12(k13_, k2_) (F2(k1_, k3_) - G2(k1_, k3_))
                                                            + dγ/nu12(k23_, k1_) (F2(k2_, k3_) - G2(k2_, k3_)) )

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2
                dZ2'/dnu12(k1_, k2_) = b1 dF2/dnu12(k1_, k2_) + d((mu12)^2 f G2(k1_, k2_))/dnu12
                                        - (b1 a2Ga - c2Ga) / 2 dγ/dnu12(k1_, k2_)

        and

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

                    =>  dk12/nu12 = k1 k2 / k12
                    =>  dmu12/dnu12 = - mu12/k12 dk12/dnu12
                    =>  d(mu1i k1i) = 0

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

                    =>  dk123/nu12 = k1 k2 / k123
                    =>  dmu123/dnu12 = - mu123/k123 dk123/dnu12
                    =>  d(mu123 k123) = 0

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h3Kernels[2];
    double dh3Kernels[2];

    double h2Kernels[3][2];
    double dh2Kernels[3][2];

    double dgammaSin[3];
    double gammaSum[3];
    double dgammaSum[3];

    double z1Kernels[3];

    double dz2Kernels[3];

    /* Additional variables */
    double dk[3];
    double dmu[3];
    double dnu[3];

    double dk123 = kernVarZ3 -> k[0] * kernVarZ3 -> k[1] / kernVarZ3 -> k[6];
    double dmu123 = - kernVarZ3 -> mu[6] * dk123 / kernVarZ3 -> k[6];


    /* Calculate kernels */

    /* F3(k1_, k2_, k3_) + G3(k1_, k2_, k3_) */
    kernels_btst_h3(kernVar, h3Kernels);

    /* dF3/dnu12(k1_, k2_, k3_) + dG3/dnu12(k1_, k2_, k3_) */
    kernels_btst_dh3_nu(kernVar, dh3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            dk[index] = (j == 0) ? kernVarZ3 -> k[i] * kernVarZ3 -> k[j + 1] / kernVarZ3 -> k[index + 3] : 0.;
            dmu[index] = (j == 0) ? - kernVarZ3 -> mu[index + 3] * dk[index] / kernVarZ3 -> k[index + 3] : 0.;
            dnu[index] = (j == 0) ? - kernVarZ3 -> nu[index + 3] * dk[index] / kernVarZ3 -> k[index + 3] : kernVarZ3 -> k[i] / kernVarZ3 -> k[index + 3];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            if (j == 0) kernels_btst_dh2_nu(kernVar, dh2Kernels[index]);
            else {dh2Kernels[index][0] = 0; dh2Kernels[index][1] = 0;}

            dgammaSin[index] = (j == 0) ? kernels_dgamma_nu(kernVarZ2) : 0.;

            /* Z2'(k1_, k2_; k3_) */
            dz2Kernels[index] = (j == 0) ? bias -> b1 * dh2Kernels[index][0]
                                            + rsd -> f * kernVarZ3 -> mu[index + 3] * (2. * dmu[index] * h2Kernels[index][1] + kernVarZ3 -> mu[index + 3] * dh2Kernels[index][1])
                                            - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * dgammaSin[index]
                                          : 0.;


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);
            dgammaSum[index] = kernels_dgamma_nu(kernVarZ2) * dnu[index];

            index++;
          }
      }

    /* dZ3/dnu12(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * kernVarZ3 -> mu[6] * (2. * dmu123 * h3Kernels[1] + kernVarZ3 -> mu[6] * dh3Kernels[1])

               + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                  * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0] )

               + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                  * ( ((dmu[0] - kernVarZ3 -> mu[3] * dk[0] / kernVarZ3 -> k[3]) / kernVarZ3 -> k[3] * h2Kernels[0][1] + kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * dh2Kernels[0][1]) * z1Kernels[2] )

               + bias -> b2 / 3. * ( dh2Kernels[0][0] )

               - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 3. * ( dgammaSum[0] * h2Kernels[0][0] + gammaSum[0] * dh2Kernels[0][0]
                                                                     + dgammaSum[1] * h2Kernels[1][0]
                                                                     + dgammaSum[2] * h2Kernels[2][0] )

               - 2. * bias -> bGam3 / 3. * ( dgammaSum[0] * (h2Kernels[0][0] - h2Kernels[0][1]) + gammaSum[0] * (dh2Kernels[0][0] - dh2Kernels[0][1])
                                          + dgammaSum[1] * (h2Kernels[1][0] - h2Kernels[1][1])
                                          + dgammaSum[2] * (h2Kernels[2][0] - h2Kernels[2][1]) );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_mu(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to mu1 (first variable):

                dZ3/dmu1(k1_, k2_, k3_) = dZ3/dmu1(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                        = d(mu123^2)/dmu1 f G3(k1_, k2_, k3_)
                                           + f d(mu123 k123)/dmu1 / 3 * (mu1/k1 Z2'(k2_, k3_) + mu2/k2 Z2'(k1_, k3_) + mu3/k3 Z2'(k1_, k2_))
                                               + f mu123 k123 / 3 * (1/k1 Z2'(k2_, k3_) + mu2/k2 dZ2'/dmu1(k1_, k3_) + mu3/k3 dZ2'/dmu1(k1_, k2_))
                                           + f d(mu123 k123)/dmu1 / 3 * (mu12/k12 G2(k1_, k2_) Z1(k3_) + mu23/k23 G2(k2_, k3_) Z1(k1_) + mu31/k31 G2(k3_, k1_) Z1(k2_))
                                               + f mu123 k123 / 3 * (d(mu12/k12)/dmu1 G2(k1_, k2_) Z1(k3_) + d(mu31/k31)/dmu1 G2(k3_, k1_) Z1(k2_) + mu23/k23 G2(k2_, k3_) dZ1/dmu1(k1_))

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2
                dZ2'/dmu1(k1_, k2_) = d((mu12)^2)/dmu1 f G2(k1_, k2_)
                                        + f d(mu123 k123)/dmu1 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_))
                                        + f mu123 k123 / 4 * (1/k1 Z1(k2_) + mu2/k2 dZ1/dmu1(k1_))

        and

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

                    =>  dmu1j/dmu1 = k1 / k1j

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

                    =>  dmu123/dmu1 = k1 / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h3Kernels[2];

    double h2Kernels[3][2];

    double gammaSin[3];

    double z1Kernels[3];
    double dz1Kernels[3];

    double z2Kernels[3];
    double dz2Kernels[3];

    /* Additional variables */
    double dmu[3];

    double dmu123 = kernVarZ3 -> k[0] / kernVarZ3 -> k[6];


    /* Calculate kernels */

    /* F3(k1_, k2_, k3_) + G3(k1_, k2_, k3_) */
    kernels_btst_h3(kernVar, h3Kernels);

    /* Z1(k1_) + dZ1/dmu1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);
    dz1Kernels[0] = kernels_dz1_mu(kernVar);

    /* Z1(k2_) + dZ1/dmu1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);
    dz1Kernels[1] = 0.;

    /* Z1(k3_) + dZ1/dmu1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);
    dz1Kernels[2] = 0.;


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            dmu[index] = kernVarZ2 -> k[0] / kernVarZ3 -> k[index + 3];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            z2Kernels[index] = bias -> b1 * h2Kernels[index][0] + rsd -> f * (kernVarZ3 -> mu[index + 3]*kernVarZ3 -> mu[index + 3]) * h2Kernels[index][1]
                                + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * (kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1]
                                                                                                + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i])
                                - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * gammaSin[index]
                                + bias -> b2 / 2.;


            dz2Kernels[index] = rsd -> f * kernVarZ3 -> k[i] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1] + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i] );
            dz2Kernels[index] += (i == 0) ? 2. * rsd -> f * kernVarZ3 -> mu[index + 3] * dmu[index] * h2Kernels[index][1]
                                                + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * ( 1. / kernVarZ2 -> k[0] * z1Kernels[j + 1]
                                                                                                             + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * dz1Kernels[i] )
                                          : 0.;


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            index++;
          }
      }

    /* dZ3/dmu1(k1_, k2_, k3_) */
    double dz3Kernel = 2. * rsd -> f * kernVarZ3 -> mu[6] * dmu123 * h3Kernels[1]

                       + rsd -> f * kernVarZ3 -> k[0] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * z2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * z2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * z2Kernels[2] )
                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] + 1. / kernVarZ3 -> k[0] * z2Kernels[2] )

                       + rsd -> f * kernVarZ3 -> k[0] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * h2Kernels[0][1] * z1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * h2Kernels[1][1] * z1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * z1Kernels[0] )
                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( dmu[0] / kernVarZ3 -> k[3] * h2Kernels[0][1] * z1Kernels[2]
                            + dmu[1] / kernVarZ3 -> k[4] * h2Kernels[1][1] * z1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * dz1Kernels[0] );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz3_a2ga(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to a^(2)_γ:

                dZ3/d(a^(2)_γ)(k1_, k2_, k3_) = dZ3/d(a^(2)_γ)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(a^(2)_γ)(k1_, k2_, k3_) + mu123^2 f dG3/d(a^(2)_γ)(k1_, k2_, k3_)
                                                 + f mu123 k123 / 3 * (mu1/k1 dZ2'/d(a^(2)_γ)(k2_, k3_) + mu2/k2 dZ2'/d(a^(2)_γ)(k1_, k3_) + mu3/k3 dZ2'/d(a^(2)_γ)(k1_, k2_))
                                                 + f mu123 k123 / 3 * (mu12/k12 dG2/d(a^(2)_γ)(k1_, k2_) Z1(k3_) + mu23/k23 dG2/d(a^(2)_γ)(k2_, k3_) Z1(k1_) + mu31/k31 dG2/d(a^(2)_γ)(k3_, k1_) Z1(k2_))
                                                 - b1 / 3 * (γ(k12_, k3_) F2(k1_, k2_) + γ(k23_, k1_) F2(k2_, k3_) + γ(k31_, k2_) F2(k3_, k1_))
                                                     - (b1 a2Ga - c2Ga) / 3 * (γ(k12_, k3_) dF2/d(a^(2)_γ)(k1_, k2_) + γ(k23_, k1_) dF2/d(a^(2)_γ)(k2_, k3_) + γ(k31_, k2_) dF2/d(a^(2)_γ)(k3_, k1_))
                                                 - 2 bGam3 / 3 * ( γ(k12_, k3_) (dF2/d(a^(2)_γ)(k1_, k2_) - dG2/d(a^(2)_γ)(k1_, k2_))
                                                                 + γ(k23_, k1_) (dF2/d(a^(2)_γ)(k2_, k3_) - dG2/d(a^(2)_γ)(k2_, k3_))
                                                                 + γ(k31_, k2_) (dF2/d(a^(2)_γ)(k3_, k1_) - dG2/d(a^(2)_γ)(k3_, k1_)) )

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2
                dZ2'/d(a^(2)_γ)(k1_, k2_) = b1 dF2/d(a^(2)_γ)(k1_, k2_) + (mu12)^2 f dG2/d(a^(2)_γ)(k1_, k2_) - b1 / 2 γ(k1_, k2_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h2Kernels[3][2];
    double dh2Kernels[3][2];

    double gammaSin[3];
    double gammaSum[3];

    double z1Kernels[3];

    double dz2Kernels[3];


    /* Calculate kernels */

    /* dF3/da^(2)_γ(k1_, k2_, k3_) + dG3/da^(2)_γ(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_a2ga(kernVar, dh3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);
            kernels_btst_dh2_a2ga(kernVar, dh2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            dz2Kernels[index] = bias -> b1 * dh2Kernels[index][0] + rsd -> f * (kernVarZ3 -> mu[index + 3] * kernVarZ3 -> mu[index + 3]) * dh2Kernels[index][1]
                                - bias -> b1 / 2. * gammaSin[index];


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);


            index++;
          }
      }

    /* dZ3/(a^(2)_γ)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1]

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] )

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * dh2Kernels[0][1] * z1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * dh2Kernels[1][1] * z1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * dh2Kernels[2][1] * z1Kernels[0] )

                       - bias -> b1 / 3. * ( gammaSum[0] * h2Kernels[0][0]
                                          + gammaSum[1] * h2Kernels[1][0]
                                          + gammaSum[2] * h2Kernels[2][0] )

                       - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 3. * ( gammaSum[0] * dh2Kernels[0][0]
                                                                         + gammaSum[1] * dh2Kernels[1][0]
                                                                         + gammaSum[2] * dh2Kernels[2][0] )

                       - 2. * bias -> bGam3 / 3. * ( gammaSum[0] * (dh2Kernels[0][0] - dh2Kernels[0][1])
                                                  + gammaSum[1] * (dh2Kernels[1][0] - dh2Kernels[1][1])
                                                  + gammaSum[2] * (dh2Kernels[2][0] - dh2Kernels[2][1]) );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_d2ga(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to d^(2)_γ:

                dZ3/d(d^(2)_γ)(k1_, k2_, k3_) = dZ3/d(d^(2)_γ)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(d^(2)_γ)(k1_, k2_, k3_) + mu123^2 f dG3/d(d^(2)_γ)(k1_, k2_, k3_)
                                                 + f mu123 k123 / 3 * (mu1/k1 dZ2'/d(d^(2)_γ)(k2_, k3_) + mu2/k2 dZ2'/d(d^(2)_γ)(k1_, k3_) + mu3/k3 dZ2'/d(d^(2)_γ)(k1_, k2_))
                                                 + f mu123 k123 / 3 * (mu12/k12 dG2/d(d^(2)_γ)(k1_, k2_) Z1(k3_) + mu23/k23 dG2/d(d^(2)_γ)(k2_, k3_) Z1(k1_) + mu31/k31 dG2/d(d^(2)_γ)(k3_, k1_) Z1(k2_))
                                                 - (b1 a2Ga - c2Ga) / 3 * (γ(k12_, k3_) dF2/d(d^(2)_γ)(k1_, k2_) + γ(k23_, k1_) dF2/d(d^(2)_γ)(k2_, k3_) + γ(k31_, k2_) dF2/d(d^(2)_γ)(k3_, k1_))
                                                 - 2 bGam3 / 3 * ( γ(k12_, k3_) (dF2/d(d^(2)_γ)(k1_, k2_) - dG2/d(d^(2)_γ)(k1_, k2_))
                                                                 + γ(k23_, k1_) (dF2/d(d^(2)_γ)(k2_, k3_) - dG2/d(d^(2)_γ)(k2_, k3_))
                                                                 + γ(k31_, k2_) (dF2/d(d^(2)_γ)(k3_, k1_) - dG2/d(d^(2)_γ)(k3_, k1_)) )

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2.
                dZ2'/d(d^(2)_γ)(k1_, k2_) = b1 dF2/d(d^(2)_γ)(k1_, k2_) + (mu12)^2 f dG2/d(d^(2)_γ)(k1_, k2_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h2Kernels[3][2];
    double dh2Kernels[3][2];

    double gammaSin[3];
    double gammaSum[3];

    double z1Kernels[3];

    double dz2Kernels[3];


    /* Calculate kernels */

    /* dF3/dd^(2)_γ(k1_, k2_, k3_) + dG3/dd^(2)_γ(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_d2ga(kernVar, dh3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);
            kernels_btst_dh2_d2ga(kernVar, dh2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            dz2Kernels[index] = bias -> b1 * dh2Kernels[index][0] + rsd -> f * (kernVarZ3 -> mu[index + 3] * kernVarZ3 -> mu[index + 3]) * dh2Kernels[index][1]
                                - bias -> b1 / 2. * gammaSin[index];


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);


            index++;
          }
      }

    /* dZ3/(d^(2)_γ)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1]

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] )

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * dh2Kernels[0][1] * z1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * dh2Kernels[1][1] * z1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * dh2Kernels[2][1] * z1Kernels[0] )

                       - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 3. * ( gammaSum[0] * dh2Kernels[0][0]
                                                                         + gammaSum[1] * dh2Kernels[1][0]
                                                                         + gammaSum[2] * dh2Kernels[2][0] )

                       - 2. * bias -> bGam3 / 3. * ( gammaSum[0] * (dh2Kernels[0][0] - dh2Kernels[0][1])
                                                  + gammaSum[1] * (dh2Kernels[1][0] - dh2Kernels[1][1])
                                                  + gammaSum[2] * (dh2Kernels[2][0] - dh2Kernels[2][1]) );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_h(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to h:

                dZ3/d(h)(k1_, k2_, k3_) = dZ3/d(h)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(h)(k1_, k2_, k3_) + mu123^2 f dG3/d(h)(k1_, k2_, k3_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed9 */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct (not needed) */
//    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];


    /* Calculate kernels */

    /* dF3/dh(k1_, k2_, k3_) + dG3/dh(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_h(kernVar, dh3Kernels);

    /* dZ3/(h)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz3_a3gaa(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to a^(3)_γa:

                dZ3/d(a^(3)_γa)(k1_, k2_, k3_) = dZ3/d(a^(3)_γa)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(a^(3)_γa)(k1_, k2_, k3_) + mu123^2 f dG3/d(a^(3)_γa)(k1_, k2_, k3_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed9 */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct (not needed) */
//    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];


    /* Calculate kernels */

    /* dF3/d(a^(3)_γa)(k1_, k2_, k3_) + dG3/d(a^(3)_γa)(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_a3gaa(kernVar, dh3Kernels);

    /* dZ3/(a^(3)_γa)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_d3gaa(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to d^(3)_γa:

                dZ3/d(d^(3)_γa)(k1_, k2_, k3_) = dZ3/d(d^(3)_γa)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(d^(3)_γa)(k1_, k2_, k3_) + mu123^2 f dG3/d(d^(3)_γa)(k1_, k2_, k3_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed9 */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct (not needed) */
//    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];


    /* Calculate kernels */

    /* dF3/d(d^(3)_γa)(k1_, k2_, k3_) + dG3/d(d^(3)_γa)(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_d3gaa(kernVar, dh3Kernels);

    /* dZ3/(d^(3)_γa)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_a3gab(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to a^(3)_γb:

                dZ3/d(a^(3)_γb)(k1_, k2_, k3_) = dZ3/d(a^(3)_γb)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(a^(3)_γb)(k1_, k2_, k3_) + mu123^2 f dG3/d(a^(3)_γb)(k1_, k2_, k3_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */


    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed9 */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct (not needed) */
//    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];


    /* Calculate kernels */

    /* dF3/d(a^(3)_γb)(k1_, k2_, k3_) + dG3/d(a^(3)_γb)(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_a3gab(kernVar, dh3Kernels);

    /* dZ3/(a^(3)_γb)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_d3gab(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to d^(3)_γb:

                dZ3/d(d^(3)_γb)(k1_, k2_, k3_) = dZ3/d(d^(3)_γb)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                              = b1 dF3/d(d^(3)_γb)(k1_, k2_, k3_) + mu123^2 f dG3/d(d^(3)_γb)(k1_, k2_, k3_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed9 */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct (not needed) */
//    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];


    /* Calculate kernels */

    /* dF3/d(d^(3)_γb)(k1_, k2_, k3_) + dG3/d(d^(3)_γb)(k1_, k2_, k3_) */
    double dh3Kernels[2];
    kernels_btst_dh3_d3gab(kernVar, dh3Kernels);

    /* dZ3/(d^(3)_γb)(k1_, k2_, k3_) */
    double dz3Kernel = bias -> b1 * dh3Kernels[0] + rsd -> f * (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * dh3Kernels[1];

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


/*  ------------------------------------------------------------------------------------------------------  */


double kernels_dz3_b1(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to b1:

                dZ3/db1(k1_, k2_, k3_) = dZ3/db1(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                       = F3(k1_, k2_, k3_)
                                          + f mu123 k123 / 3 * (mu1/k1 dZ2'/db1(k2_, k3_) + mu2/k2 dZ2'/db1(k1_, k3_) + mu3/k3 dZ2'/db1(k1_, k2_))
                                          + f mu123 k123 / 3 * (mu12/k12 G2(k1_, k2_) dZ1/db1(k3_) + mu23/k23 G2(k2_, k3_) dZ1/db1(k1_) + mu31/k31 G2(k3_, k1_) dZ1/db1(k2_))
                                          - a^(2)_γ / 3 * (γ(k12_, k3_) F2(k1_, k2_) + γ(k23_, k1_) F2(k2_, k3_) + γ(k31_, k2_) F2(k3_, k1_))

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2
                dZ2'/db1(k1_, k2_) = F2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 dZ1/db1(k2_) + mu2/k2 dZ1/db1(k1_)) - a2Ga / 2 γ(k1_, k2_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h3Kernels[2];

    double h2Kernels[3][2];

    double gammaSin[3];
    double gammaSum[3];

    double dz1Kernels[3];

    double dz2Kernels[3];


    /* Calculate kernels */

    /* F3(k1_, k2_, k3_) + G3(k1_, k2_, k3_) */
    kernels_btst_h3(kernVar, h3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    dz1Kernels[0] = kernels_dz1_b1(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    dz1Kernels[1] = kernels_dz1_b1(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    dz1Kernels[2] = kernels_dz1_b1(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            dz2Kernels[index] = h2Kernels[index][0]
                                + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * dz1Kernels[j + 1]
                                                                                           + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * dz1Kernels[i] )
                                - btst -> a2Ga / 2. * gammaSin[index];


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);

            index++;
          }
      }

    /* dZ3/db1(k1_, k2_, k3_) */
    double dz3Kernel = h3Kernels[0]

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] )

                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * h2Kernels[0][1] * dz1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * h2Kernels[1][1] * dz1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * dz1Kernels[0] )

                       - btst -> a2Ga / 3. * ( gammaSum[0] * h2Kernels[0][0]
                                                 + gammaSum[1] * h2Kernels[1][0]
                                                 + gammaSum[2] * h2Kernels[2][0] );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_c2ga(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to c^(2)_γ:

                dZ3/d(c^(2)_γ)(k1_, k2_, k3_) = dZ3/d(c^(2)_γ)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                       = f mu123 k123 / 3 * (mu1/k1 dZ2'/d(c^(2)_γ)(k2_, k3_) + mu2/k2 dZ2'/d(c^(2)_γ)(k1_, k3_) + mu3/k3 dZ2'/d(c^(2)_γ)(k1_, k2_))
                                          + 1 / 3 * (γ(k12_, k3_) F2(k1_, k2_) + γ(k23_, k1_) F2(k2_, k3_) + γ(k31_, k2_) F2(k3_, k1_))

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2
                dZ2'/d(c^(2)_γ)(k1_, k2_) = 1 / 2 γ(k1_, k2_)

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_rsd_t *rsd = kernVar -> rsd;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed) */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h2Kernels[3][2];

    double gammaSin[3];
    double gammaSum[3];

    double dz2Kernels[3];


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            dz2Kernels[index] = 1. / 2. * gammaSin[index];


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);

            index++;
          }
      }

    /* dZ3/d(c^(2)_γ)(k1_, k2_, k3_) */
    double dz3Kernel = rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                         * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                           + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                           + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] )

                     + 1. / 3. * ( gammaSum[0] * h2Kernels[0][0]
                                 + gammaSum[1] * h2Kernels[1][0]
                                 + gammaSum[2] * h2Kernels[2][0] );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_bgam3(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to b_Γ3:

                dZ3/d(b_Γ3)(k1_, k2_, k3_) = dZ3/d(b_Γ3)(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                           = - 2 / 3 * ( γ(k12_, k3_) (F2(k1_, k2_) - G2(k1_, k2_)) + γ(k23_, k1_) (F2(k2_, k3_) - G2(k2_, k3_)) + γ(k31_, k2_) (F2(k3_, k1_) - G2(k3_, k1_)) )

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

    */

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct (not needed) */
//    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h2Kernels[3][2];

    double gammaSum[3];


    /* Calculate kernels */

    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);


            /*

                    γ(k12_, k3_)

                for k12, k3, nu12,3 given in the following order: (k[3],k[2],nu[3]) -> (k[4],k[1],nu[4]) -> (k[5],k[0],nu[5])

            */

            /* Variables (only need nu12,3 for γ(k12_, k3_)) */
            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[3 + index];

            /* Kernels */
            gammaSum[index] = kernels_gamma(kernVarZ2);

            index++;
          }
      }

    /* dZ3/d(b_Γ3)(k1_, k2_, k3_) */
    double dz3Kernel = -2. / 3. * ( gammaSum[0] * (h2Kernels[0][0] - h2Kernels[0][1])
                                  + gammaSum[1] * (h2Kernels[1][0] - h2Kernels[1][1])
                                  + gammaSum[2] * (h2Kernels[2][0] - h2Kernels[2][1]) );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}


double kernels_dz3_f(kern_t *kernVar)
{
    /*

        This function calculates the derivative of the third order biased density kernel with respect to f:

                dZ37df(k1_, k2_, k3_) = dZ3/df(k1, k2, k3, nu12, nu23, nu31, mu1, mu2, mu3)
                                      = mu123^2 G3(k1_, k2_, k3_)
                                         + mu123 k123 / 3 * (mu1/k1 Z2'(k2_, k3_) + mu2/k2 Z2'(k1_, k3_) + mu3/k3 Z2'(k1_, k2_))
                                             + f mu123 k123 / 3 * (mu1/k1 dZ2'/df(k2_, k3_) + mu2/k2 dZ2'/df(k1_, k3_) + mu3/k3 dZ2'/df(k1_, k2_))
                                         + mu123 k123 / 3 * (mu12/k12 G2(k1_, k2_) Z1(k3_) + mu23/k23 G2(k2_, k3_) Z1(k1_) + mu31/k31 G2(k3_, k1_) Z1(k2_))
                                             + f mu123 k123 / 3 * (mu12/k12 G2(k1_, k2_) dZ1(k3_)/df + mu23/k23 G2(k2_, k3_) dZ1(k1_)/df + mu31/k31 G2(k3_, k1_) dZ1(k2_)/df)

        where

                Z2'(k1_, k2_) = b1 F2(k1_, k2_) + (mu12)^2 f G2(k1_, k2_) + f mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) - (b1 a2Ga - c2Ga) / 2 γ(k1_, k2_) + b2 / 2
                dZ2'/df(k1_, k2_) = (mu12)^2 G2(k1_, k2_) + mu123 k123 / 4 * (mu1/k1 Z1(k2_) + mu2/k2 Z1(k1_)) + f mu123 k123 / 4 * (mu1/k1 dZ1/df(k2_) + mu2/k2 dZ1/df(k1_))

        where

                kij_ = ki_ + kj_

                    =>  kij^2 = |kij_|^2 = ki^2 + kj^2 + 2 ki kj nuij
                    =>  muij = kij_.s_ / kij = (ki mui + kj muj) / kij

        and

                k123_ = k1_ + k2_ + k3_

                    =>  k123^2 = |k1_ + k2_ + k3_|^2 = k1^2 + k2^2 + k3^2 + 2 k1 k2 nu12 + 2 k2 k3 nu23 + 2 k3 k1 nu31
                    =>  mu123 = k123_.s_ / k123 = (k1 mu1 + k2 mu2 + k3 mu3) / k123

    */

    /* Fiducials */
    fid_bias_t *bias = kernVar -> bias;
    fid_rsd_t *rsd = kernVar -> rsd;
    fid_btst_t *btst = kernVar -> btst;

    /* Compute all variables */
    bool computedWork = (kernVar -> computeWork) ? !_kernels_work_var(kernVar, 3) : false;

    /* First order kernel struct */
    kern_t *kernVarZ1 = ((kern_t**) kernVar -> kernWork)[0];

    /* Second order kernel struct */
    kern_t *kernVarZ2 = ((kern_t**) kernVar -> kernWork)[1];

    /* Third order kernel struct */
    kern_t *kernVarZ3 = ((kern_t**) kernVar -> kernWork)[2];

    /* Declare kernels */
    double h2Kernels[3][2];

    double gammaSin[3];

    double z1Kernels[3];
    double dz1Kernels[3];

    double z2Kernels[3];
    double dz2Kernels[3];


    /* Calculate kernels */

    /* F3(k1_, k2_, k3_) + G3(k1_, k2_, k3_) */
    double h3Kernels[2];
    kernels_btst_h3(kernVar, h3Kernels);

    /* Z1(k1_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[0];
    z1Kernels[0] = kernels_z1(kernVar);
    dz1Kernels[0] = kernels_dz1_f(kernVar);

    /* Z1(k2_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[1];
    z1Kernels[1] = kernels_z1(kernVar);
    dz1Kernels[1] = kernels_dz1_f(kernVar);

    /* Z1(k3_) */
    kernVarZ1 -> mu[0] = kernVarZ3 -> mu[2];
    z1Kernels[2] = kernels_z1(kernVar);
    dz1Kernels[2] = kernels_dz1_f(kernVar);


    size_t index = 0;

    for (size_t i = 0; i < kernVarZ2 -> kernOrder; i++)
      {
        for (size_t j = i; j < kernVarZ2 -> kernOrder; j++)
          {
            /*

                    F2(k1_, k2_) + G2(k1_, k2_) + γ(k1_, k2_)

                for k1, k2, nu12 given in the following order: (k[0],k[1],nu[0]) -> (k[0],k[2],nu[1]) -> (k[1],k[2],nu[2])

            */

            /* Variables */
            kernVarZ2 -> k[0] = kernVarZ3 -> k[i];
            kernVarZ2 -> mu[0] = kernVarZ3 -> mu[i];

            kernVarZ2 -> k[1] = kernVarZ3 -> k[j + 1];
            kernVarZ2 -> mu[1] = kernVarZ3 -> mu[j + 1];

            kernVarZ2 -> nu[0] = kernVarZ3 -> nu[index];

            /* Kernels */
            kernels_btst_h2(kernVar, h2Kernels[index]);

            gammaSin[index] = kernels_gamma(kernVarZ2);

            /* Z2'(k1_, k2_; k3_) */
            z2Kernels[index] = bias -> b1 * h2Kernels[index][0] + rsd -> f * (kernVarZ3 -> mu[index + 3] * kernVarZ3 -> mu[index + 3]) * h2Kernels[index][1]
                            + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1]
                                                                                       + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i] )
                            - (bias -> b1 * btst -> a2Ga - bias -> c2Ga) / 2. * gammaSin[index]
                            + bias -> b2 / 2.;

            dz2Kernels[index] = (kernVarZ3 -> mu[index + 3] * kernVarZ3 -> mu[index + 3]) * h2Kernels[index][1]
                                + kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * z1Kernels[j + 1]
                                                                                + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * z1Kernels[i] )
                                + rsd -> f * kernVarZ3 -> k[6] * kernVarZ3 -> mu[6] / 4. * ( kernVarZ2 -> mu[0] / kernVarZ2 -> k[0] * dz1Kernels[j + 1]
                                                                                           + kernVarZ2 -> mu[1] / kernVarZ2 -> k[1] * dz1Kernels[i] );

            index++;
          }
      }

    /* dZ3/df(k1_, k2_, k3_) */
    double dz3Kernel = (kernVarZ3 -> mu[6] * kernVarZ3 -> mu[6]) * h3Kernels[1]

                       + kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * z2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * z2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * z2Kernels[2] )
                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[2] / kernVarZ3 -> k[2] * dz2Kernels[0]
                            + kernVarZ3 -> mu[1] / kernVarZ3 -> k[1] * dz2Kernels[1]
                            + kernVarZ3 -> mu[0] / kernVarZ3 -> k[0] * dz2Kernels[2] )

                       + kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * h2Kernels[0][1] * z1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * h2Kernels[1][1] * z1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * z1Kernels[0] )
                       + rsd -> f * kernVarZ3 -> mu[6] * kernVarZ3 -> k[6] / 3.
                          * ( kernVarZ3 -> mu[3] / kernVarZ3 -> k[3] * h2Kernels[0][1] * dz1Kernels[2]
                            + kernVarZ3 -> mu[4] / kernVarZ3 -> k[4] * h2Kernels[1][1] * dz1Kernels[1]
                            + kernVarZ3 -> mu[5] / kernVarZ3 -> k[5] * h2Kernels[2][1] * dz1Kernels[0] );

    /* Reset computeWork */
    kernVar -> computeWork = computedWork;

    return dz3Kernel;
}




/*  ------------------------------------------------------------------------------------------------------  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
/*  ------------------------------------------------------------------------------------------------------  */

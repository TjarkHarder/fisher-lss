#include "../include/flss.h"


int main()
{

    /**  Initialise Global Parameters and Global Functions  **/

    flss_ini();


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */

    /**  Set the Fiducials  **/

    /* Change the file (in the '/fiducials/' directory): The files must have the same structure as the examples */
    fiducials_set_file('P', "pk.dat"); // Set the file for the linear power spectrum
    fiducials_set_file('F', "fiducials.dat"); // Set the file for constant / redshift dependent variables

    /* Read and setup the fiducial interpolation functions */
    fiducials_setup();

    /* Remove fiducials by setting their functions to zero */
    _fidCtrxC2_ = _zeroFunc_;
    _fidCtrxC4_ = _zeroFunc_;

    /* Remove sets of fiducials by setting flags to false */
//    _fidInclBias_ = false; // Bias
//    _fidInclRSD_ = false; // Redshift space distortions
//    _fidInclCtr_ = false; // Counter terms
//    _fidInclSn_ = false; // Shot noise

/*  ----------------------------------------------------  */

    /**  Number of Threads in Parallel Regions  **/

    flss_set_threads(10);

/*  ----------------------------------------------------  */

    /**  Set the Random Number Generator's Seed  **/

    flss_set_random_seed(42);

/*  ----------------------------------------------------  */

    /**  LOOP INTEGRATION  **/

    /* Power Spectrum */

      {
        /* One-Loop Integration Parameters */
        intgrt_t *intgrt1Loop = spec_info_get_loop_integrate(_idSpecPnl_, 0);

        /* Use Divonne Routine */
        integrate_set_routine(intgrt1Loop, _idIntgrtVegas_);

        /* Adjust VEGAS Parameters */
        intgrt_vegas_t *vegas1Loop = integrate_get_vegas(intgrt1Loop);

        integrate_vegas_set_wucalls(vegas1Loop, 1000); // Number of warm up calls
        integrate_vegas_set_mcalls(vegas1Loop, 10000); // Minimum number of main calls
        integrate_vegas_set_ucalls(vegas1Loop, 1.5); // Update number of main calls after failed iteration

        integrate_vegas_set_mloop(vegas1Loop, 10); // Maximum number of iterations

        integrate_vegas_set_chisqerr(vegas1Loop, 0.5); // Maximum Χ^2 err
        integrate_vegas_set_relerr(vegas1Loop, 1.e-3); // Precision goal

        integrate_vegas_set_verbose(vegas1Loop, 0); // 0 : No output, 1 : Information on convergence of integral

        /* Reduce Number of Loops */
        spec_info_set_loop_order(_idSpecPnl_, 1); // Replace 1 by 0 to only use the tree-level power spectrum
      }


    /* Power Spectrum Derivatives */

      {
        /* One-Loop Integration Parameters */
        intgrt_t *intgrt1Loop = spec_info_get_loop_integrate(_idSpecDPnl_, 0);

        /* Use Divonne Routine */
        integrate_set_routine(intgrt1Loop, _idIntgrtVegas_);

        /* Adjust VEGAS Parameters */
        intgrt_vegas_t *vegas1Loop = integrate_get_vegas(intgrt1Loop);

        integrate_vegas_set_wucalls(vegas1Loop, 10000);
        integrate_vegas_set_mcalls(vegas1Loop, 50000);
        integrate_vegas_set_ucalls(vegas1Loop, 1);

        integrate_vegas_set_mloop(vegas1Loop, 10);

        integrate_vegas_set_chisqerr(vegas1Loop, 0.5);
        integrate_vegas_set_relerr(vegas1Loop, 1.e-3);

        integrate_vegas_set_verbose(vegas1Loop, 0);

        /* Reduce Number of Loops */
        spec_info_set_loop_order(_idSpecDPnl_, 1); // Replace 1 by 0 to only use the tree-level power spectrum derivatives
      }



    /**  BIN-AVERAGE INTEGRATION  **/

    /* Average Flags */

    flss_set_avr_shape_flag(_idAvrLine_, false);
    flss_set_avr_shape_flag(_idAvrTri_, false);

    /* Line Bin-Volume */

      {
        /* Set the shape-bin-volume function */
        avr_set_shape_vol_func(_idAvrLineVol_, _idAvrLineVolFuncAna_);
      }


    /* Triangle Bin-Volume */

      {
        /* Set the shape-bin-volume function */
        avr_set_shape_vol_func(_idAvrTriVol_, _idAvrTriVolFuncNum_);

        /* Average Integration Parameters */
        intgrt_t *intgrt = avr_get_integrate(_idAvrTriVol_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrt, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonne = integrate_get_divonne(intgrt);

        integrate_divonne_set_epsabs(divonne, 0.);
        integrate_divonne_set_epsrel(divonne, 1.e-4);

        integrate_divonne_set_maxeval(divonne, 1000000);
      }


    /* Line Bin-Average Integration */

      {
        /* Average Integration Parameters */
        intgrt_t *intgrt = avr_get_integrate(_idAvrLine_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrt, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonne = integrate_get_divonne(intgrt);

        integrate_divonne_set_epsabs(divonne, 0.);
        integrate_divonne_set_epsrel(divonne, 1.e-4);

        integrate_divonne_set_maxeval(divonne, 1000000);
      }


    /* Triangle Bin-Average Integration */

      {
        /* Average Integration Parameters */
        intgrt_t *intgrt = avr_get_integrate(_idAvrTri_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrt, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonne = integrate_get_divonne(intgrt);

        integrate_divonne_set_epsabs(divonne, 0.);
        integrate_divonne_set_epsrel(divonne, 1.e-4);

        integrate_divonne_set_maxeval(divonne, 1000000);
      }


    /* PP Covariance Matrix Integration */

      {
        /* Average PP Gaussian Covariance Matrix Integration Parameters */
        intgrt_t *intgrtGauss = avr_get_integrate(_idAvrCovPPGauss_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrtGauss, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonneGauss = integrate_get_divonne(intgrtGauss);

        integrate_divonne_set_epsabs(divonneGauss, 0.);
        integrate_divonne_set_epsrel(divonneGauss, 1.e-3);

        integrate_divonne_set_maxeval(divonneGauss, 1000000);


        /* Average PP Non-Gaussian Covariance Matrix Integration Parameters */
        intgrt_t *intgrtNGauss = avr_get_integrate(_idAvrCovPPNGauss_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrtNGauss, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonneNGauss = integrate_get_divonne(intgrtNGauss);

        integrate_divonne_set_epsabs(divonneNGauss, 0.);
        integrate_divonne_set_epsrel(divonneNGauss, 1.e-3);

        integrate_divonne_set_key2(divonneNGauss, 47);

        integrate_divonne_set_maxpass(divonneNGauss, 50);

        integrate_divonne_set_maxchisq(divonneNGauss, 10.);
        integrate_divonne_set_mindeviation(divonneNGauss, 0.25);

        integrate_divonne_set_maxeval(divonneNGauss, 1000000);


        /* Average PP Non-Gaussian (Infinitesimal Limit) Covariance Matrix Integration Parameters */
        intgrt_t *intgrtNGaussInf = avr_get_integrate(_idAvrCovPPNGaussInf_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrtNGaussInf, _idIntgrtCQUAD_);

        /* Adjust QAGS Parameters */
        intgrt_cquad_t *cquadNGaussInf = integrate_get_cquad(intgrtNGaussInf);

        integrate_cquad_set_abserr(cquadNGaussInf, 0.);
        integrate_cquad_set_relerr(cquadNGaussInf, 1.e-3);
        integrate_cquad_set_limit(cquadNGaussInf, 100000);
      }


    /* BB Covariance Matrix Integration */

      {
        /* Average BB Gaussian Covariance Matrix Integration Parameters */
        intgrt_t *intgrtGauss = avr_get_integrate(_idAvrCovBBGauss_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrtGauss, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonneGauss = integrate_get_divonne(intgrtGauss);

        integrate_divonne_set_epsabs(divonneGauss, 0.);
        integrate_divonne_set_epsrel(divonneGauss, 1.e-2);

        integrate_divonne_set_key2(divonneGauss, 47);

        integrate_divonne_set_maxpass(divonneGauss, 100);

        integrate_divonne_set_maxchisq(divonneGauss, 10.);
        integrate_divonne_set_mindeviation(divonneGauss, 0.25);

        integrate_divonne_set_maxeval(divonneGauss, 10000000);


        /* Average BB Non-Gaussian Covariance Matrix Integration Parameters */
        intgrt_t *intgrtNGauss = avr_get_integrate(_idAvrCovBBNGauss_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrtNGauss, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonneNGauss = integrate_get_divonne(intgrtGauss);

        integrate_divonne_set_epsabs(divonneNGauss, 0.);
        integrate_divonne_set_epsrel(divonneNGauss, 1.e-2);

        integrate_divonne_set_key2(divonneNGauss, 47);

        integrate_divonne_set_maxpass(divonneNGauss, 100);

        integrate_divonne_set_maxchisq(divonneNGauss, 10.);
        integrate_divonne_set_mindeviation(divonneNGauss, 0.25);

        integrate_divonne_set_maxeval(divonneNGauss, 10000000);
      }


    /* PB Covariance Matrix Integration */

      {
        /* Average PB Non-Gaussian Covariance Matrix Integration Parameters */
        intgrt_t *intgrtNGauss = avr_get_integrate(_idAvrCovPBNGauss_);

        /* Change the Integration Routine */
        integrate_set_routine(intgrtNGauss, _idIntgrtDivonne_);

        /* Adjust Divonne Parameters */
        intgrt_divonne_t *divonneNGauss = integrate_get_divonne(intgrtNGauss);

        integrate_divonne_set_epsabs(divonneNGauss, 0.);
        integrate_divonne_set_epsrel(divonneNGauss, 1.e-2);

        integrate_divonne_set_key2(divonneNGauss, 47);

        integrate_divonne_set_maxpass(divonneNGauss, 100);

        integrate_divonne_set_maxchisq(divonneNGauss, 10.);
        integrate_divonne_set_mindeviation(divonneNGauss, 0.25);

        integrate_divonne_set_maxeval(divonneNGauss, 10000000);
      }


    /**  SHAPE SAMPLING  **/

    /* Line Sampling */

    if (false)
      {
        /* File name of line sample in '/output/sample/' directory */
        const char *fileName = "lines.dat";

        /* sample_arg_t struct for k sampling */
        sample_arg_t *sampleArgK = sample_arg_new();
        sample_arg_set_min(sampleArgK, 0.01);
        sample_arg_set_max(sampleArgK, 0.25);
        sample_arg_set_step(sampleArgK, 0.01);

        /* sample_arg_t struct for mu sampling */
        sample_arg_t *sampleArgMu = sample_arg_new();
        sample_arg_set_min(sampleArgMu, -0.95);
        sample_arg_set_max(sampleArgMu, 0.95);
        sample_arg_set_step(sampleArgMu, 0.1);

        /* sample_shape_t struct for sampling of lines */
        sample_shape_t *sampleShapeLine;

        if (true)
          {
            sampleShapeLine = sample_shape_line(sampleArgK, sampleArgMu, fileName); // Sample lines and output result to file (if fileName is NULL nothing is output)
          }

        else
          {
            sampleShapeLine = sample_shape_input(fileName, sampleArgK, sampleArgMu, NULL);  // Read sample from file (must have ngaussect sampleArgs)
          }

        /* Store sample of lines in global _sampleShapeLine_ variable (also accessible via flss_get_sample_shape(_idShapeLine_)) */
        flss_set_sample_shape(_idShapeLine_, sampleShapeLine);

        /* Free memory */
        sampleArgK = sample_arg_free(sampleArgK);
        sampleArgMu = sample_arg_free(sampleArgMu);
        sampleShapeLine = sample_shape_free(sampleShapeLine);
      }


    /* Triangle Sampling */

    if (false)
      {
        /* File name of triangle sample in '/output/sample/' directory */
        const char *fileName = "triangles.dat";

        /* sample_arg_t struct for k sampling */
        sample_arg_t *sampleArgK = sample_arg_new();
        sample_arg_set_min(sampleArgK, 0.01);
        sample_arg_set_max(sampleArgK, 0.1);
        sample_arg_set_step(sampleArgK, 0.01);

        /* sample_arg_t struct for mu sampling */
        sample_arg_t *sampleArgMu = sample_arg_new();
        sample_arg_set_min(sampleArgMu, -0.95);
        sample_arg_set_max(sampleArgMu, 0.95);
        sample_arg_set_step(sampleArgMu, 0.1);

        /* sample_shape_t struct for sampling of triangles */
        sample_shape_t *sampleShapeTri;

        if (true)
          {
            sampleShapeTri = sample_shape_tri(sampleArgK, sampleArgMu, fileName);  // Sample triangles and output result to file (if fileName is NULL nothing is output)
          }

        else
          {
            sampleShapeTri = sample_shape_input(fileName, sampleArgK, sampleArgMu, NULL); // Read sample from file (must have ngaussect sampleArgs)
          }

        /* Store sample of trangles in global _sampleShapeTri_ variable (also accessible via flss_get_sample_shape(_idShapeTri_)) */
        flss_set_sample_shape(_idShapeTri_, sampleShapeTri);

        /* Free memory */
        sampleArgK = sample_arg_free(sampleArgK);
        sampleArgMu = sample_arg_free(sampleArgMu);
        sampleShapeTri = sample_shape_free(sampleShapeTri);
      }


    /* Redshift Sampling */

    if (false)
      {
        /* sample_arg_t struct for z sampling */
        sample_arg_t *sampleArgZ = sample_arg_new();
        sample_arg_set_min(sampleArgZ, 0.7);
        sample_arg_set_max(sampleArgZ, 1.9);
        sample_arg_set_step(sampleArgZ, 0.2);

        /* Sample the redshift linearly for the chosen sampleArgs */
        sample_raw_t *sampleRawZ = sample_raw_lin(sampleArgZ);

        /* Store sample of redshifts in global _sampleRedshift_ variable (also accessible via flss_get_sample_redshift()) */
        flss_set_sample_redshift(sampleRawZ);

        /* Free memory */
        sampleArgZ = sample_arg_free(sampleArgZ);
        sampleRawZ = sample_raw_free(sampleRawZ);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /**  Interpolate the Power Spectrum from File  **/

    if (false)
      {
        /* Set the source file of the power spectrum in the 'output/spec/' directory */
        spec_in_pnl_set_file("pnl.dat");

        /* Setup the _specPnl_ function */
        spec_in_pnl_setup();
      }


    /**  Interpolate the Power Spectrum Derivatives from File (need this for fast )  **/

    if (false)
      {
        /* Set the source file of the power spectrum in the 'output/spec/' directory */
        spec_in_dpnl_set_file("dpnl.dat");

        /* Setup the _specDPnl_ functions */
        spec_in_dpnl_setup();

        /* Remove some derivatives by setting their functions to zero; If data for some derivatives was missing from the file, the functions are set to zero automatically */
//        _specDPnlA3GaB_ = _zeroFunc_;
//        _specDPnlD3GaB_ = _zeroFunc_;
      }


/*  ----------------------------------------------------  */


    /**  Read the PP Covariance Matrix  **/

    if (false)
      {
        /* Set the source file of the pp covariance matrix in the '/output/cov/' directory */
        cov_in_pp_set_file("covpp.mat");

        /* Set the label of the pp covariance matrix */
        cov_in_pp_set_label("'covpp-gaussian'");

        /* Setup the _covPP_ and _covPPMat_ functions; TODO (Note): _covPP_ only works when no matrix has been setup */
        cov_in_pp_setup(NULL);
      }


    /**  Read the BB Covariance Matrix  **/

    if (false)
      {
        /* Set the source file of the bb covariance matrix in the '/output/cov/' directory */
        cov_in_bb_set_file("covbb.mat");

        /* Set the label of the bb covariance matrix */
        cov_in_bb_set_label("'covbb-gaussian'");

        /* Setup the _covBB_ and _covBBMat_ functions; TODO (Note): _covBB_ only works when no matrix has been setup */
        cov_in_bb_setup(NULL);
      }


    /**  Read the PB Covariance Matrix  **/

    if (false)
      {
        /* Set the source file of the pb covariance matrix in the '/output/cov/' directory */
        cov_in_pb_set_file("covpb.mat");

        /* Set the label of the pb covariance matrix */
        cov_in_pb_set_label("'covpb-gaussian'");

        /* Setup the _covPB_ and _covPBMat_ functions; TODO (Note): _covPB_ only works when no matrix has been setup */
        cov_in_pb_setup(NULL, NULL);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /* Calculate the non-linear power spectrum directly */

    if (false)
      {
        /* spec_arg_t struct */
        spec_arg_t *specArg = spec_arg_new(_idSpecPnl_);

        /* kern_t struct */
        kern_t *kern = spec_arg_get_kern(specArg);

        /* Set the temporal variable */
        kernels_set_z(kern, 0.9);

        /* Set the spatial variables */
        kernels_set_k(kern, 0, 0.01); // 0'th vertex' length
        kernels_set_k(kern, 1, 0.01); // 1'st vertex' length

        kernels_set_mu(kern, 0, 0.05); // 0'th vertex' orientation (w.r.t. the line-of-sight)
        kernels_set_mu(kern, 1, -0.05); // 1'st vertex' orientation (w.r.t. the line-of-sight)

        kernels_set_nu(kern, 0, 1, -1.); // angle between 0'th vertex and 1'st vertex

        /* Set the bin-widths */
        spec_arg_set_dk(specArg, 0.01);
        spec_arg_set_dmu(specArg, 0.1);

        /* Set the bin-volume */
        double resultBinVol[3];
        avr_shape_line_vol(specArg, resultBinVol);

        spec_arg_set_volume(specArg, resultBinVol[0]);

        printf("bin-vol: %e %e\n", resultBinVol[0], resultBinVol[1]);
        printf("bin-vol: %e\n", avr_shape_vol_direct(avr_shape_line_vol, specArg));
        printf("bin-vol (analytic): %e\n\n", avr_shape_vol_direct(avr_shape_line_vol_ana, specArg));

        /* Calculate the power spectrum */
        double result[3];
        spec_pnl_full(specArg, result);

        printf("pnl: %e %e\n", result[0], result[1]);
        printf("pnl: %e\n", _specPnl_(specArg, NULL));

        /* Calculate the bin-averaged power spectrum (requires interpolation of power spectrum) */
//        double resultAvr[3];
//        avr_shape_line(_specPnl_, specArg, resultAvr);
//
//        printf("pnl-avr: %e %e\n", resultAvr[0], resultAvr[1]);
//        printf("pnl-avr: %e\n\n", avr_shape_direct(avr_shape_line, _specPnl_, specArg));

        /* Free memory */
        specArg = spec_arg_free(specArg);
      }


/*  ----------------------------------------------------  */


    /* Calculate the tree-level bispectrum directly */

    if (false)
      {
        /* spec_arg_t struct */
        spec_arg_t *specArg = spec_arg_new(_idSpecBtr_);

        /* kern_t struct */
        kern_t *kern = spec_arg_get_kern(specArg);

        /* Set the temporal variable */
        size_t redshiftIndex = 0; // Index of the redshift in the sample
        kernels_set_z(kern, sample_raw_get_value(_sampleRedshift_, redshiftIndex));

        /* Set the spatial variables */
        size_t shapeIndex = 0; // Index of the triangle in the sample
        shape_t *shape = sample_shape_get_shape(_sampleShapeTri_, shapeIndex);

        kernels_set_k(kern, 0, shape_get_vertex_length(shape, 0)); // 0'th vertex' length
        kernels_set_k(kern, 1, shape_get_vertex_length(shape, 1)); // 1'st vertex' length
        kernels_set_k(kern, 2, shape_get_vertex_length(shape, 2)); // 2'nd vertex' length

        kernels_set_mu(kern, 0, shape_get_vertex_orientation(shape, 0)); // 0'th vertex' orientation (w.r.t. line-of-sight)
        kernels_set_mu(kern, 1, shape_get_vertex_orientation(shape, 1)); // 1'st vertex' orientation (w.r.t. line-of-sight)
        kernels_set_mu(kern, 2, shape_get_vertex_orientation(shape, 2)); // 2'nd vertex' orientation (w.r.t. line-of-sight)

        kernels_set_nu(kern, 0, 1, shape_get_vertex_angle(shape, 0, 1)); // angle between 0'th vertex and 1'st vertex
        kernels_set_nu(kern, 0, 2, shape_get_vertex_angle(shape, 0, 2)); // angle between 0'th vertex and 2'nd vertex
        kernels_set_nu(kern, 1, 2, shape_get_vertex_angle(shape, 1, 2)); // angle between 1'st vertex and 2'nd vertex

        /* Set the bin-widths */
        spec_arg_set_dk(specArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeTri_))));
        spec_arg_set_dmu(specArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeTri_))));

        /* Set the bin-volume */
        spec_arg_set_volume(specArg, shape_get_volume(shape));

        /* Calculate the bispectrum (requires interpolation of power spectrum) */
        double result[3];
        spec_btr_full(specArg, result);

        printf("btr: %e %e\n", result[0], result[1]);
        printf("btr: %e\n", _specBtr_(specArg, NULL));

        /* Calculate the bin-averaged power spectrum (requires interpolation of power spectrum) */
        double resultAvr[3];
        avr_shape_tri(_specBtr_, specArg, resultAvr);

        printf("btr-avr: %e %e\n", resultAvr[0], resultAvr[1]);
        printf("btr-avr: %e\n\n", avr_shape_direct(avr_shape_tri, _specBtr_, specArg));

        /* Free memory */
        specArg = spec_arg_free(specArg);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /* Calculate the PP covariance matrix directly */

    if (false)
      {
        /* cov_arg_t struct */
        cov_arg_t *covArg = cov_arg_new(_idCovPP_);

        /* spec_arg_t struct */
        spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

        /* kern_t struct */
        kern_t *kern = spec_arg_get_kern(specArg);

        /* Set the temporal variable */
        size_t redshiftIndex = 0; // Index of the redshift in the sample
        kernels_set_z(kern, sample_raw_get_value(_sampleRedshift_, redshiftIndex));

        /* Set spatial variables */
        size_t shapeIndex1 = 0; // Index of first line
        shape_t *shape1 = sample_shape_get_shape(_sampleShapeLine_, shapeIndex1);
        cov_arg_set_shape(covArg, shape1, 0); // Set first line

        size_t shapeIndex2 = 0; // Index of second line
        shape_t *shape2 = sample_shape_get_shape(_sampleShapeLine_, shapeIndex2);
        cov_arg_set_shape(covArg, shape2, 1); // Set second line

        kernels_set_k(kern, 0, shape_get_vertex_length(shape1, 0)); // First line's 0'th vertex' length
        kernels_set_k(kern, 1, shape_get_vertex_length(shape1, 1)); // First line's 1'st vertex' length
        kernels_set_k(kern, 2, shape_get_vertex_length(shape2, 0)); // Second line's 0'th vertex' length
        kernels_set_k(kern, 3, shape_get_vertex_length(shape2, 1)); // Second line's 1'st vertex' length

        kernels_set_mu(kern, 0, shape_get_vertex_orientation(shape1, 0)); // First line's 0'th vertex' orientation
        kernels_set_mu(kern, 1, shape_get_vertex_orientation(shape1, 1)); // First line's 1'st vertex' orientation
        kernels_set_mu(kern, 2, shape_get_vertex_orientation(shape2, 0)); // Second line's 0'th vertex' orientation
        kernels_set_mu(kern, 3, shape_get_vertex_orientation(shape2, 1)); // Second line's 1'st vertex' orientation

        kernels_set_nu(kern, 0, 1, shape_get_vertex_angle(shape1, 0, 1)); // First line's angle between 0'th vertex and 1'st vertex
        kernels_set_nu(kern, 2, 3, shape_get_vertex_angle(shape2, 0, 1)); // Second line's angle between 0'th vertex and 1'st vertex

        /* Set bin widths */
        cov_arg_set_dk(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeLine_))), 0); // Set first line's dk
        cov_arg_set_dk(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeLine_))), 1); // Set second lines's dk

        cov_arg_set_dmu(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeLine_))), 0); // Set first line's dmu
        cov_arg_set_dmu(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeLine_))), 1); // Set second line's dmu


        /* PP Covariance Matrix in the Infinitesimal Limit (requires interpolation of power spectrum) */
        flss_set_avr_shape_flag(_idShapeLine_, false); // Disable bin-averaging

        double resultGauss[3], resultNGauss[3], resultFull[3];

        cov_pp_gauss(covArg, resultGauss);
        cov_pp_ngauss(covArg, resultNGauss);
        cov_pp_full(covArg, resultFull);

        printf("covpp-gaussian: %e %e\n", resultGauss[0], resultGauss[1]);
        printf("covpp-non-gaussian: %e %e\n", resultNGauss[0], resultNGauss[1]);
        printf("covpp-full: %e %e\n", resultFull[0], resultFull[1]);

        /* PP Covariance Matrix for Finite Bins (requires interpolation of power spectrum) */
        flss_set_avr_shape_flag(_idShapeLine_, true); // Enable bin-averaging

        double resultGaussAvr[3], resultNGaussAvr[3], resultFullAvr[3];

        cov_pp_gauss(covArg, resultGaussAvr);
        cov_pp_ngauss(covArg, resultNGaussAvr);
        cov_pp_full(covArg, resultFullAvr);

        printf("covpp-gaussian-avr: %e %e\n", resultGaussAvr[0], resultGaussAvr[1]);
        printf("covpp-non-gaussian-avr: %e %e\n", resultNGaussAvr[0], resultNGaussAvr[1]);
        printf("covpp-full-avr: %e %e\n", resultFullAvr[0], resultFullAvr[1]);

        /* Free memory */
        covArg = cov_arg_free(covArg);
      }


/*  ----------------------------------------------------  */


    /* Calculate the BB covariance matrix directly */

    if (false)
      {
        /* cov_arg_t struct */
        cov_arg_t *covArg = cov_arg_new(_idCovBB_);

        /* spec_arg_t struct */
        spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

        /* kern_t struct */
        kern_t *kern = spec_arg_get_kern(specArg);

        /* Set the temporal variable */
        size_t redshiftIndex = 0; // Index of the redshift in the sample
        kernels_set_z(kern, sample_raw_get_value(_sampleRedshift_, redshiftIndex));

        /* Set spatial variables */
        size_t shapeIndex1 = 0; // Index of first triangle
        shape_t *shape1 = sample_shape_get_shape(_sampleShapeTri_, shapeIndex1);
        cov_arg_set_shape(covArg, shape1, 0); // Set first triangle

        size_t shapeIndex2 = 0; // Index of second triangle
        shape_t *shape2 = sample_shape_get_shape(_sampleShapeTri_, shapeIndex2);
        cov_arg_set_shape(covArg, shape2, 1); // Set second triangle

        kernels_set_k(kern, 0, shape_get_vertex_length(shape1, 0)); // First triangle's 0'th vertex' length
        kernels_set_k(kern, 1, shape_get_vertex_length(shape1, 1)); // First triangle's 1'st vertex' length
        kernels_set_k(kern, 2, shape_get_vertex_length(shape1, 2)); // First triangle's 2'nd vertex' length
        kernels_set_k(kern, 3, shape_get_vertex_length(shape2, 0)); // Second triangle's 0'th vertex' length
        kernels_set_k(kern, 4, shape_get_vertex_length(shape2, 1)); // Second triangle's 1'st vertex' length
        kernels_set_k(kern, 5, shape_get_vertex_length(shape2, 2)); // Second triangle's 2'nd vertex' length

        kernels_set_mu(kern, 0, shape_get_vertex_orientation(shape1, 0)); // First triangle's 0'th vertex' orientation
        kernels_set_mu(kern, 1, shape_get_vertex_orientation(shape1, 1)); // First triangle's 1'st vertex' orientation
        kernels_set_mu(kern, 2, shape_get_vertex_orientation(shape1, 2)); // First triangle's 2'ns vertex' orientation
        kernels_set_mu(kern, 3, shape_get_vertex_orientation(shape2, 0)); // Second triangle's 0'th vertex' orientation
        kernels_set_mu(kern, 4, shape_get_vertex_orientation(shape2, 1)); // Second triangle's 1'st vertex' orientation
        kernels_set_mu(kern, 5, shape_get_vertex_orientation(shape2, 2)); // Second triangle's 2'nd vertex' orientation

        kernels_set_nu(kern, 0, 1, shape_get_vertex_angle(shape1, 0, 1)); // First triangle's angle between 0'th vertex and 1'st vertex
        kernels_set_nu(kern, 0, 2, shape_get_vertex_angle(shape1, 0, 2)); // First triangle's angle between 0'th vertex and 2'nd vertex
        kernels_set_nu(kern, 1, 2, shape_get_vertex_angle(shape1, 1, 2)); // First triangle's angle between 1'st vertex and 2'nd vertex
        kernels_set_nu(kern, 3, 4, shape_get_vertex_angle(shape2, 0, 1)); // Second triangle's angle between 0'th vertex and 1'st vertex
        kernels_set_nu(kern, 3, 5, shape_get_vertex_angle(shape2, 0, 2)); // Second triangle's angle between 0'th vertex and 2'nd vertex
        kernels_set_nu(kern, 4, 5, shape_get_vertex_angle(shape2, 1, 2)); // Second triangle's angle between 1'st vertex and 2'nd vertex

        /* Set bin widths */
        cov_arg_set_dk(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeTri_))), 0); // Set first triangle's dk
        cov_arg_set_dk(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeTri_))), 1); // Set second triangle's dk

        cov_arg_set_dmu(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeTri_))), 0); // Set first triangle's dmu
        cov_arg_set_dmu(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeTri_))), 1); // Set second triangle's dmu


        /* BB Covariance Matrix in the Infinitesimal Limit (requires interpolation of power spectrum) */
        flss_set_avr_shape_flag(_idShapeTri_, false); // Disable bin-averaging

        double resultGauss[3], resultNGauss[3], resultFull[3];

        cov_bb_gauss(covArg, resultGauss);
        cov_bb_ngauss(covArg, resultNGauss);
        cov_bb_full(covArg, resultFull);

        printf("covbb-gaussian: %e %e\n", resultGauss[0], resultGauss[1]);
        printf("covbb-non-gaussian: %e %e\n", resultNGauss[0], resultNGauss[1]);
        printf("covbb-full: %e %e\n", resultFull[0], resultFull[1]);

        /* BB Covariance Matrix for Finite Bins (requires interpolation of power spectrum) */
        flss_set_avr_shape_flag(_idShapeTri_, true); // Enable bin-averaging

        double resultGaussAvr[3], resultNGaussAvr[3], resultFullAvr[3];

        cov_bb_gauss(covArg, resultGaussAvr);
        cov_bb_ngauss(covArg, resultNGaussAvr);
        cov_bb_full(covArg, resultFullAvr);

        printf("covbb-gaussian-avr: %e %e\n", resultGaussAvr[0], resultGaussAvr[1]);
        printf("covbb-non-gaussian-avr: %e %e\n", resultNGaussAvr[0], resultNGaussAvr[1]);
        printf("covbb-full-avr: %e %e\n", resultFullAvr[0], resultFullAvr[1]);

        /* Free memory */
        covArg = cov_arg_free(covArg);
      }


/*  ----------------------------------------------------  */


    /* Calculate the PB covariance matrix directly / indirectly */

    if (false)
      {
        /* cov_arg_t struct */
        cov_arg_t *covArg = cov_arg_new(_idCovPB_);

        /* spec_arg_t struct */
        spec_arg_t *specArg = cov_arg_get_spec_arg(covArg);

        /* kern_t struct */
        kern_t *kern = spec_arg_get_kern(specArg);

        /* Set the temporal variable */
        size_t redshiftIndex = 0; // Index of the redshift in the sample
        kernels_set_z(kern, sample_raw_get_value(_sampleRedshift_, redshiftIndex));

        /* Set spatial variables */
        size_t shapeIndex1 = 0; // Index of line
        shape_t *shape1 = sample_shape_get_shape(_sampleShapeLine_, shapeIndex1);
        cov_arg_set_shape(covArg, shape1, 0); // Set line

        size_t shapeIndex2 = 0; // Index of triangle
        shape_t *shape2 = sample_shape_get_shape(_sampleShapeTri_, shapeIndex2);
        cov_arg_set_shape(covArg, shape2, 1); // Set triangle

        kernels_set_k(kern, 0, shape_get_vertex_length(shape1, 0)); // Line's 0'th vertex' length
        kernels_set_k(kern, 1, shape_get_vertex_length(shape1, 1)); // Line's 1'st vertex' length
        kernels_set_k(kern, 2, shape_get_vertex_length(shape2, 0)); // Triangle's 0'th vertex' length
        kernels_set_k(kern, 3, shape_get_vertex_length(shape2, 1)); // Triangle's 1'st vertex' length
        kernels_set_k(kern, 4, shape_get_vertex_length(shape2, 2)); // Triangle's 2'nd vertex' length

        kernels_set_mu(kern, 0, shape_get_vertex_orientation(shape1, 0)); // Line's 0'th vertex' orientation
        kernels_set_mu(kern, 1, shape_get_vertex_orientation(shape1, 1)); // Line's 1'st vertex' orientation
        kernels_set_mu(kern, 2, shape_get_vertex_orientation(shape2, 0)); // Triangle's 0'th vertex' orientation
        kernels_set_mu(kern, 3, shape_get_vertex_orientation(shape2, 1)); // Triangle's 1'st vertex' orientation
        kernels_set_mu(kern, 4, shape_get_vertex_orientation(shape2, 2)); // Triangle's 2'nd vertex' orientation

        kernels_set_nu(kern, 0, 1, shape_get_vertex_angle(shape1, 0, 1)); // Line's angle between 0'th vertex and 1'st vertex
        kernels_set_nu(kern, 2, 3, shape_get_vertex_angle(shape2, 0, 1)); // Triangle's angle between 0'th vertex and 1'st vertex
        kernels_set_nu(kern, 2, 4, shape_get_vertex_angle(shape2, 0, 2)); // Triangle's angle between 0'th vertex and 2'nd vertex
        kernels_set_nu(kern, 3, 4, shape_get_vertex_angle(shape2, 1, 2)); // Triangle's angle between 1'st vertex and 2'nd vertex

        /* Set bin widths */
        cov_arg_set_dk(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeLine_))), 0); // Set lines's dk
        cov_arg_set_dk(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_length(_sampleShapeTri_))), 1); // Set triangle's dk

        cov_arg_set_dmu(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeLine_))), 0); // Set line's dmu
        cov_arg_set_dmu(covArg, sample_arg_get_step(sample_raw_get_sample_arg(sample_shape_get_sample_raw_orientation(_sampleShapeTri_))), 1); // Set triangle's dmu


        /* PB Covariance Matrix in the Infinitesimal Limit (requires interpolation of power spectrum) */
        flss_set_avr_shape_flag(_idShapeLine_, false); // Disable bin-averaging
        flss_set_avr_shape_flag(_idShapeTri_, false); // Disable bin-averaging

        double resultGauss[3], resultNGauss[3], resultFull[3];

        cov_pb_gauss(covArg, resultGauss);
        cov_pb_ngauss(covArg, resultNGauss);
        cov_pb_full(covArg, resultFull);

        printf("covpb-gaussian: %e %e\n", resultGauss[0], resultGauss[1]);
        printf("covpb-non-gaussian: %e %e\n", resultNGauss[0], resultNGauss[1]);
        printf("covpb-full: %e %e\n", resultFull[0], resultFull[1]);

        /* PB Covariance Matrix for Finite Bins (requires interpolation of power spectrum) */
        flss_set_avr_shape_flag(_idShapeLine_, true); // Enable bin-averaging
        flss_set_avr_shape_flag(_idShapeTri_, true); // Enable bin-averaging

        double resultGaussAvr[3], resultNGaussAvr[3], resultFullAvr[3];

        cov_pb_gauss(covArg, resultGaussAvr);
        cov_pb_ngauss(covArg, resultNGaussAvr);
        cov_pb_full(covArg, resultFullAvr);

        printf("covpb-gaussian-avr: %e %e\n", resultGaussAvr[0], resultGaussAvr[1]);
        printf("covpb-non-gaussian-avr: %e %e\n", resultNGaussAvr[0], resultNGaussAvr[1]);
        printf("covpb-full-avr: %e %e\n", resultFullAvr[0], resultFullAvr[1]);

        /* Free memory */
        covArg = cov_arg_free(covArg);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /**  Calculate the Power Spectrum  **/

    if (false)
      {
        /* Get the Default Parameters  */
        spec_dat_t *pnlParams = spec_dat_pnl_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        spec_out_t *out = spec_dat_get_out(pnlParams);

        /* Remove labels */
//        spec_out_rm_label(out, "'pnl-tree'");
//        spec_out_rm_label(out, "'pnl-p22'");
//        spec_out_rm_label(out, "'pnl-p31'");
//        spec_out_rm_label(out, "'pnl-ctr'");
//        spec_out_rm_label(out, "'pnl-sn'");
//        spec_out_rm_label(out, "'pnl-full'");

        /* Set the output file in the '/output/spec/' directory; if NULL pointer is passed, nothing will be written to file */
        spec_out_set_file(out, "pnl.dat");
//        spec_out_set_file(out, NULL);

        /* Set the precision of the text file; Overrides latest 'spec_out_set_binary()' call and sets binary flag to false */
        spec_out_set_precision(out, 16);

        /* Set the flag determining if the file should be in text or in binary mode */
        spec_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = spec_dat_get_flags(pnlParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the Power Spectrum */

        dat_t *datPnl = spec_dat_poly(pnlParams);
        datPnl = dat_free(datPnl); // Free memory if not further needed


        /* Free memory */
        pnlParams = spec_dat_free(pnlParams);
      }


    /**  Calculate the Power Spectrum Derivatives  **/

    if (false)
      {
        /* Get the Default Parameters  */
        spec_dat_t *dpnlParams = spec_dat_dpnl_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        spec_out_t *out = spec_dat_get_out(dpnlParams);

        /* Remove labels */
//        spec_out_rm_label(out, "'dpnl-a2ga'");
//        spec_out_rm_label(out, "'dpnl-d2ga'");
//        spec_out_rm_label(out, "'dpnl-a3gaa'");
        spec_out_rm_label(out, "'dpnl-a3gab'"); // Does not enter at one-loop order in P
//        spec_out_rm_label(out, "'dpnl-d3gaa'");
        spec_out_rm_label(out, "'dpnl-d3gab'"); // Does not enter at one-loop order in P

//        spec_out_rm_label(out, "'dpnl-b1'");
//        spec_out_rm_label(out, "'dpnl-b2'");
//        spec_out_rm_label(out, "'dpnl-c2ga'");
        spec_out_rm_label(out, "'dpnl-bgam3'"); // Does not enter at one-loop order in P

//        spec_out_rm_label(out, "'dpnl-f'");
//        spec_out_rm_label(out, "'dpnl-sigv'");

//        spec_out_rm_label(out, "'dpnl-d'");
//        spec_out_rm_label(out, "'dpnl-h'");

//        spec_out_rm_label(out, "'dpnl-c0'");
        spec_out_rm_label(out, "'dpnl-c2'"); // Not included
        spec_out_rm_label(out, "'dpnl-c4'"); // Not included

//        spec_out_rm_label(out, "'dpnl-psn'");
        spec_out_rm_label(out, "'dpnl-bsn1'"); // Does not enter
        spec_out_rm_label(out, "'dpnl-bsn2'"); // Does not enter

//        spec_out_rm_label(out, "'dpnl-pk'");

        /* Set the logarithmic derivative flags of the spectra (i.e., multiply derivative with non-zero fiducial) (should not be set to true if also the respective logarithmic flags for the Fisher matrices are set to true : multiplied by fiducial twice...) */
//        spec_out_set_log(out, "a2ga", true);
//        spec_out_set_log(out, "d2ga", true);
//        spec_out_set_log(out, "a3gaa", true);
//        spec_out_set_log(out, "d3gaa", true);
//
//        spec_out_set_log(out, "b1", true);
//        spec_out_set_log(out, "b2", true);
//        spec_out_set_log(out, "c2ga", true);
//
//        spec_out_set_log(out, "f", true);
//        spec_out_set_log(out, "sigv", true);
//
//        spec_out_set_log(out, "c0", true);
//
//        spec_out_set_log(out, "pk", true);

        /* Set the output file in the '/output/spec/' directory; if NULL pointer is passed, nothing will be written to file */
        spec_out_set_file(out, "dpnl.dat");

        /* Set the precision of the text file; Overrides latest 'spec_out_set_binary()' call and sets binary flag to false */
        spec_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        spec_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = spec_dat_get_flags(dpnlParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the Power Spectrum Derivatives */

        dat_t *datDPnl = spec_dat_poly(dpnlParams);
        datDPnl = dat_free(datDPnl); // Free memory if not further needed


        /* Free memory */
        dpnlParams = spec_dat_free(dpnlParams);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /**  Calculate the PP Covariance Matrix  **/

    if (false)
      {
        /* Get the Default Parameters  */
        cov_mat_t *covppParams = cov_mat_pp_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        cov_out_t *out = cov_mat_get_out(covppParams);

        /* Remove labels */
//        cov_out_rm_label(out, "'covpp-gaussian'");
        cov_out_rm_label(out, "'covpp-non-gaussian'");
//        cov_out_rm_label(out, "'covpp-full'");

        /* Set the output file in the '/output/cov/' directory; if NULL pointer is passed, nothing will be written to file */
        cov_out_set_file(out, "covpp.mat");

        /* Set the precision of the text file; Overrides latest 'cov_out_set_binary()' call and sets binary flag to false */
        cov_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        cov_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = cov_mat_get_flags(covppParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the PP Covariance Matrix */

        mats_t *matsCovPP = cov_mat_poly(covppParams);
        matsCovPP = mats_free_full(matsCovPP); // Free memory if not further needed


        /* Free memory */
        covppParams = cov_mat_free(covppParams);
      }


    /**  Calculate the BB Covariance Matrix  **/

    if (false)
      {
        /* Get the Default Parameters  */
        cov_mat_t *covbbParams = cov_mat_bb_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        cov_out_t *out = cov_mat_get_out(covbbParams);

        /* Remove labels */
//        cov_out_rm_label(out, "'covbb-gaussian'");
        cov_out_rm_label(out, "'covbb-non-gaussian'");
//        cov_out_rm_label(out, "'covbb-full'");

        /* Set the output file in the '/output/cov/' directory; if NULL pointer is passed, nothing will be written to file */
        cov_out_set_file(out, "covbb.mat");

        /* Set the precision of the text file; Overrides latest 'cov_out_set_binary()' call and sets binary flag to false */
        cov_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        cov_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = cov_mat_get_flags(covbbParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the BB Covariance Matrix */

        mats_t *matsCovBB = cov_mat_poly(covbbParams);
        matsCovBB = mats_free_full(matsCovBB); // Free memory if not further needed


        /* Free memory */
        covbbParams = cov_mat_free(covbbParams);
      }


    /**  Calculate the PB Covariance Matrix  **/

    if (false)
      {
        /* Get the Default Parameters  */
        cov_mat_t *covpbParams = cov_mat_pb_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        cov_out_t *out = cov_mat_get_out(covpbParams);

        /* Remove labels */
//        cov_out_rm_label(out, "'covpb-gaussian'");
        cov_out_rm_label(out, "'covpb-non-gaussian'");
//        cov_out_rm_label(out, "'covpb-full'");

        /* Set the output file in the '/output/cov/' directory; if NULL pointer is passed, nothing will be written to file */
        cov_out_set_file(out, "covpb.mat");

        /* Set the precision of the text file; Overrides latest 'cov_out_set_binary()' call and sets binary flag to false */
        cov_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        cov_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = cov_mat_get_flags(covpbParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the PB Covariance Matrix */

        mats_t *matsCovPB = cov_mat_poly(covpbParams);
        matsCovPB = mats_free_full(matsCovPB); // Free memory if not further needed


        /* Free memory */
        covpbParams = cov_mat_free(covpbParams);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /**  Calculate the PP Fisher Matrix  **/

    if (false)
      {
        /* Get the Default Parameters  */
        fish_mat_t *fishppParams = fish_mat_pp_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        fish_out_t *out = fish_mat_get_out(fishppParams);

        /* Remove labels */
//        fish_out_rm_label(out, "a2ga");
//        fish_out_rm_label(out, "d2ga");
//        fish_out_rm_label(out, "a3gaa");
        fish_out_rm_label(out, "a3gab"); // Does not enter at one-loop order
//        fish_out_rm_label(out, "d3gaa");
        fish_out_rm_label(out, "d3gab"); // Does not enter at one-loop order

//        fish_out_rm_label(out, "b1");
//        fish_out_rm_label(out, "b2");
//        fish_out_rm_label(out, "c2ga");
        fish_out_rm_label(out, "bgam3"); // Does not enter at one-loop order

//        fish_out_rm_label(out, "f");
//        fish_out_rm_label(out, "sigv");

//        fish_out_rm_label(out, "d");
//        fish_out_rm_label(out, "h");

//        fish_out_rm_label(out, "c0");
        fish_out_rm_label(out, "c2"); // Not included
        fish_out_rm_label(out, "c4"); // Not included

//        fish_out_rm_label(out, "psn");
        fish_out_rm_label(out, "bsn1"); // Does not enter
        fish_out_rm_label(out, "bsn2"); // Does not enter

//        fish_out_rm_label(out, "pk");

        /* Set the logarithmic derivative flags of the spectra (i.e., multiply Fisher matrix with non-zero fiducial) */
        fish_out_set_log(out, "a2ga", true);
        fish_out_set_log(out, "d2ga", true);
        fish_out_set_log(out, "a3gaa", true);
        fish_out_set_log(out, "d3gaa", true);

        fish_out_set_log(out, "b1", true);
        fish_out_set_log(out, "b2", true);
        fish_out_set_log(out, "c2ga", true);

        fish_out_set_log(out, "f", true);
        fish_out_set_log(out, "sigv", true);

        fish_out_set_log(out, "c0", true);

        fish_out_set_log(out, "pk", true);

        /* Set the output file in the '/output/fish/' directory; if NULL pointer is passed, nothing will be written to file */
        fish_out_set_file(out, "fishpp.mat");

        /* Set the precision of the text file; Overrides latest 'fish_out_set_binary()' call and sets binary flag to false */
        fish_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        fish_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = fish_mat_get_flags(fishppParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the PP Fisher Matrix */

        mats_t *matsFishPP = fish_mat_poly(fishppParams);
        matsFishPP = mats_free_full(matsFishPP); // Free memory if not further needed


        /* Free memory */
        fishppParams = fish_mat_free(fishppParams);
      }


    /**  Calculate the BB Fisher Matrix  **/

    if (false)
      {
        /* Get the Default Parameters  */
        fish_mat_t *fishbbParams = fish_mat_bb_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        fish_out_t *out = fish_mat_get_out(fishbbParams);

        /* Remove labels */
//        fish_out_rm_label(out, "a2ga");
//        fish_out_rm_label(out, "d2ga");
//        fish_out_rm_label(out, "a3gaa");
        fish_out_rm_label(out, "a3gab"); // Does not enter at tree-level in B
//        fish_out_rm_label(out, "d3gaa");
        fish_out_rm_label(out, "d3gab"); // Does not enter at tree-level in B

//        fish_out_rm_label(out, "b1");
//        fish_out_rm_label(out, "b2");
//        fish_out_rm_label(out, "c2ga");
        fish_out_rm_label(out, "bgam3"); // Does not enter at tree-level in B

//        fish_out_rm_label(out, "f");
//        fish_out_rm_label(out, "sigv");

//        fish_out_rm_label(out, "d");
//        fish_out_rm_label(out, "h");

//        fish_out_rm_label(out, "c0");
        fish_out_rm_label(out, "c2"); // Not included
        fish_out_rm_label(out, "c4"); // Not included

        fish_out_rm_label(out, "psn"); // Does not enter in B
//        fish_out_rm_label(out, "bsn1");
//        fish_out_rm_label(out, "bsn2");

//        fish_out_rm_label(out, "pk");

        /* Set the logarithmic derivative flags of the spectra (i.e., multiply Fisher matrix with non-zero fiducial) */
        fish_out_set_log(out, "a2ga", true);
        fish_out_set_log(out, "d2ga", true);
        fish_out_set_log(out, "a3gaa", true);
        fish_out_set_log(out, "d3gaa", true);

        fish_out_set_log(out, "b1", true);
        fish_out_set_log(out, "b2", true);
        fish_out_set_log(out, "c2ga", true);

        fish_out_set_log(out, "f", true);
        fish_out_set_log(out, "sigv", true);

        fish_out_set_log(out, "c0", true);

        fish_out_set_log(out, "pk", true);

        /* Set the output file in the '/output/fish/' directory; if NULL pointer is passed, nothing will be written to file */
        fish_out_set_file(out, "fishbb.mat");

        /* Set the precision of the text file; Overrides latest 'fish_out_set_binary()' call and sets binary flag to false */
        fish_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        fish_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = fish_mat_get_flags(fishbbParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the BB Fisher Matrix */

        mats_t *matsFishBB = fish_mat_poly(fishbbParams);
        matsFishBB = mats_free_full(matsFishBB); // Free memory if not further needed


        /* Free memory */
        fishbbParams = fish_mat_free(fishbbParams);
      }


    /**  Calculate the PB Fisher Matrix  **/

    if (false)
      {
        /* Get the Default Parameters  */
        fish_mat_t *fishpbParams = fish_mat_pb_default();


        /* Adjust the Output Parameters */

        /* Get the output struct */
        fish_out_t *out = fish_mat_get_out(fishpbParams);

        /* Remove labels */
//        fish_out_rm_label(out, "a2ga");
//        fish_out_rm_label(out, "d2ga");
//        fish_out_rm_label(out, "a3gaa");
        fish_out_rm_label(out, "a3gab"); // Does not enter at one-loop order in P and at tree-level in B
//        fish_out_rm_label(out, "d3gaa");
        fish_out_rm_label(out, "d3gab"); // Does not enter at one-loop order in P and at tree-level in B

//        fish_out_rm_label(out, "b1");
//        fish_out_rm_label(out, "b2");
//        fish_out_rm_label(out, "c2ga");
        fish_out_rm_label(out, "bgam3"); // Does not enter at one-loop order in P and at tree-level in B

//        fish_out_rm_label(out, "f");
//        fish_out_rm_label(out, "sigv");

//        fish_out_rm_label(out, "d");
//        fish_out_rm_label(out, "h");

//        fish_out_rm_label(out, "c0");
        fish_out_rm_label(out, "c2"); // Not included
        fish_out_rm_label(out, "c4"); // Not included

//        fish_out_rm_label(out, "psn");
//        fish_out_rm_label(out, "bsn1");
//        fish_out_rm_label(out, "bsn2");

//        fish_out_rm_label(out, "pk");

        /* Set the logarithmic derivative flags of the spectra (i.e., multiply Fisher matrix with non-zero fiducial) */
        fish_out_set_log(out, "a2ga", true);
        fish_out_set_log(out, "d2ga", true);
        fish_out_set_log(out, "a3gaa", true);
        fish_out_set_log(out, "d3gaa", true);

        fish_out_set_log(out, "b1", true);
        fish_out_set_log(out, "b2", true);
        fish_out_set_log(out, "c2ga", true);

        fish_out_set_log(out, "f", true);
        fish_out_set_log(out, "sigv", true);

        fish_out_set_log(out, "c0", true);

        fish_out_set_log(out, "pk", true);

        /* Set the output file in the '/output/fish/' directory; if NULL pointer is passed, nothing will be written to file */
        fish_out_set_file(out, "fishpb.mat");

        /* Set the precision of the text file; Overrides latest 'fish_out_set_binary()' call and sets binary flag to false */
        fish_out_set_precision(out, 16);

        /* Set if the file should be in text or in binary */
        fish_out_set_binary(out, true);


        /* Adjust the Print Flags */

        /* Get the print flag struct */
        print_flags_t *flags = fish_mat_get_flags(fishpbParams);

        /* Set the flags to true or false */
        print_flags_set_main(flags, true);


        /* Calculate the PB Fisher Matrix */

        mats_t *matsFishPB = fish_mat_poly(fishpbParams); // Get the Fisher matrix via a UDU^T decomposition of the P+B covariance matrix; faulty for gaussian + non-gaussian P+B covariance matrix
//        mats_t *matsFishPB = fish_mat_poly_eig(fishpbParams); // Get the Fisher matrix via an eigendecomposition of the P+B covariance matrix; discards negative eigenvalues

        matsFishPB = mats_free_full(matsFishPB); // Free memory if not further needed


        /* Free memory */
        fishpbParams = fish_mat_free(fishpbParams);
      }


/*  ----------------------------------------------------  */
/*  ----------------------------------------------------  */


    /**  Free Memory  **/

    flss_free();


    return 0;
}

# rOMT_spdup
This repository contains the modified version of the rOMT algorithm part in rOMT repository. <br />
We speed up the original algorithm and include the version of handling 2D data. <br />
Go to Inverse -> GNblock_u.m for editing history<br />

For more info about the theory and details about the project, please go to https://github.com/xinan-nancy-chen/rOMT.

# Sample cases for demonstration
Run ``driver_CAA.m`` which contains a sample data case with default paramters. This is a healthy rat brain DCE-MRI data. It takes about 2 hours and 45 minutes to run locally with 2.6 GHz Intel Core i7 and 16G memory on MacOS, while the original version before improvement (https://github.com/xinan-nancy-chen/rOMT) took about 37 hours on a CPU cluster with 40 cores.<br />

The results of Lagrangian Speed Map (without QuickBundle) and Velocity Flux Vectors will pop up automatically, both ran and visualized in Matlab_R2019b.<br />


## (1) We compare this new version with the previous version

The previous version takes the output of current loop as the input of the next loop, and here we call it the *dependent* version. However, we have developed a method to run each loop *independently*, thus makes parallelization realistic. <br />

Next we show an example of a healthy rat brain data tagged as 'C294'. <br />

For the independent (new) version (cfg.reinitR = 1), <br />

<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Speed/C294_LagSpeed_E31_53.png" width="470" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagPathlines_E31_53.png" width="470" /> 
</p>

<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagFluxVector_E31_53.png" width="470" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagAdvDiffVector_E31_53.png" width="470" /> 
</p>

Top-left: Lagrangian speed map; Top-right: Lagrangian pathlines; Bottom-left: Velocity flux vectors; Bottom-right: Classification of flux vectors into advection and diffusion.<br />

For the dependent (previous) version (cfg.reinitR = 0), <br />
<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Speed/C294_LagSpeed_E31_53.png" width="470" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagPathlines_E31_53.png" width="470" /> 
</p>

<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagFluxVector_E31_53.png" width="470" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagAdvDiffVector_E31_53.png" width="470" /> 
</p>

Even though the dependent version gives smoother pathlines and higher intensity of activities, but by comparing the results with the input brain images (12 frames shown as follows), the independent version actually fits the data best in terms of the spatial distribution of activities. <br />

<p float="left">
  <img src="test_results/C294/C294_InputData_E31_53_t_1.png" width="190" />
  <img src="test_results/C294/C294_InputData_E31_53_t_2.png" width="190" /> 
  <img src="test_results/C294/C294_InputData_E31_53_t_3.png" width="190" />
  <img src="test_results/C294/C294_InputData_E31_53_t_4.png" width="190" /> 
  <img src="test_results/C294/C294_InputData_E31_53_t_5.png" width="190" />
  <img src="test_results/C294/C294_InputData_E31_53_t_6.png" width="190" /> 
  <img src="test_results/C294/C294_InputData_E31_53_t_7.png" width="190" />
  <img src="test_results/C294/C294_InputData_E31_53_t_8.png" width="190" /> 
  <img src="test_results/C294/C294_InputData_E31_53_t_9.png" width="190" />
  <img src="test_results/C294/C294_InputData_E31_53_t_10.png" width="190" /> 
  <img src="test_results/C294/C294_InputData_E31_53_t_11.png" width="190" />
  <img src="test_results/C294/C294_InputData_E31_53_t_12.png" width="190" /> 
</p>

This comparison illustrates the viability of upgrading the dependent into an indenpendent algorithm. <br />

## We add the smoothing of velocity fields in post-processing trying to reduce the noise

Having identified that the indenpendent version gives noisier and less smooth results which is within expectation by introducing new noise in each loop, we add a  smoothing step to the velocity fields during post-processing. <br />

Since our algorithm is dynamic, the smoothing can be divided into two categories: <br />

smoothing the velocity field in the <br />
(1) time space (whos intensity is controlled by paramter Svt) <br />
(2) spatial space (controlled by paramter Svs) <br />

According to testing, tuning on Svs is way more sensitive than on Svt. <br />

Here we present two example rat brain cases, one healthy 'C294' and one with CAA (Cerebral Amyloid Angiopathy) 'C371' to show the effect of smoothing. <br />

### The healthy case. <br />
<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Speed/C294_LagSpeed_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Speed/C294_LagSpeed_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Speed/C294_LagSpeed_E31_53.png" width="300" /> <br />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagPathlines_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Vectors/C294_LagPathlines_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Vectors/C294_LagPathlines_E31_53.png" width="300" /> <br />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagFluxVector_E31_53.png" width="300" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Vectors/C294_LagFluxVector_E31_53.png" width="300" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Vectors/C294_LagFluxVector_E31_53.png" width="300" /> <br />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagAdvDiffVector_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Vectors/C294_LagAdvDiffVector_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Vectors/C294_LagAdvDiffVector_E31_53.png" width="300" /> <br />
</p>

The first column is witout any smoothing on velocities. The middle column is smoothed at Svt = 5, Svs = 1; The third column is smoothed at Svt = 10, Svs = 5.<br />

### The diseased case. <br />
<p float="left">
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061321/Speed/C371_LagSpeed_E31_53.png" width="300" /> 
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Speed/C371_LagSpeed_E31_53.png" width="300" /> 
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Speed/C371_LagSpeed_E31_53.png" width="300" /> <br />
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061321/Vectors/C371_LagPathlines_E31_53.png" width="300" /> 
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Vectors/C371_LagPathlines_E31_53.png" width="300" /> 
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Vectors/C371_LagPathlines_E31_53.png" width="300" /> <br />
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061321/Vectors/C371_LagFluxVector_E31_53.png" width="300" />
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Vectors/C371_LagFluxVector_E31_53.png" width="300" />
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Vectors/C371_LagFluxVector_E31_53.png" width="300" /> <br />
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061321/Vectors/C371_LagAdvDiffVector_E31_53.png" width="300" /> 
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set002_061421/Vectors/C371_LagAdvDiffVector_E31_53.png" width="300" /> 
  <img src="test_results/C371/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_061421/Vectors/C371_LagAdvDiffVector_E31_53.png" width="300" /> <br />
</p>

The first column is witout any smoothing on velocities. The middle column is smoothed at Svt = 5, Svs = 1; The third column is smoothed at Svt = 10, Svs = 5.<br />


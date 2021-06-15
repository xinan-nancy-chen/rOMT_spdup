# rOMT_spdup
This repository contains the modified version of the rOMT algorithm part in rOMT repository. <br />
We speed up the original algorithm and include the version of handling 2D data. <br />
Go to Inverse -> GNblock_u.m for editing history<br />

For more info about the theory and details about the project, please go to https://github.com/xinan-nancy-chen/rOMT

# Sample cases for demonstration
Run ``driver_CAA.m`` which contains a sample data case with default paramters. This is a healthy rat brain DCE-MRI data. It takes about 2 hours and 45 minutes to run on a local computer with i7 and 16G memory while the original version before improvement took about 24 hours.<br />

The results of Lagrangian Speed Map (without QuickBundle) and Velocity Flux Vectors will pop up automatically, both ran and visualized in Matlab_R2019b.<br />


## (1) We compare this new version with the previous version

The previous version takes the output of current loop as the input of the next loop, and here we call it the *dependent* version. However, we have developed a method to run each loop *independently*, thus makes parallelization realistic. <br />

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

Even though the dependent version gives smoother pathlines and higher intensity of activities, but by comparing the results with the input brain images, the independent version actually fits the data best. <br />


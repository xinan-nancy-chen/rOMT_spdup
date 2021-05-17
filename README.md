# rOMT_spdup
This repository contains the modified version of the rOMT algorithm part in rOMT repository. <br />
We speed up the original algorithm and include the version of handling 2D data. <br />
Go to Inverse -> GNblock_u.m for editing history<br />

# Sample case for demonstration
Run ``driver_CAA.m`` which contains a sample data case with default paramters. This is a healthy rat brain DCE-MRI data. It takes about 2 hours and 45 minutes to run on a local computer with i7 and 16G memory while the original version before improvement took about 24 hours.<br />

The results of Lagrangian Speed Map (without QuickBundle) and Velocity Flux Vectors will pop up automatically, both ran and visualized in Matlab_R2019b.<br />

For the independent version, 
![LagSpeed1](⁨test_results⁩/C294⁩/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/⁨LPPA_set001_051721⁩/C294_LagSpeed_E31_53.png)
![LagVector1](C294_LagFluxVector_E31_52.png)

![LagSpeed](C294_LagSpeed_E31_52.png)
![LagVector](C294_LagFluxVector_E31_52.png)

For more info about the theory and details about the project, please go to https://github.com/xinan-nancy-chen/rOMT

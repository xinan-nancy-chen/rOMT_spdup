# Introduction of rOMT_spdup

The regularized optimal mass transport (rOMT) problem can be described as follows. Given the initial mass distribution function $\rho_0(x)\geqslant0$ and the final one $\rho_1(x)\geqslant0$ defined on a bounded region $\Omega\subseteq\mathbb{R}^3$, one solves

$$\underset{\rho,v}{\text{min}}\quad \int_0^T\int_{\Omega}\left\lVert v(t,x)\right\rVert^2\rho(t,x)dx dt $$
subject to

$$\frac{\partial\rho}{\partial t} + \nabla\cdot(\rho v) = \sigma\Delta\rho, $$

$$\rho(0,x) = \rho_0(x), \quad\rho(T,x) = \rho_1(x)$$

where a temporal dimension $t\in[0,T]$ is added to the transport process. In the above expression, $\rho(t,x)$ is the dynamic density function; $v(t,x)$ is the velocity field defining the flow from $\rho_0$ to $\rho_1$; constant $\sigma>0$ is the diffusion coefficient.

This repository contains the modified version of algorithm in rOMT repository (https://github.com/xinan-nancy-chen/rOMT). We <br />
(1) Upgraded the previous code and realized 91% percent runtime reduction; <br />
(2) Offered the option to run multiple rOMT loops in parallel to further cut down on the runtime; <br />
(3) Offered the option to smoothen pathlines for the parallelized version in post-processing; <br />
(4) Included the version of handling <em>2D</em> data in rOMT code, in addition to <em>3D</em> data.

Go to Inverse -> GNblock_u.m for editing history<br />

For more info about the theory and details about the project, please go to https://github.com/xinan-nancy-chen/rOMT or

> -- <cite>[Visualizing fluid flows via regularized optimal mass transport with applications to neuroscience][1]</cite>,

[1]: https://arxiv.org/abs/2201.07307

Contact Xinan Chen at chenx7@mskcc.org for questions.

# Sample cases for demonstration

## (A) Gaussian Spheres
Run ``driver_gauss.m`` which contains a synthetic geometric dataset with default paramters. It took about 26 minutes on a 2.6 GHz Intel Core i7-9750H, 16G RAM, running macOS Mojave (version 10.14.6) with MATLAB 2019b.<br />

The inputs are 5 successive <em>3D</em> Gaussian Spheres (shown as follows, colormap='jet'). <br />

<p float="left">
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/gauss2_data_E01.png" width="190" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/gauss2_data_E02.png" width="190" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/gauss2_data_E03.png" width="190" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/gauss2_data_E04.png" width="190" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/gauss2_data_E05.png" width="190" />
</p>

The Lagrangian results: <em>Speed Map</em> (without QuickBundle), <em>Pe Map</em> , <em>Pathlines</em>, <em>Speed-lines</em>, <em>Péclet-lines</em> and <em>Velocity Flux Vectors</em>, will pop up automatically.<br />
<p float="left">
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/Pathlines/gauss2_LagPathlines_E01_05.png" width="300" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/Pathlines/gauss2_LagSpdlines_E01_05.png" width="300" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/Pathlines/gauss2_LagPelines_E01_05.png" width="300" />
</p>
<em>Lagrangian Results. Left: Pathlines; Middle: Speed-lines; Right: Péclet-lines</em>.<br /><br />
<p float="left">
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/Pathlines/gauss2_LagFluxVector_E01_05.png" width="300" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/Speed/gauss2_LagSpeed_E01_05.png" width="300" />
  <img src="test_results/gauss2/diff_2e3_tj_1_dt_0.4_nt_10_ti_1_tf_4_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth0_rreinit1_source0_pcg60/LPPA_set001_103021/Speed/gauss2_LagPe_E01_05.png" width="300" />
</p>
<em>Lagrangian Results. Left: Velocity flux vectors; Middle: Speed map; Right: Pe map</em>.<br />

## (B) Rat Brain MRI
Run ``driver_CAA.m`` which contains a sample data case with default paramters. This is DCE-MRI data from a healthy rat brain. It took about 2 hours and 45 minutes to run unparalleled locally with 2.6 GHz Intel Core i7 and 16G RAM on MacOS, while the original version before improvement took about 37 hours on a CPU cluster with 40 cores. If run in paralled, it took about 24 minutes on the cluster with the same configuration. <br />

The inputs are 12 successive <em>3D</em> images within a masked region (shown as follows). <br />

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

The Lagrangian results: <em>Speed Map</em> (without QuickBundle), <em>Pathlines</em> and <em>Velocity Flux Vectors</em>, will pop up automatically, all ran and visualized with Matlab_R2019b.<br />

Note that if run unparalleled, we put the final interpolated image from the previous loop into the next loop as the initial image. If run in parallel, we use the original input images as initial images in each loop. In ``driver_CAA.m``, by setting cfg.reinitR = 0, it will give the unparallel version, and 1 for the parallel version. The latter will give 10-fold faster results, which may however result in unsmooth pathlines.

### (1) We compare unparallel and parallel results

Next we show an example of a healthy rat brain data tagged as 'C294', comparing the Lagrangian results of unparalleled and parallel code. <br />

### For the unparallel version (cfg.reinitR = 0), <br />
<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Speed/C294_LagSpeed_E31_53.png" width="300" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagPathlines_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit0_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagFluxVector_E31_53.png" width="300" />
</p>
<em>Lagrangian Results. Left: speed map; Middle: Lagrangian pathlines; Right: Velocity flux vectors</em>.<br />

<br />

### For the parallel version (cfg.reinitR = 1), <br />

<p float="left">
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Speed/C294_LagSpeed_E31_53.png" width="300" />
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagPathlines_E31_53.png" width="300" /> 
  <img src="test_results/C294/diff_2e3_tj_2_dt_0.4_nt_10_ti_31_tf_51_uini_0_beta_0.0001_R_gamma_0.008_dtri1_rsmooth1_rreinit1_source0_dilate3_pcg60/LPPA_set001_051721/Vectors/C294_LagFluxVector_E31_53.png" width="300" />
</p>



As the above figures exhbit, both algorithms give the approximately same distributions of speed maps and the overall same directions of flows, while the unparallel version gives smoother pathlines that penetrate deeper.  This comparison illustrates the viability of upgrading the unparallel into parallel algorithm if short runtime is emphasized. <br />

### We add the smoothening of velocity fields in post-processing

Having identified that the parallel version gives noisier and less smooth results which is within expectation by constantly introducing new noise in each loop, we therefore offer an option of adding a smoothing step to the velocity fields during post-processing. <br />

Since our algorithm is dynamic, the smoothing can be divided into two categories: <br />

smoothening the velocity field in the <br />
(1) time space (whose intensity is controlled by paramter Svt) <br />
(2) spatial space (controlled by paramter Svs) <br />

According to testing, tuning on Svs is way more sensitive than on Svt. Here we present two example rat brain cases, one healthy 'C294' and one with CAA (Cerebral Amyloid Angiopathy) 'C371' to show the effect of smoothing. <br />

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
</p>

The first column is without any smoothing on velocities. The middle column is smoothed at Svt = 5, Svs = 1; The third column is smoothed at Svt = 10, Svs = 5.<br />

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
</p>

The first column is witout any smoothing on velocities. The middle column is smoothed at Svt = 5, Svs = 1; The third column is smoothed at Svt = 10, Svs = 5.<br />


Assessing the drivers of the syphilis epidemic among men-who-have-sex-with-men

Suraj Balakrishna et al.

Programs

1.	STM9_f.R – MAIN PROGRAM: to model the syphilis epidemic among MSM in Switzerland
2.	STM8_c.R – Counterfactual scenario: Impact of no change in transmission exposure (nsP) during the observation period. 
3.	STM8_d.R – Counterfactual scenario: Impact of syphilis screening frequency in MSM with HIV diagnosis and with occasional partners 
4.	STM8_e.R – Counterfactual scenario: Impact of syphilis screening frequency in MSM without HIV diagnosis and with occasional partners 
5.	STM3_c_test.R – Sensitivity analysis: ‘nsCAI’ as transmission risk instead of ‘nsP’
6.	STM7_h.R – Sensitivity analysis: Assume same transmission rate for MSM with and without HIV diagnosis
7.	STM7_e.R – Sensitivity analysis: Assume same transmission rate for MSM with and without HIV diagnosis but taking into account a potential underreporting of syphilis cases among MSM without HIV diagnosis
8.	STM7_g1.R – Sensitivity analysis: Impact of a possible overestimation of transmission risk among MSM without HIV diagnosis assuming overestimation to be 25%
9.	STM7_g2.R – Sensitivity analysis: Impact of a possible overestimation of transmission risk among MSM without HIV diagnosis assuming overestimation to be 50%
10.	STM7_g3.R – Sensitivity analysis: Impact of a possible overestimation of transmission risk among MSM without HIV diagnosis assuming overestimation to be 75%
11.	STM7_d.R – Sensitivity analysis: Infectiousness during the latent stage to be 1% of that in primary and secondary stage of syphilis
12.	STM7_d2.R – Sensitivity analysis: Infectiousness during the latent stage to be 10% of that in primary and secondary stage of syphilis 
13.	STM8_sens_mod_d.R – Counterfactual scenario of the sensitivity analysis: Impact of syphilis screening frequency in MSM with HIV diagnosis and with occasional partners assuming that infectiousness during the latent stage to be 1% of that in primary and secondary stage of syphilis
14.	STM8_sens_mod_d2.R – Counterfactual scenario of the sensitivity analysis: Impact of syphilis screening frequency in MSM with HIV diagnosis and with occasional partners assuming that infectiousness during the latent stage to be 10% of that in primary and secondary stage of syphilis
15.	STM9_b.R – Sensitivity analysis: Model fitting through maximization of the likelihood by assuming Poisson distributed incident cases of syphilis.
16.	STM9_loglik_modfit.R – Plotting the model fitting by minimizing SSR and maximizing loglikelihood
17.	STM9_gelman.R – Model diagnostics: Gelman and Rubin's convergence diagnostic

ODE Models

1.	STM9_b_model_birth.cpp – ODE model written in C++ compatible with STM9_f.R, STM3_c_test.R, STM7_h.R, STM7_g1.R, STM7_g2.R, STM7_g3.R, and STM9_b.R 
2.	STM7_e.cpp – ODE model written in C++ compatible with STM7_e.R
3.	STM7_d_cpp – ODE model written in C++ compatible with STM7_d.R and STM7_d2.R


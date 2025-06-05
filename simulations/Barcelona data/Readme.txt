Code used to run the simulations described in the paper:

Unbiased estimates using temporally aggregated outcome data in time series analysis: generalization to different outcomes, exposures and types of aggregation

First, one needs to load the functions to fit the models by sourcing the following file:

- agr_functions_v36_cpp.r

The following files need to be located in the same folder:

- agr_logLik_cpp.cpp
- d1func_aggr_cpp.cpp
- d1score_ind_cpp.cpp
- d2func_aggr_cpp.cpp

Second, one needs to run the following 4 files, which simulate responses and fit the
D|D, W|D, M|D and Dow|D models to the simulated data. Running these has a high computational cost (days).

- For mortality and temperature: 01_bcn_sim_data_mort_temp_cpp_v14.r
- For hospitalizations and temperature: 02_bcn_sim_data_hosp_temp_cpp_v14_i.r
- For mortality and air pollution: 03_bcn_sim_mort_ap_0101_dow_holi_i.r
- For hospitalizations and air pollution: 04_bcn_sim_hosp_ap_0101_1pct_dow_holi_i.r

These files use the following data files, which contain the exposure data and the parameters used to simulate the outcomes:

- bcn_ts_exposures.RData
- parameters_mort_temp.RData
- parameters_hosp_temp.RData
- parameters_mort_ap_01_dow_holi.RData
- parameters_hosp_ap_01_dow_holi.RData

After running the files, results will be saved into different RData files, with names: 

- results_sim_bcn_v14_mort_temp_1.RData, results_sim_bcn_v14_mort_temp_2.RData, ...
- results_sim_bcn_v14_hosp_temp_1.RData, results_sim_bcn_v14_hosp_temp_1.RData, ...
- results_sim_bcn_v14_mort_ap_dow_holi_0101_1.RData, results_sim_bcn_v14_mort_ap_dow_holi_0101_2.RData, ...
- results_sim_bcn_v14_hosp_ap_dow_holi_0101_1.RData, results_sim_bcn_v14_hosp_ap_dow_holi_0101_2.RData, ...

For the simulations based on N/10, we conduct the simulations using a similar process.
Running files:

- For mortality and temperature: 05_bcn_sim_data_mort10_temp_cpp_v14_i.r
- For hospitalizations and temperature: 06_bcn_sim_data_hosp10_temp_cpp_v14_i.r
- For mortality and air pollution: 07_bcn_sim_mort10_ap_0101_dow_holi_i.r
- For hospitalizations and air pollution: 08_bcn_sim_hosp10_ap_0101_1pct_dow_holi_i.r

will create the results files:

- results_sim_bcn10_v14_mort_temp_1.RData, results_sim_bcn10_v14_mort_temp_2.RData, ...
- results_sim_bcn10_v14_hosp_temp_1.RData, results_sim_bcn10_v14_hosp_temp_2.RData, ...
- results_sim_bcn10_v14_mort_ap_dow_holi_0101_1.RData, results_sim_bcn10_v14_mort_ap_dow_holi_0101_2.RData, ... 
- results_sim_bcn10_v14_hosp_ap_dow_holi_0101_1pct_1.RData, results_sim_bcn10_v14_hosp_ap_dow_holi_0101_1pct_2.RData, ...

Third, one needs to run the following file, which takes the results of the models across simulations and calculates the properties: bias, rmse, coverage, power. 
It also formats the results in a convenient way to do the plots:

- 09_process_results_0101_dow_holi_v14_4.r

Finally, one needs to run the following file, which produces the plots:

- 10_plots_bcn_v14_0101_dow_holi_5.r


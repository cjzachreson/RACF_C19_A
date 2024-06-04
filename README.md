# RACF_C19_A
code for modelling COVID-19 outbreaks in RACF environments with an ABM. 

## Overview
The parent folder, `full_model`, contains two subdirectories: 
`population_generator`, and `outbreak_simulator`. The subdirectory named `outbreak_simulator` contains all necessary input data and scripts to reproduce the figures of the paper entitled "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities" by Zachreson et al., however the population generator is included as well for reproducibility of methods. 

### outbreak_simulator
The outbreak simulator performs stochastic infectious disease transmission simulations within the populations specified by the input files contained in the subdirectory `outbreak_simulator/input`.

The simulations described in the paper listed above can be performed by calling the main functions contained in the files: 

1. `Main_RACF_OB_v9.jl` (main test results)
2. `Main_RACF_OB_v10_LD_unmitigated.jl` (supplemental results for scenarios with only lockdown policies)
3. `Main_RACF_R0_Hom_v9.jl` (homogeneous R0 calibration)
4. `Main_RACF_R0_v9.jl` (structured population R0 calibration)
5. `Main_RACF_FS_v9.jl` (structured population final size statistics)

Initialisation and calls to these functions can be performed by modifying the file `full_run.jl`, which contains the `include()` commands needed to initialise the model's submodules, which should be called in order: 

1. `header_RACF.jl`
2. `Setup_RACF_v9.jl`
3. `Networks_RACF_v9.jl`
4. `Diseases_RACF_v9.jl`
5. `Agents_RACF_v9.jl`
6. `Facility_Structure_v9.jl`
7. `Outbreak_Response_RACF_v9.jl`
8. `Transmission_Dynamics_v9.jl`

The `main` file of interest can then be included, and the corresponding `main_` function can be called (note that these main functions have different names depending on the main file included for the corresponding set of simulations)

Raw simulation output is recorded in the designated output folder: 

1. `Main_RACF_OB_v9.jl` $\rightarrow$  output_v9_OB_test
2. `ain_RACF_OB_v10_LD_unmitigated.jl` $\rightarrow$   output_OB_LD_only
3. `Main_RACF_R0_Hom_v9.jl` $\rightarrow$  output_v9_R0_hom_test_3
4. `Main_RACF_R0_v9.jl` $\rightarrow$  output_v9_R0_test_3
5. `Main_RACF_FS_v9.jl` $\rightarrow$  output_v9_FS_test_3

Analysis of summary statistics can be performed by running the MATLAB script `compute_sumary_stats.m`, `compute_summary_stats_LD_only.m`, or `FS_dist_vs_R0.m`, which record output in the designated subdirectory: `outbreak_simulator/analysis/output_summary_stats`

### population_generator
The subdirectory `population_generator` contains input data and scripts used to generate the text files that are read in to the ABM to initialise the population. This is performed in two steps, which correspond to subdirectories `facility_generator_step_1` and `prepare_agent_files_step_2`. 

#### facility_generator_step_1
The subdirectory `facility_generator_step_1` contains two main scripts, one for generating homogeneous populations (used for computing basic reproductive ratios and testing final size output for different transmission rates) and one for generating structured populations to represent contact patterns in Residential Aged Care Facilities (RACFs). Running the scripts produces output files from the input data contained in the `.csv` files named `homogeneous_facility_characteristics.csv` and `hypothetical_facility_characteristics.csv`, for which the latter parameterises the generator for RACF facility populations and the former parameterises the generator for homogeneous populations. After running the scripts, the output folder will be updated to contain `.csv` files with the population details required for initialising the ABM. 

#### prepare_agent_files_step_2

After generating the agent files, they must be formatted for read-in to the ABM outbreak simulator, which is coded in Julia language. The scripts `Prepare_agent_files_homogeneous_facility_output_v1.m` and `Prepare_agent_files_hypothetical_facility_output_v5.m` perform this function for the homogeneous and RACF agent files, respectively. The output will be storred in the designated subdirectory e.g. `output\agents_model_ready\...` After generating these output files, they should be copied into the designated input folder contained in  the `outbreak_simulator` subdirectory. 





#scripts for setup of utilities for RACF model:

#TODO: open and compose config file, could be implemented as a string. 

macro Name(arg)
    string(arg)
end

# TODO: config should list facility id and outbreak onset date for population generator 
config_str = "parameter variable name, human-readable description, value \n" # initialises an empty string 

# add control parameter values to config:
# response delay  
    description = "delay between outbreak declaration and implementation of infection control (days)"
    config_line = (@Name(delay_val_i) * ", " * description * ", " * "$delay_val_i" )
    config_str = config_str*config_line*"\n"

# compliance with scheduled testing 
    description = "staff compliance with scheduled testing (probability of compliance with each test)"
    config_line = (@Name(test_compliance_j) * ", " * description * ", " * "$test_compliance_j" )
    config_str = config_str*config_line*"\n"


# setup random number generator: 
global seed_infections = 13 #seed 13 outbreak as of 2022_08_29
    description = "random number seed for infection dynamics"
    config_line = (@Name(seed_infections) * ", " * description * ", " * "$seed_infections" )
    config_str = config_str*config_line*"\n"


global rng_infections = MersenneTwister(seed_infections)


global seed_testing = 1 #seed 1 testing as of 2022_08_29
    description = "random number seed for testing"
    config_line = (@Name(seed_testing) * ", " * description * ", " * "$seed_testing" )
    config_str = config_str*config_line*"\n"


global rng_testing = MersenneTwister(seed_testing)

#NOTE: another seed for network rewiring and contact selection. 
global seed_contacts = 1
    description = "random number seed for contact dynamics"
    config_line = (@Name(seed_contacts) * ", " * description * ", " * "$seed_contacts" )
    config_str = config_str*config_line*"\n"


global rng_contacts = MersenneTwister(seed_contacts)


global seed_immunity = 1
    description = "random number seed for immunity status (neut titers), [only used if IMMUNITY_FROM_DIST = true]"
    config_line = (@Name(seed_immunity) * ", " * description * ", " * "$seed_immunity" )
    config_str = config_str*config_line*"\n"

global rng_immunity = MersenneTwister(seed_immunity)


function apply_seed_offset(seed_offset::Int64)

    global seed_contacts += seed_offset
    global seed_testing += seed_offset 
    global seed_infections += seed_offset 
    global seed_immunity += seed_offset 

    global rng_contacts = MersenneTwister(seed_contacts)
    global rng_testing = MersenneTwister(seed_testing)
    global rng_infections = MersenneTwister(seed_infections)
    global rng_immunity = MersenneTwister(seed_immunity)


end

#define discrete time step
global dt = 0.1
    description = "discrete time step"
    config_line = (@Name(dt) * ", " * description * ", " * "$dt" )
    config_str = config_str*config_line*"\n"

config_str = config_str*"\n"*"*** contact rates *** \n"

# contact rate per resident per day 
# 2022 09 27 : reducing contact rate 
# 2022 09 28 : contact rate now stratified by needs of residents (aligning assumptions from contact weighting to contact rates)
global contact_rate_per_resident_per_day = 8.0#10.0
    description = "average number of contacts sampled per resident on each day (high needs, reg. needs is 1/3 of this)"
    config_line = (@Name(contact_rate_per_resident_per_day) * ", " * description * ", " * "$contact_rate_per_resident_per_day" )
    config_str = config_str*config_line*"\n"

# rate of resient background contacts per day (i.e., in shared communal spaces)
global bkg_contact_rate_per_resident_per_day = 5.0
    description = "average number of random resident-resident contacts on each day"
    config_line = (@Name(bkg_contact_rate_per_resident_per_day) * ", " * description * ", " * "$bkg_contact_rate_per_resident_per_day" )
    config_str = config_str*config_line*"\n"

config_str = config_str*"\n"*"***relative contact durations*** \n"
# some global control parameters: 
# relative contact strength for high-needs residents etc: 
global w_worker_worker = 1.0
    description = "between workers"
    config_line = (@Name(w_worker_worker) * ", " * description * ", " * "$w_worker_worker" )
    config_str = config_str*config_line*"\n"

global w_reg_needs = 2.0
    description = "between workers and residents with regular need levels"
    config_line = (@Name(w_reg_needs) * ", " * description * ", " * "$w_reg_needs" )
    config_str = config_str*config_line*"\n"

global w_high_needs = 6.0 #2.0 # factor by which to increase contact strength between workers and high-needs residents  
    description = "between workers and residents with high needs levels"
    config_line = (@Name(w_high_needs) * ", " * description * ", " * "$w_high_needs" )
    config_str = config_str*config_line*"\n"

global w_res_same_room = 10.0
    description = "between residents sharing the same room"
    config_line = (@Name(w_res_same_room) * ", " * description * ", " * "$w_res_same_room" )
    config_str = config_str*config_line*"\n"

global w_res_diff_room = 1.0 
    description = "between residents in different rooms"
    config_line = (@Name(w_res_diff_room) * ", " * description * ", " * "$w_res_diff_room" )
    config_str = config_str*config_line*"\n"

# NOTE: as of 2022 08 30, w_res_diff_room is not used because interactions between
# residents in different rooms is handled by random sampling at a fixed rate. 


config_str = config_str*"\n"*"***\n"


#flag for toggling on and off vaccine-derived protection: 
#NOTE: 2022 09 25 : moving immunity flags from setup to preamble. 
#global IMMUNITY = false
description = "boolean flag for whether to implement vaccine-derived protection"
config_line = (@Name(IMMUNITY) * ", " * description * ", " * "$IMMUNITY" )
config_str = config_str*config_line*"\n"

#global IMMUNITY_FROM_DIST = true
description = "boolean flag for whether to draw immunity status from file or generate it for each run"
config_line = (@Name(IMMUNITY_FROM_DIST) * ", " * description * ", " * "$IMMUNITY_FROM_DIST" )
config_str = config_str*config_line*"\n"



# if IMMUNITY_FROM_DIST is true, define the cohort-specific distributions 
# in this implementation these are the normal distributions from which
# ln(neut) values are drawn: 
# NOTE: this is not robust code, in future versions the distributions 
# can be setup in a way that's more flexible. 
if IMMUNITY_FROM_DIST
    # these parameters were estimated for the 
    # exemplar facilities studied in this specific 
    # set of experiments. These are not universal parameters. 
    # see the population generator optional step 5
    # immunity analysis. 
    
    global log_mu_res = -1.47
        description = "mean of natural log neuts (residents)"
        config_line = (@Name(log_mu_res) * ", " * description * ", " * "$log_mu_res" )
        config_str = config_str*config_line*"\n"
    global log_sig_res = 1.34
        description = "standard deviation of natural log neuts (residents)"
        config_line = (@Name(log_sig_res) * ", " * description * ", " * "$log_sig_res" )
        config_str = config_str*config_line*"\n"

    global log_mu_wG = -1.79
        description = "mean of natural log neuts (general staff)"
        config_line = (@Name(log_mu_wG) * ", " * description * ", " * "$log_mu_wG" )
        config_str = config_str*config_line*"\n"
    global log_sig_wG = 1.43
        description = "standard deviation of natural log neuts (general staff)"
        config_line = (@Name(log_sig_wG) * ", " * description * ", " * "$log_sig_wG" )
        config_str = config_str*config_line*"\n"

    global log_mu_wM = -1.64
        description = "mean of natural log neuts (medical staff)"
        config_line = (@Name(log_mu_wM) * ", " * description * ", " * "$log_mu_wM" )
        config_str = config_str*config_line*"\n"
    global log_sig_wM = 1.51
        description = "standard deviation of natural log neuts (medical staff)"
        config_line = (@Name(log_sig_wM) * ", " * description * ", " * "$log_sig_wM" )
        config_str = config_str*config_line*"\n"
end


# flag for toggling infection control during active outbreaks (false means unmitigated)
#global OUTBREAK_CONTROL = true moving this to MAIN for convenience. 

#note the parameters below are not used currently so I'm not adding them to the config. 
# is an active outbreak declared? 
global ACTIVE_OUTBREAK = false #initialised as false 

# global parameter, probability of testing an agent each day 
global p_test_per_day = 0.2 # only used if testing is not sub-divided by agent type 

# NOTE: implementing scheduled tests for workers, 2022 09 12
# see Outbreak_response.jl for test scheduling. 
# each of the scheduled tests will be subject to the probabilities 
# noted below (note that residents have these probabilities applied 
# every day, under the baseline or outbreak conditions ) 

config_str = config_str*"\n"*"***testing probabilities*** \n"

#NOTE: this is now assigned from control parameter for test compliance. 
global p_test_per_day_workers_baseline = test_compliance_j#1.0
    description = "compliance probability for worker tests (no outbreak)"
    config_line = (@Name(p_test_per_day_workers_baseline) * ", " * description * ", " * "$p_test_per_day_workers_baseline" )
    config_str = config_str*config_line*"\n"

global p_test_per_day_residents_baseline = 0.05
    description = "robability of a resident testing if no outbreak"
    config_line = (@Name(p_test_per_day_residents_baseline) * ", " * description * ", " * "$p_test_per_day_residents_baseline" )
    config_str = config_str*config_line*"\n"

#NOTE: this is now assigned from control parameter for test compliance. 
global p_test_per_day_workers_outbreak = test_compliance_j#1.0
    description = "compliance probability for worker tests (outbreak)"
    config_line = (@Name(p_test_per_day_workers_outbreak) * ", " * description * ", " * "$p_test_per_day_workers_outbreak" )
    config_str = config_str*config_line*"\n"

global p_test_per_day_residents_outbreak = 0.5
    description = "probability of a resident testing during an outbreak"
    config_line = (@Name(p_test_per_day_residents_outbreak) * ", " * description * ", " * "$p_test_per_day_residents_outbreak" )
    config_str = config_str*config_line*"\n"

global p_test_if_symptomatic = 1.0
    description = "probability of a symptomatic individual testing if on site"
    config_line = (@Name(p_test_if_symptomatic) * ", " * description * ", " * "$p_test_if_symptomatic" )
    config_str = config_str*config_line*"\n"


config_str = config_str*"\n"*"*** pairwise efficacy of infection control ***  \n"

#efficacy of infection control for pairwise contacts: 
# triggered by active outbreaks 
global eff_IC_worker_worker = 0.5
    description = "between workers"
    config_line = (@Name(eff_IC_worker_worker) * ", " * description * ", " * "$eff_IC_worker_worker" )
    config_str = config_str*config_line*"\n"

global eff_IC_worker_resident = 0.9
    description = "between workers and residents"
    config_line = (@Name(eff_IC_worker_resident) * ", " * description * ", " * "$eff_IC_worker_resident" )
    config_str = config_str*config_line*"\n"

global eff_IC_resident_resident = 0.2
    description = "between residents and residents"
    config_line = (@Name(eff_IC_resident_resident) * ", " * description * ", " * "$eff_IC_resident_resident" )
    config_str = config_str*config_line*"\n"


config_str = config_str*"\n"*"*** outbreak response *** \n"

# factor by which background interaction rates are reduced during active outbreak 
global RESIDENT_LOCKDOWN = false  
    description = "flag for whether contact rates between non-isolated residents will be reduced"
    config_line = (@Name(RESIDENT_LOCKDOWN) * ", " * description * ", " * "$RESIDENT_LOCKDOWN" )
    config_str = config_str*config_line*"\n"

global resident_lockdown_efficacy = 0.5 # representing discretionary changes 
    description = "fraction by which background contact rates are reduced for residents during outbreaks"
    config_line = (@Name(resident_lockdown_efficacy) * ", " * description * ", " * "$resident_lockdown_efficacy" )
    config_str = config_str*config_line*"\n"

# duration for which a worker is removed from the facility if they test positive 
global WORKER_CASE_ISOLATION = false
    description = "flag whether or not to furlough staff who test positive"
    config_line = (@Name(WORKER_CASE_ISOLATION) * ", " * description * ", " * "$WORKER_CASE_ISOLATION" )
    config_str = config_str*config_line*"\n"

global RESIDENT_CASE_ISOLATION = false # 2022 09 12
    description = "flag whether or not to isolate residents who test positive"
    config_line = (@Name(RESIDENT_CASE_ISOLATION) * ", " * description * ", " * "$RESIDENT_CASE_ISOLATION" )
    config_str = config_str*config_line*"\n"

global resident_isolation_efficacy = 0.0
    description = "reduction in background contact rate for isolated residents"
    config_line = (@Name(resident_isolation_efficacy) * ", " * description * ", " * "$resident_isolation_efficacy" )
    config_str = config_str*config_line*"\n"

global removal_period = 0.0 # changed from 14 to 7 on 2022 09 12
    description = "period of time (days) for isolation or furlough"
    config_line = (@Name(removal_period) * ", " * description * ", " * "$removal_period" )
    config_str = config_str*config_line*"\n"

#2022 09 21 setup variables moved from Main: 
# flag for toggling infection control during active outbreaks (false means unmitigated)
global OUTBREAK_CONTROL = false
    description = "flag whether to apply any outbreak control measures"
    config_line = (@Name(OUTBREAK_CONTROL) * ", " * description * ", " * "$OUTBREAK_CONTROL" )
    config_str = config_str*config_line*"\n"


config_str = config_str*"\n"*"***\n"
#flag telling the system to write the transmission tree and detections to file: 
global write_transmission_tree_flag = false
    description = "flag telling the system to write the transmission tree to file"
    config_line = (@Name(write_transmission_tree_flag) * ", " * description * ", " * "$write_transmission_tree_flag" )
    config_str = config_str*config_line*"\n"

config_str = config_str*"\n"*"*** type of index case *** \n"
# global control parameter: worker index case vs. resident index case (if both are true it can be either)
global resident_index_case = true
    description = "worker index case vs. resident index case (if both are true it can be either)"
    config_line = (@Name(resident_index_case) * ", " * description * ", " * "$resident_index_case" )
    config_str = config_str*config_line*"\n"

global worker_index_case = true
    description = "worker index case vs. resident index case (if both are true it can be either)"
    config_line = (@Name(worker_index_case) * ", " * description * ", " * "$worker_index_case" )
    config_str = config_str*config_line*"\n"



# paths input files (TODO: define functions for modifying these from list in csv):

config_str = config_str*"\n"*"***population data location***\n"
config_line = data_dirname 
config_str = config_str*config_line*"\n"
#data_dirname = pwd() * "\\network_constructor_output\\OB_2"


# write config str to output file. (note - applies global variable output_dir_OB)
config_file = open("$(output_dir_OB)\\config.txt", "w")
write(config_file, config_str)
close(config_file)


E_list_fname = data_dirname * "\\Edge_lists_test.csv"
residents_fname = data_dirname * "\\residents_for_Julia_test.csv"
workers_G_fname = data_dirname * "\\workers_G_for_Julia_test.csv"
workers_M_fname = data_dirname * "\\workers_M_for_Julia_test.csv"
N_list_fname = data_dirname * "\\Neighbour_lists_test.csv"


# conversion of input files to dataframe objects for read-in: 
E_list_str = DataFrame(CSV.File(E_list_fname, delim = "\t"))
N_lists_str = DataFrame(CSV.File(N_list_fname, delim = "\t"))
residents_str = DataFrame(CSV.File(residents_fname, delim = "\t"))
workers_G_str = DataFrame(CSV.File(workers_G_fname, delim = "\t"))
workers_M_str = DataFrame(CSV.File(workers_M_fname, delim = "\t"))

# set up the timeline mapping numeric weeks to dates 
# this is used when interpreting numeric week values from
# immunity status inputs. 


#abstract type definitions 
abstract type Agent_T end
abstract type Agents_T end 
abstract type Contact_T end 
abstract type Edge_T end 
abstract type E_list_T end 
abstract type E_list_temporal_T end 
abstract type N_list_T end 
abstract type N_list_temporal_multigraph_T end 
abstract type Room_T end 
abstract type Rooms_T end 
abstract type Disease_T end
abstract type Infection_T end
abstract type Immunity_Profile_T end
abstract type Vaccine_T end
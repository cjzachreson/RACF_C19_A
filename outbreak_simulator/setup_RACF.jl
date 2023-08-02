
#scripts for setup of utilities for RACF model:

# setup random number generator: 
global seed_infections = 13 #seed 13 outbreak as of 2022_08_29
global rng_infections = MersenneTwister(seed_infections)

global seed_testing = 1 #seed 13 outbreak as of 2022_08_29
global rng_testing = MersenneTwister(seed_testing)

#define discrete time step
global dt = 0.1


# some global control parameters: 
# relative contact strength for high-needs residents etc: 
global w_worker_worker = 1.0
global w_reg_needs = 2.0
global w_high_needs = 6.0 #2.0 # factor by which to increase contact strength between workers and high-needs residents  
global w_res_same_room = 10.0
global w_res_diff_room = 1.0 
# NOTE: as of 2022 08 30, w_res_diff_room is not used because interactions between
# residents in different rooms is handled by random sampling at a fixed rate. 

#flag for toggling on and off vaccine-derived protection: 
global IMMUNITY = false 

# flag for toggling infection control during active outbreaks (false means unmitigated)
global OUTBREAK_CONTROL = true
# is an active outbreak declared? 
global ACTIVE_OUTBREAK = false #initialised as false 

    # global parameter, probability of testing an agent each day 
    global p_test_per_day = 0.2 # only used if testing is not sub-divided by agent type 

    global p_test_per_day_workers_baseline = 0.2
    global p_test_per_day_residents_baseline = 0.05

    global p_test_per_day_workers_outbreak = 0.8
    global p_test_per_day_residents_outbreak = 0.5

    global p_test_if_symptomatic = 1.0


    #efficacy of infection control for pairwise contacts: 
    # triggered by active outbreaks 
    global eff_IC_worker_worker = 0.5
    global eff_IC_worker_resident = 0.9
    global eff_IC_resident_resident = 0.2

# factor by which background interaction rates are reduced 
global RESIDENT_LOCKDOWN = true 
    global resident_lockdown_efficacy = 0.9

# duration for which a worker is removed from the facility if they test positive 
global WORKER_CASE_ISOLATION = true  
    global removal_period = 14.0


#

# paths input files (TODO: define functions for modifying these from list in csv):

#data_dirname = pwd() * "\\network_constructor_output\\OB_2"

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
abstract type Disease_T end
abstract type Infection_T end
abstract type Immunity_T end
abstract type Vaccine_T end
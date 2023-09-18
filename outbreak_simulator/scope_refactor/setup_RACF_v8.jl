

module Setup_RACF

using Random 
using DataFrames
using CSV

#abstract type definitions 
abstract type Agent_T end
abstract type Agents_T end 

abstract type Room_T end 
abstract type Rooms_T end 

abstract type Config_T end
abstract type Pop_Input_T end

#scripts for setup of utilities for RACF model:

macro Name(arg)
    string(arg)
end

mutable struct population_input<:Pop_Input_T
    residents_str::DataFrame
    workers_G_str::DataFrame
    workers_M_str::DataFrame
    population_input() = new()
end

mutable struct run_configuration <: Config_T
    config_str::String

    R0::Float64

    delay_infection_control::Float64

    test_compliance_staff::Float64

    seed_infections::Int64
    rng_infections::AbstractRNG

    seed_testing::Int64
    rng_testing::AbstractRNG

    seed_contacts::Int64
    rng_contacts::AbstractRNG

    seed_immunity::Int64
    rng_immunity::AbstractRNG 

    dt::Float64

    contact_rate_per_resident_per_day::Float64
    bkg_contact_rate_per_resident_per_day::Float64

    w_worker_worker::Float64
    w_reg_needs::Float64
    w_high_needs::Float64
    w_res_same_room::Float64
    w_res_diff_room::Float64

    immunity::Bool
    immunity_from_dist::Bool
    uniform_immunity::Bool 

    outbreak_control::Bool
    active_outbreak::Bool

    log_mu_res::Float64
    log_sig_res::Float64
    log_mu_wG::Float64
    log_sig_wG::Float64
    log_mu_wM::Float64
    log_sig_wM::Float64

    p_test_per_day::Float64

    p_test_per_day_workers_baseline::Float64
    p_test_per_day_residents_baseline::Float64
    p_test_per_day_workers_outbreak::Float64
    p_test_per_day_residents_outbreak::Float64
    p_test_if_symptomatic::Float64

    eff_IC_worker_worker::Float64
    eff_IC_worker_resident::Float64
    eff_IC_resident_resident::Float64

    resident_lockdown::Bool
    resident_lockdown_efficacy::Float64
    worker_case_isolation::Bool
    resident_case_isolation::Bool
    resident_isolation_efficacy::Float64
    removal_period::Float64

    write_transmission_tree_flag::Bool

    resident_index_case::Bool
    worker_index_case::Bool

    pop_data_dirname::String
    residents_fname::String
    workers_G_fname::String
    workers_M_fname::String

    run_configuration() = new()# default initialiser is defined below. 
end


function setup_run_default!(config::Config_T, data_dirname::String)

    # setting up some utilities: 
    config_str = "parameter variable name, human-readable description, value \n" # initialises an empty string 

        #reproductive ratio: 
        R0 = 6.0 # default 
        description = "Reproductive ratio, calibrated to beta and timestep"
        config_line = (@Name(R0) * ", " * description * ", " * "$R0" )
        config_str = config_str*config_line*"\n"
        # assign to struct property. 
        config.R0 = R0

        # add control parameter values to config:
        # response delay  
        delay_infection_control = 2.0 # default 
            description = "delay between outbreak declaration and implementation of infection control (days)"
            config_line = (@Name(delay_infection_control) * ", " * description * ", " * "$delay_infection_control" )
            config_str = config_str*config_line*"\n"
            # assign to struct property. 
            config.delay_infection_control = delay_infection_control

        # compliance with scheduled testing 
        test_compliance_staff = 1.0 # default 
            description = "staff compliance with scheduled testing (probability of compliance with each test)"
            config_line = (@Name(test_compliance_staff) * ", " * description * ", " * "$test_compliance_staff" )
            config_str = config_str*config_line*"\n"
            config.test_compliance_staff = test_compliance_staff


        # setup random number generator: 
        seed_infections = 13 
            description = "random number seed for infection dynamics"
            config_line = (@Name(seed_infections) * ", " * description * ", " * "$seed_infections" )
            config_str = config_str*config_line*"\n"
            config.seed_infections = seed_infections 



        config.rng_infections = MersenneTwister(seed_infections)

        seed_testing = 1 #seed 1 testing as of 2022_08_29
            description = "random number seed for testing"
            config_line = (@Name(seed_testing) * ", " * description * ", " * "$seed_testing" )
            config_str = config_str*config_line*"\n"
            config.seed_testing = seed_testing


        config.rng_testing = MersenneTwister(seed_testing)

        #NOTE: another seed for network rewiring and contact selection. 
        seed_contacts = 1
            description = "random number seed for contact dynamics"
            config_line = (@Name(seed_contacts) * ", " * description * ", " * "$seed_contacts" )
            config_str = config_str*config_line*"\n"
            config.seed_contacts = seed_contacts


        config.rng_contacts = MersenneTwister(seed_contacts)

        seed_immunity = 1
            description = "random number seed for immunity status (neut titers), [only used if IMMUNITY_FROM_DIST = true]"
            config_line = (@Name(seed_immunity) * ", " * description * ", " * "$seed_immunity" )
            config_str = config_str*config_line*"\n"
            config.seed_immunity = seed_immunity

        config.rng_immunity = MersenneTwister(seed_immunity)

        #define discrete time step
        dt = 0.1
            description = "discrete time step"
            config_line = (@Name(dt) * ", " * description * ", " * "$dt" )
            config_str = config_str*config_line*"\n"
            config.dt = dt

        
        # subsection for contact rates. 
    config_str = config_str*"\n"*"*** contact rates *** \n"

        # contact rate per resident per day 
        contact_rate_per_resident_per_day = 8.0#10.0
            description = "average number of contacts sampled per resident on each day (high needs, reg. needs is 1/3 of this)"
            config_line = (@Name(contact_rate_per_resident_per_day) * ", " * description * ", " * "$contact_rate_per_resident_per_day" )
            config_str = config_str*config_line*"\n"
            config.contact_rate_per_resident_per_day = contact_rate_per_resident_per_day

        # rate of resient background contacts per day (i.e., in shared communal spaces)
        bkg_contact_rate_per_resident_per_day = 5.0
            description = "average number of random resident-resident contacts on each day"
            config_line = (@Name(bkg_contact_rate_per_resident_per_day) * ", " * description * ", " * "$bkg_contact_rate_per_resident_per_day" )
            config_str = config_str*config_line*"\n"
            config.bkg_contact_rate_per_resident_per_day = bkg_contact_rate_per_resident_per_day

    config_str = config_str*"\n"*"***relative contact durations*** \n"
        # some global control parameters: 
        # relative contact strength for high-needs residents etc: 
        w_worker_worker = 1.0
            description = "between workers"
            config_line = (@Name(w_worker_worker) * ", " * description * ", " * "$(w_worker_worker)" )
            config_str = config_str*config_line*"\n"
            config.w_worker_worker = w_worker_worker

        w_reg_needs = 2.0
            description = "between workers and residents with regular need levels"
            config_line = (@Name(w_reg_needs) * ", " * description * ", " * "$(w_reg_needs)" )
            config_str = config_str*config_line*"\n"
            config.w_reg_needs = w_reg_needs

        w_high_needs = 6.0 #2.0 # factor by which to increase contact strength between workers and high-needs residents  
            description = "between workers and residents with high needs levels"
            config_line = (@Name(w_high_needs) * ", " * description * ", " * "$(w_high_needs)" )
            config_str = config_str*config_line*"\n"
            config.w_high_needs = w_high_needs

        w_res_same_room = 10.0
            description = "between residents sharing the same room"
            config_line = (@Name(w_res_same_room) * ", " * description * ", " * "$(w_res_same_room)" )
            config_str = config_str*config_line*"\n"
            config.w_res_same_room = w_res_same_room

        w_res_diff_room = 1.0 
            description = "between residents in different rooms"
            config_line = (@Name(w_res_diff_room) * ", " * description * ", " * "$(w_res_diff_room)" )
            config_str = config_str*config_line*"\n"
            config.w_res_diff_room = w_res_diff_room

        # NOTE: as of 2022 08 30, w_res_diff_room is not used because interactions between
        # residents in different rooms is handled by random sampling at a fixed rate. 


    config_str = config_str*"\n"*"***\n"

        #flag for toggling on and off vaccine-derived protection: 
        immunity = false #default
            description = "boolean flag for whether to implement vaccine-derived protection"
            config_line = (@Name(immunity) * ", " * description * ", " * "$(immunity)" )
            config_str = config_str*config_line*"\n"
            config.immunity = immunity


        immunity_from_dist = true
            description = "boolean flag for whether to draw immunity status from file or generate it for each run"
            config_line = (@Name(immunity_from_dist) * ", " * description * ", " * "$immunity_from_dist" )
            config_str = config_str*config_line*"\n"
            config.immunity_from_dist = immunity_from_dist



            # if IMMUNITY_FROM_DIST is true, define the cohort-specific distributions 
            # in this implementation these are the normal distributions from which
            # ln(neut) values are drawn: 
            # NOTE: this is not robust code, in future versions the distributions 
            # can be setup in a way that's more flexible. 
        
        uniform_immunity = true #default 
            description = "boolean flag true if all agent types have immunity drawn from the same distribution"
            config_line = (@Name(uniform_immunity) * ", " * description * ", " * "$uniform_immunity" )
            config_str = config_str*config_line*"\n"
            config.uniform_immunity = uniform_immunity

        # flag for toggling infection control during active outbreaks (false means unmitigated)
        outbreak_control = true 
            description = "boolean flag true outbreak control measures are implemented"
            config_line = (@Name(outbreak_control) * ", " * description * ", " * "$outbreak_control" )
            config_str = config_str*config_line*"\n"
            config.outbreak_control = outbreak_control

        #note the parameters below are not used currently so I'm not adding them to the config. 
        # is an active outbreak declared? 
        active_outbreak = false #initialised as false 
            description = "boolean flag true outbreak has been detected and has not resolved"
            config_line = (@Name(active_outbreak) * ", " * description * ", " * "$active_outbreak" )
            config_str = config_str*config_line*"\n"
            config.active_outbreak = active_outbreak

    config_str = config_str*"\n"*"***immunity distribution***\n"

        log_mu_res = 0.0 #default 
            description = "mean of natural log neuts (residents)"
            config_line = (@Name(log_mu_res) * ", " * description * ", " * "$log_mu_res" )
            config_str = config_str*config_line*"\n"
            config.log_mu_res = log_mu_res

        log_sig_res = 1.5 #default 
            description = "standard deviation of natural log neuts (residents)"
            config_line = (@Name(log_sig_res) * ", " * description * ", " * "$log_sig_res" )
            config_str = config_str*config_line*"\n"
            config.log_sig_res = log_sig_res

        log_mu_wG = 0.0
            description = "mean of natural log neuts (general staff)"
            config_line = (@Name(log_mu_wG) * ", " * description * ", " * "$log_mu_wG" )
            config_str = config_str*config_line*"\n"
            config.log_mu_wG = log_mu_wG

        log_sig_wG = 1.5
            description = "standard deviation of natural log neuts (general staff)"
            config_line = (@Name(log_sig_wG) * ", " * description * ", " * "$log_sig_wG" )
            config_str = config_str*config_line*"\n"
            config.log_sig_wG = log_sig_wG

        log_mu_wM = 0.0
            description = "mean of natural log neuts (medical staff)"
            config_line = (@Name(log_mu_wM) * ", " * description * ", " * "$log_mu_wM" )
            config_str = config_str*config_line*"\n"
            config.log_mu_wM = log_mu_wM

        log_sig_wM = 1.5
            description = "standard deviation of natural log neuts (medical staff)"
            config_line = (@Name(log_sig_wM) * ", " * description * ", " * "$log_sig_wM" )
            config_str = config_str*config_line*"\n"
            config.log_sig_wM = log_sig_wM 


    config_str = config_str*"\n"*"***testing probabilities*** \n"

        # global parameter, probability of testing an agent each day 
        p_test_per_day = 0.2 # only used if testing is not sub-divided by agent type 
            description = "prob. of testing an agent each day, only used if testing is not stratified by agent type"
            config_line = (@Name(p_test_per_day) * ", " * description * ", " * "$p_test_per_day" )
            config_str = config_str*config_line*"\n"
            config.p_test_per_day = p_test_per_day

            # NOTE: implementing scheduled tests for workers, 2022 09 12
            # see Outbreak_response.jl for test scheduling. 
            # each of the scheduled tests will be subject to the probabilities 
            # noted below (note that residents have these probabilities applied 
            # every day, under the baseline or outbreak conditions ) 

        

        #NOTE: this is now assigned from control parameter for test compliance. 
        p_test_per_day_workers_baseline = config.test_compliance_staff #default 
            description = "compliance probability for worker tests (no outbreak)"
            config_line = (@Name(p_test_per_day_workers_baseline) * ", " * description * ", " * "$p_test_per_day_workers_baseline" )
            config_str = config_str*config_line*"\n"
            config.p_test_per_day_workers_baseline = p_test_per_day_workers_baseline

        p_test_per_day_residents_baseline = 0.05
            description = "robability of a resident testing if no outbreak"
            config_line = (@Name(p_test_per_day_residents_baseline) * ", " * description * ", " * "$p_test_per_day_residents_baseline" )
            config_str = config_str*config_line*"\n"
            config.p_test_per_day_residents_baseline = p_test_per_day_residents_baseline

        #NOTE: this is now assigned from control parameter for test compliance. 
        p_test_per_day_workers_outbreak = config.test_compliance_staff #default 
            description = "compliance probability for worker tests (outbreak)"
            config_line = (@Name(p_test_per_day_workers_outbreak) * ", " * description * ", " * "$p_test_per_day_workers_outbreak" )
            config_str = config_str*config_line*"\n"
            config.p_test_per_day_workers_outbreak = p_test_per_day_workers_outbreak

        p_test_per_day_residents_outbreak = 0.5 #default
            description = "probability of a resident testing during an outbreak"
            config_line = (@Name(p_test_per_day_residents_outbreak) * ", " * description * ", " * "$p_test_per_day_residents_outbreak" )
            config_str = config_str*config_line*"\n"
            config.p_test_per_day_residents_outbreak = p_test_per_day_residents_outbreak

        p_test_if_symptomatic = 1.0
            description = "probability of a symptomatic individual testing if on site"
            config_line = (@Name(p_test_if_symptomatic) * ", " * description * ", " * "$p_test_if_symptomatic" )
            config_str = config_str*config_line*"\n"
            config.p_test_if_symptomatic = p_test_if_symptomatic


    
    config_str = config_str*"\n"*"*** pairwise efficacy of infection control ***  \n"

        #efficacy of infection control for pairwise contacts: 
        # triggered by active outbreaks 
        #NOTE: including global modifier (2022 10 24, range [0, 1], [set max for each pairwise efficacy value here])
        eff_IC_worker_worker = 0.5 #* PPE_efficacy_relative_to_default
            description = "between workers"
            config_line = (@Name(eff_IC_worker_worker) * ", " * description * ", " * "$eff_IC_worker_worker" )
            config_str = config_str*config_line*"\n"
            config.eff_IC_worker_worker = eff_IC_worker_worker

        eff_IC_worker_resident = 0.9 #* PPE_efficacy_relative_to_default
            description = "between workers and residents"
            config_line = (@Name(eff_IC_worker_resident) * ", " * description * ", " * "$eff_IC_worker_resident" )
            config_str = config_str*config_line*"\n"
            config.eff_IC_worker_resident = eff_IC_worker_resident

        eff_IC_resident_resident = 0.2 #* PPE_efficacy_relative_to_default
        description = "between residents and residents"
        config_line = (@Name(eff_IC_resident_resident) * ", " * description * ", " * "$eff_IC_resident_resident" )
        config_str = config_str*config_line*"\n"
        config.eff_IC_resident_resident = eff_IC_resident_resident


    config_str = config_str*"\n"*"*** outbreak response *** \n"

        # factor by which background interaction rates are reduced during active outbreak 
        resident_lockdown = true 
            description = "flag for whether contact rates between non-isolated residents will be reduced"
            config_line = (@Name(resident_lockdown) * ", " * description * ", " * "$resident_lockdown" )
            config_str = config_str*config_line*"\n"
            config.resident_lockdown = resident_lockdown

        resident_lockdown_efficacy = 0.5 # representing discretionary changes 
            description = "fraction by which background contact rates are reduced for residents during outbreaks"
            config_line = (@Name(resident_lockdown_efficacy) * ", " * description * ", " * "$resident_lockdown_efficacy" )
            config_str = config_str*config_line*"\n"
            config.resident_lockdown_efficacy = resident_lockdown_efficacy

        # duration for which a worker is removed from the facility if they test positive 
        worker_case_isolation = true
            description = "flag whether or not to furlough staff who test positive"
            config_line = (@Name(worker_case_isolation) * ", " * description * ", " * "$worker_case_isolation" )
            config_str = config_str*config_line*"\n"
            config.worker_case_isolation = worker_case_isolation

        resident_case_isolation = true # 2022 09 12
            description = "flag whether or not to isolate residents who test positive"
            config_line = (@Name(resident_case_isolation) * ", " * description * ", " * "$resident_case_isolation" )
            config_str = config_str*config_line*"\n"
            config.resident_case_isolation = resident_case_isolation

        resident_isolation_efficacy = 0.90
            description = "reduction in background contact rate for isolated residents"
            config_line = (@Name(resident_isolation_efficacy) * ", " * description * ", " * "$resident_isolation_efficacy" )
            config_str = config_str*config_line*"\n"
            config.resident_isolation_efficacy = resident_isolation_efficacy


        removal_period = 7.0 # changed from 14 to 7 on 2022 09 12
            description = "period of time (days) for isolation or furlough"
            config_line = (@Name(removal_period) * ", " * description * ", " * "$removal_period" )
            config_str = config_str*config_line*"\n"
            config.removal_period = removal_period

    config_str = config_str*"\n"*"***\n"
    #flag telling the system to write the transmission tree and detections to file: 
    write_transmission_tree_flag = false
        description = "flag telling the system to write the transmission tree to file"
        config_line = (@Name(write_transmission_tree_flag) * ", " * description * ", " * "$write_transmission_tree_flag" )
        config_str = config_str*config_line*"\n"
        config.write_transmission_tree_flag = write_transmission_tree_flag

    config_str = config_str*"\n"*"*** type of index case *** \n"
    # global control parameter: worker index case vs. resident index case (if both are true it can be either)
    resident_index_case = true
        description = "worker index case vs. resident index case (if both are true it can be either)"
        config_line = (@Name(resident_index_case) * ", " * description * ", " * "$resident_index_case" )
        config_str = config_str*config_line*"\n"
        config.resident_index_case = resident_index_case

    worker_index_case = true
        description = "worker index case vs. resident index case (if both are true it can be either)"
        config_line = (@Name(worker_index_case) * ", " * description * ", " * "$worker_index_case" )
        config_str = config_str*config_line*"\n"
        config.worker_index_case = worker_index_case



    # paths input files (TODO: define functions for modifying these from list in csv):

    config_str = config_str*"\n"*"***population data location***\n"
        config_line = data_dirname 
        config_str = config_str*config_line*"\n"
        config.pop_data_dirname = data_dirname

    #E_list_fname = data_dirname * "\\Edge_lists_test.csv"
    config.residents_fname = data_dirname * "\\residents_for_Julia_test.csv"
    config.workers_G_fname = data_dirname * "\\workers_G_for_Julia_test.csv"
    config.workers_M_fname = data_dirname * "\\workers_M_for_Julia_test.csv"
    #N_list_fname = data_dirname * "\\Neighbour_lists_test.csv"

    config.config_str = config_str

end

function update_config_record!(config::Config_T)

    config_str = "parameter variable name, human-readable description, value \n" # initialises an empty string 

        #reproductive ratio: 
        R0 = config.R0
            description = "Reproductive ratio, calibrated to beta and timestep"
            config_line = (@Name(R0) * ", " * description * ", " * "$R0" )
            config_str = config_str*config_line*"\n"

        # add control parameter values to config:
        # response delay  
        delay_infection_control = config.delay_infection_control
            description = "delay between outbreak declaration and implementation of infection control (days)"
            config_line = (@Name(delay_infection_control) * ", " * description * ", " * "$delay_infection_control" )
            config_str = config_str*config_line*"\n"

        # compliance with scheduled testing 
        test_compliance_staff = config.test_compliance_staff
            description = "staff compliance with scheduled testing (probability of compliance with each test)"
            config_line = (@Name(test_compliance_staff) * ", " * description * ", " * "$test_compliance_staff" )
            config_str = config_str*config_line*"\n"


        # setup random number generator: 
        seed_infections = config.seed_infections
            description = "random number seed for infection dynamics"
            config_line = (@Name(seed_infections) * ", " * description * ", " * "$seed_infections" )
            config_str = config_str*config_line*"\n"

        seed_testing = config.seed_testing
            description = "random number seed for testing"
            config_line = (@Name(seed_testing) * ", " * description * ", " * "$seed_testing" )
            config_str = config_str*config_line*"\n"


        #NOTE: another seed for network rewiring and contact selection. 
        seed_contacts = config.seed_contacts
            description = "random number seed for contact dynamics"
            config_line = (@Name(seed_contacts) * ", " * description * ", " * "$seed_contacts" )
            config_str = config_str*config_line*"\n"

        seed_immunity = config.seed_immunity
            description = "random number seed for immunity status (neut titers), [only used if IMMUNITY_FROM_DIST = true]"
            config_line = (@Name(seed_immunity) * ", " * description * ", " * "$seed_immunity" )
            config_str = config_str*config_line*"\n"

        #define discrete time step
        dt = config.dt
            description = "discrete time step"
            config_line = (@Name(dt) * ", " * description * ", " * "$dt" )
            config_str = config_str*config_line*"\n"

        
        # subsection for contact rates. 
    config_str = config_str*"\n"*"*** contact rates *** \n"

        # contact rate per resident per day 
        contact_rate_per_resident_per_day = config.contact_rate_per_resident_per_day
            description = "average number of contacts sampled per resident on each day (high needs, reg. needs is 1/3 of this)"
            config_line = (@Name(contact_rate_per_resident_per_day) * ", " * description * ", " * "$contact_rate_per_resident_per_day" )
            config_str = config_str*config_line*"\n"

        # rate of resient background contacts per day (i.e., in shared communal spaces)
        bkg_contact_rate_per_resident_per_day = config.bkg_contact_rate_per_resident_per_day
            description = "average number of random resident-resident contacts on each day"
            config_line = (@Name(bkg_contact_rate_per_resident_per_day) * ", " * description * ", " * "$bkg_contact_rate_per_resident_per_day" )
            config_str = config_str*config_line*"\n"

    config_str = config_str*"\n"*"***relative contact durations*** \n"
        # some global control parameters: 
        # relative contact strength for high-needs residents etc: 
        w_worker_worker = config.w_worker_worker
            description = "between workers"
            config_line = (@Name(w_worker_worker) * ", " * description * ", " * "$(w_worker_worker)" )
            config_str = config_str*config_line*"\n"

        w_reg_needs = config.w_reg_needs
            description = "between workers and residents with regular need levels"
            config_line = (@Name(w_reg_needs) * ", " * description * ", " * "$(w_reg_needs)" )
            config_str = config_str*config_line*"\n"

        w_high_needs = config.w_high_needs 
            description = "between workers and residents with high needs levels"
            config_line = (@Name(w_high_needs) * ", " * description * ", " * "$(w_high_needs)" )
            config_str = config_str*config_line*"\n"

        w_res_same_room = config.w_res_same_room
            description = "between residents sharing the same room"
            config_line = (@Name(w_res_same_room) * ", " * description * ", " * "$(w_res_same_room)" )
            config_str = config_str*config_line*"\n"

        w_res_diff_room = config.w_res_diff_room
            description = "between residents in different rooms"
            config_line = (@Name(w_res_diff_room) * ", " * description * ", " * "$(w_res_diff_room)" )
            config_str = config_str*config_line*"\n"

        # NOTE: as of 2022 08 30, w_res_diff_room is not used because interactions between
        # residents in different rooms is handled by random sampling at a fixed rate. 


    config_str = config_str*"\n"*"***\n"

        #flag for toggling on and off vaccine-derived protection: 
        immunity = config.immunity
            description = "boolean flag for whether to implement vaccine-derived protection"
            config_line = (@Name(immunity) * ", " * description * ", " * "$(immunity)" )
            config_str = config_str*config_line*"\n"


        immunity_from_dist = config.immunity_from_dist
            description = "boolean flag for whether to draw immunity status from file or generate it for each run"
            config_line = (@Name(immunity_from_dist) * ", " * description * ", " * "$immunity_from_dist" )
            config_str = config_str*config_line*"\n"


            # if IMMUNITY_FROM_DIST is true, define the cohort-specific distributions 
            # in this implementation these are the normal distributions from which
            # ln(neut) values are drawn: 
            # NOTE: this is not robust code, in future versions the distributions 
            # can be setup in a way that's more flexible. 
        
        uniform_immunity = config.uniform_immunity
            description = "boolean flag true if all agent types have immunity drawn from the same distribution"
            config_line = (@Name(uniform_immunity) * ", " * description * ", " * "$uniform_immunity" )
            config_str = config_str*config_line*"\n"

        # flag for toggling infection control during active outbreaks (false means unmitigated)
        outbreak_control = config.outbreak_control
            description = "boolean flag true outbreak control measures are implemented"
            config_line = (@Name(outbreak_control) * ", " * description * ", " * "$outbreak_control" )
            config_str = config_str*config_line*"\n"

        #note the parameters below are not used currently so I'm not adding them to the config. 
        # is an active outbreak declared? 
        active_outbreak = config.active_outbreak
            description = "boolean flag true outbreak has been detected and has not resolved"
            config_line = (@Name(active_outbreak) * ", " * description * ", " * "$active_outbreak" )
            config_str = config_str*config_line*"\n"

    config_str = config_str*"\n"*"***immunity distribution***\n"

        log_mu_res = config.log_mu_res
            description = "mean of natural log neuts (residents)"
            config_line = (@Name(log_mu_res) * ", " * description * ", " * "$log_mu_res" )
            config_str = config_str*config_line*"\n"

        log_sig_res = config.log_sig_res
            description = "standard deviation of natural log neuts (residents)"
            config_line = (@Name(log_sig_res) * ", " * description * ", " * "$log_sig_res" )
            config_str = config_str*config_line*"\n"

        log_mu_wG = config.log_mu_wG
            description = "mean of natural log neuts (general staff)"
            config_line = (@Name(log_mu_wG) * ", " * description * ", " * "$log_mu_wG" )
            config_str = config_str*config_line*"\n"

        log_sig_wG = config.log_sig_wG
            description = "standard deviation of natural log neuts (general staff)"
            config_line = (@Name(log_sig_wG) * ", " * description * ", " * "$log_sig_wG" )
            config_str = config_str*config_line*"\n"
            

        log_mu_wM = config.log_mu_wM
            description = "mean of natural log neuts (medical staff)"
            config_line = (@Name(log_mu_wM) * ", " * description * ", " * "$log_mu_wM" )
            config_str = config_str*config_line*"\n"

        log_sig_wM = config.log_sig_wM
            description = "standard deviation of natural log neuts (medical staff)"
            config_line = (@Name(log_sig_wM) * ", " * description * ", " * "$log_sig_wM" )
            config_str = config_str*config_line*"\n"


    config_str = config_str*"\n"*"***testing probabilities*** \n"

        # global parameter, probability of testing an agent each day 
        p_test_per_day = config.p_test_per_day
            description = "prob. of testing an agent each day, only used if testing is not stratified by agent type"
            config_line = (@Name(p_test_per_day) * ", " * description * ", " * "$p_test_per_day" )
            config_str = config_str*config_line*"\n"
            

            # NOTE: implementing scheduled tests for workers, 2022 09 12
            # see Outbreak_response.jl for test scheduling. 
            # each of the scheduled tests will be subject to the probabilities 
            # noted below (note that residents have these probabilities applied 
            # every day, under the baseline or outbreak conditions ) 

        

        #NOTE: this is now assigned from control parameter for test compliance. 
        p_test_per_day_workers_baseline = config.p_test_per_day_workers_baseline
            description = "compliance probability for worker tests (no outbreak)"
            config_line = (@Name(p_test_per_day_workers_baseline) * ", " * description * ", " * "$p_test_per_day_workers_baseline" )
            config_str = config_str*config_line*"\n"
            

        p_test_per_day_residents_baseline = config.p_test_per_day_residents_baseline
            description = "robability of a resident testing if no outbreak"
            config_line = (@Name(p_test_per_day_residents_baseline) * ", " * description * ", " * "$p_test_per_day_residents_baseline" )
            config_str = config_str*config_line*"\n"
            

        #NOTE: this is now assigned from control parameter for test compliance. 
        p_test_per_day_workers_outbreak = config.p_test_per_day_workers_outbreak
            description = "compliance probability for worker tests (outbreak)"
            config_line = (@Name(p_test_per_day_workers_outbreak) * ", " * description * ", " * "$p_test_per_day_workers_outbreak" )
            config_str = config_str*config_line*"\n"
            

        p_test_per_day_residents_outbreak = config.p_test_per_day_residents_outbreak
            description = "probability of a resident testing during an outbreak"
            config_line = (@Name(p_test_per_day_residents_outbreak) * ", " * description * ", " * "$p_test_per_day_residents_outbreak" )
            config_str = config_str*config_line*"\n"
            

        p_test_if_symptomatic = config.p_test_if_symptomatic
            description = "probability of a symptomatic individual testing if on site"
            config_line = (@Name(p_test_if_symptomatic) * ", " * description * ", " * "$p_test_if_symptomatic" )
            config_str = config_str*config_line*"\n"
            



    config_str = config_str*"\n"*"*** pairwise efficacy of infection control ***  \n"

        #efficacy of infection control for pairwise contacts: 
        # triggered by active outbreaks 
        #NOTE: including global modifier (2022 10 24, range [0, 1], [set max for each pairwise efficacy value here])
        eff_IC_worker_worker = config.eff_IC_worker_worker
            description = "between workers"
            config_line = (@Name(eff_IC_worker_worker) * ", " * description * ", " * "$eff_IC_worker_worker" )
            config_str = config_str*config_line*"\n"
            

        eff_IC_worker_resident = config.eff_IC_worker_resident
            description = "between workers and residents"
            config_line = (@Name(eff_IC_worker_resident) * ", " * description * ", " * "$eff_IC_worker_resident" )
            config_str = config_str*config_line*"\n"
            

        eff_IC_resident_resident = config.eff_IC_resident_resident
            description = "between residents and residents"
            config_line = (@Name(eff_IC_resident_resident) * ", " * description * ", " * "$eff_IC_resident_resident" )
            config_str = config_str*config_line*"\n"
        


    config_str = config_str*"\n"*"*** outbreak response *** \n"

        # factor by which background interaction rates are reduced during active outbreak 
        resident_lockdown = config.resident_lockdown
            description = "flag for whether contact rates between non-isolated residents will be reduced"
            config_line = (@Name(resident_lockdown) * ", " * description * ", " * "$resident_lockdown" )
            config_str = config_str*config_line*"\n"

        resident_lockdown_efficacy = config.resident_lockdown_efficacy
            description = "fraction by which background contact rates are reduced for residents during outbreaks"
            config_line = (@Name(resident_lockdown_efficacy) * ", " * description * ", " * "$resident_lockdown_efficacy" )
            config_str = config_str*config_line*"\n"
            

        # duration for which a worker is removed from the facility if they test positive 
        worker_case_isolation = config.worker_case_isolation
            description = "flag whether or not to furlough staff who test positive"
            config_line = (@Name(worker_case_isolation) * ", " * description * ", " * "$worker_case_isolation" )
            config_str = config_str*config_line*"\n"
            

        resident_case_isolation = config.resident_case_isolation
            description = "flag whether or not to isolate residents who test positive"
            config_line = (@Name(resident_case_isolation) * ", " * description * ", " * "$resident_case_isolation" )
            config_str = config_str*config_line*"\n"

        resident_isolation_efficacy = config.resident_isolation_efficacy
            description = "reduction in background contact rate for isolated residents"
            config_line = (@Name(resident_isolation_efficacy) * ", " * description * ", " * "$resident_isolation_efficacy" )
            config_str = config_str*config_line*"\n"

        removal_period = config.removal_period
            description = "period of time (days) for isolation or furlough"
            config_line = (@Name(removal_period) * ", " * description * ", " * "$removal_period" )
            config_str = config_str*config_line*"\n"

        #2022 09 21 setup variables moved from Main: 
        # flag for toggling infection control during active outbreaks (false means unmitigated)
        outbreak_control = config.outbreak_control
            description = "flag whether to apply any outbreak control measures"
            config_line = (@Name(outbreak_control) * ", " * description * ", " * "$outbreak_control" )
            config_str = config_str*config_line*"\n"


    config_str = config_str*"\n"*"***\n"
        #flag telling the system to write the transmission tree and detections to file: 
        write_transmission_tree_flag = config.write_transmission_tree_flag
            description = "flag telling the system to write the transmission tree to file"
            config_line = (@Name(write_transmission_tree_flag) * ", " * description * ", " * "$write_transmission_tree_flag" )
            config_str = config_str*config_line*"\n"

    config_str = config_str*"\n"*"*** type of index case *** \n"
        # global control parameter: worker index case vs. resident index case (if both are true it can be either)
        resident_index_case = config.resident_index_case
            description = "worker index case vs. resident index case (if both are true it can be either)"
            config_line = (@Name(resident_index_case) * ", " * description * ", " * "$resident_index_case" )
            config_str = config_str*config_line*"\n"

        worker_index_case = config.worker_index_case
            description = "worker index case vs. resident index case (if both are true it can be either)"
            config_line = (@Name(worker_index_case) * ", " * description * ", " * "$worker_index_case" )
            config_str = config_str*config_line*"\n"

    # paths input files 
    config_str = config_str*"\n"*"***population data location***\n"
        config_line = config.pop_data_dirname #note this is a global variable, might not be able to see it
        config_str = config_str*config_line*"\n"

    
    config.config_str = config_str

end



function read_in_population_data!(pop_input::Pop_Input_T, config::Config_T)
    pop_input.residents_str = DataFrame(CSV.File(config.residents_fname, delim = "\t"))
    pop_input.workers_G_str = DataFrame(CSV.File(config.workers_G_fname, delim = "\t"))
    pop_input.workers_M_str = DataFrame(CSV.File(config.workers_M_fname, delim = "\t"))
    return
end


function write_config_details(config::Config_T, outdir::String)

    # write config str to output file. (note - applies global variable output_dir_fac)
    config_file = open("$(outdir)\\config.txt", "w")
    write(config_file, config.config_str)
    close(config_file)

    return

end


function set_immunity_dist!(config::Config_T)

    # default values (taken from facility-specific distributions or aggregates thereof.)
    # NOTE: larger sigma values reflect heterogeneity in timing and vaccine type in real facilities. 
    # NOTE: for reference, mu of approx. -1.5 corresponds to (for example) Pfizer dose 3 after 12 weeks of waning. 
    config.log_mu_res = -1.47
    config.log_sig_res = 1.34

    config.log_mu_wG = -1.79
    config.log_sig_wG = 1.43

    config.log_mu_wM = -1.64
    config.log_sig_wM = 1.51

    return

end


function apply_seed_offset!(config::Config_T, seed_offset::Int64)

    config.seed_contacts += seed_offset
    config.seed_testing += seed_offset 
    config.seed_infections += seed_offset 
    config.seed_immunity += seed_offset 

    config.rng_contacts = MersenneTwister(seed_contacts)
    config.rng_testing = MersenneTwister(seed_testing)
    config.rng_infections = MersenneTwister(seed_infections)
    config.rng_immunity = MersenneTwister(seed_immunity)

    return

end

end
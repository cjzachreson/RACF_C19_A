


module Diseases_RACF

using Distributions
using Random
using Dates
using Main.Setup_RACF # need to import configuration variables. 

abstract type Disease_T end
abstract type Diseases_T end 
abstract type Infection_T end
abstract type Immunity_Profile_T end
abstract type Vaccine_T end


# immutable stuct for general disease parameters, defined 
# for each different pathogen: 
mutable struct disease_params <: Disease_T

    ## beta, global transmission scaler
    beta::Float64 # 1
    
    # symptomatic fraction
    p_asymp::Float64 # 2

    #latent period 
    latent_period::Float64 #here initialised as dt (for asynchronous update) # 3

    ## incubation period, log-normal 
    inc_mu::Float64 # 4
    inc_sig::Float64 # 5
    inc_dist::LogNormal{Float64} # 6

    ## post-incubation (recovery) period 
    rec_min::Float64 # 7
    rec_max::Float64 # 8
    rec_dist::Uniform{Float64} # 9
    
    ## distribution of infectiousness 
    b_dispersion::Float64 # 10
    b_scale::Float64 # 11
    b_dist::Gamma{Float64} # 12

    Vmax::Float64 # scaling factor for growth and decline of viral load # 13
    inc_plat_fac::Float64 # proportion of incubation period in plateau viral load # 14
    rec_plat_fac::Float64 # proportion of recovery period in plateau viral load # 15

    #give the pathogen a name, for example
    name::String # 16

    #ranges of test sensitivity function parameters: 
    c_lims::Vector{Float64} # 17
    b1_lims::Vector{Float64} # 18
    b2_lims::Vector{Float64} # 19
    b3_lims::Vector{Float64} # 20

    #constructor (used below for Delta veriant SARS-CoV-2)

    # TODO: set this up as a null constructor, with default constructor defined separately
    disease_params() = new()

end

# default parameters are delta-variant like
function set_disease_params_default!(params::Disease_T, config::Setup_RACF.Config_T)

    # initialise disease (Delta variant, SARS-CoV-2) 
    # TODO: update for omicron
    #R0 = config.R0#6.0
    # note: scaling this down because I removed dt from the FoI computation 
    # in Transmission_Dynamics.jl, this was a carryover from when I was 
    # implementing this as frequency-dependent transmission in a previous project
    # in this implementation, the effect of the timestep is accounted for in the 
    # contact rate, not in the force of infection. Including a constant of 0.1 
    # here to make sure calibration is stable. 
    
    #old code: beta = R0/3.94 
    # new code:
    params.beta = (config.R0/3.94) * 0.1 # 2022 10 20
    # note this is set to R0 = 6 in config, by default. 
    # if multiple disease are used, this will need to be replaced with 
    # a vector in the Setup module. 

    #note: 3.94 is the scaling factor 
    # determined during the calibration of the hotel quarantine model
    #pre-dating this work 
    # TODO: verify calibration 
    params.p_asymp = 0.33
    params.latent_period = config.dt
    params.inc_mu = 1.62
    params.inc_sig = 0.418
    params.rec_min = 5.0
    params.rec_max = 10.0
    params.b_dispersion = 0.15
    params.Vmax = 7.0
    params.inc_plat_fac = 0.1
    params.rec_plat_fac = 0.0

    #test sensitivity function parameters: 
    params.c_lims = [1.0, 5.11]
    params.b1_lims = [0.8, 2.31]
    params.b2_lims = [1.26, 3.47]
    params.b3_lims = [1.05, 1.14]

    ## incubation period, log-normal 
    params.inc_dist = LogNormal(params.inc_mu, params.inc_sig) #inc_dist::LogNormal{Float64}

    ## post-incubation (recovery) period 
    params.rec_dist = Uniform(params.rec_min, params.rec_max) #rec_dist::Uniform{Float64}
        
    ## distribution of infectiousness 
    params.b_scale = params.beta / params.b_dispersion  #b_scale::Float64
    params.b_dist = Gamma(params.b_dispersion, params.beta/params.b_dispersion) #b_dist::Gamma{Float64}

    params.name = "Default"

end
# mutable struct for specific infection parameters, for 
# infection of an agent with a specific pathogen. 
mutable struct infection <: Infection_T 
    ## gets created during call to infect_agent!()
    ## stores properties and time-dependent characteristics of individual infections

    #pathogen object: 
    pathogen::Disease_T
    #pathogen name:
    pathogen_name::String
    # time since infection
    t_infected::Float64

    #latent period 
    t_latent::Float64

    # incubation period (time between exposure and symptom expression)
    t_inc::Float64
    # cdf of t_inc from incubation period dist, used to enforce correlations with test sensitivity and infectiousness 
    q_inc::Float64 
    # recovery period (time between peak viral load and recovery)
    t_rec::Float64
    # symptom expression (bool symptomatic or not) "Will they express symptoms?"
    symptomatic::Bool 
    # symptom expression (bool expressing symptoms or not) "are they currently expressing symptoms?"
    expressing_symptoms::Bool 

    # force of infection (i.e., baseline infection transmission rate as a function of time)
    beta_t::Float64 
    
    ## parameters determining force of infection
    #maximum infectiousness (peaks at symptom onset)
    beta_max::Float64 
    beta_min::Float64 #Def: beta_max/Vmax

    k_inc::Float64 #rate of exponential increase of infectiousness during incubation 
    k_rec::Float64 #rate of exponential decrease of infectiousness during recovery 

    # parameters determining detection probability
    # RAT detection probability, NOTE: decide whether this is an infection property, or something else... 
    # for now, I'm going to include these as infection properties, but they could potentially be implemented 
    # separately... 
    test_sensitivity::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    changepoint::Float64

    #null constructor: 
    infection() = new()
 
end

function set_infection_default!(infection::Infection_T, 
                                pathogen::Disease_T, 
                                config::Setup_RACF.Config_T)
    
    #moved these 4 up here so rng order is same as original implementation. 
    infection.t_inc = rand(config.rng_infections, pathogen.inc_dist)
    infection.q_inc = cdf(pathogen.inc_dist, infection.t_inc)
    infection.t_rec = rand(config.rng_infections, pathogen.rec_dist)#t_rec::Float64
    infection.beta_max = rand(config.rng_infections, pathogen.b_dist)#beta_max::Float64 
    #

    ## pathogen reference
    infection.pathogen = pathogen
    ## name of pathogen 
    infection.pathogen_name = pathogen.name # pathogen::String 
    infection.t_infected = 0.0 #t_infected::Float64
    ## latent period
    infection.t_latent = pathogen.latent_period # TODO: draw from distribution as with other intervals
    ## incubation period (time between exposure and symptom expression)
    #infection.t_inc = rand(config.rng_infections, pathogen.inc_dist)#t_inc::Float64
    #infection.q_inc = cdf(pathogen.inc_dist, infection.t_inc)#q_inc::Float64 NOTE: not sure if it will let me do this, need to check. *** Nope - need to pass as input to constructor. 
    ## recovery period (time between peak viral load and recovery)
    #infection.t_rec = rand(config.rng_infections, pathogen.rec_dist)#t_rec::Float64
    ## symptom expression (bool symptomatic or not) "Will they express symptoms?"
    infection.symptomatic = rand(config.rng_infections) < (1.0 - pathogen.p_asymp)#symptomatic::Bool 
    ## symptom expression (bool expressing symptoms or not) "are they currently expressing symptoms?"
    infection.expressing_symptoms = false #expressing_symptoms::Bool 
    ## force of infection (i.e., baseline infection transmission rate as a function of time)
    infection.beta_t = 0.0#beta_t::Float64 
    
    ### parameters determining force of infection
    ##maximum infectiousness (peaks at symptom onset)
    #infection.beta_max = rand(config.rng_infections, pathogen.b_dist)#beta_max::Float64 
    infection.beta_min = infection.beta_max / pathogen.Vmax #beta_min::Float64 #Def: beta_max/Vmax #[NOTE: again, not sure if Julia will let me initialise in this way (does it know the value of beta_max?)] ***

    infection.k_inc = compute_kinc(infection.t_inc, pathogen)#k_inc::Float64 #rate of exponential increase of infectiousness during incubation *** (may need to initialise these after construction)
    infection.k_rec = compute_krec(infection.t_rec, pathogen)#k_rec::Float64 #rate of exponential decrease of infectiousness during recovery ***

    ## parameters determining detection probability
    ## RAT detection probability, NOTE: decide whether this is an infection property, or something else... 
    ## for now, I'm going to include these as infection properties, but they could potentially be implemented 
    ## separately... 
    
    infection.test_sensitivity = 0.0#test_sensitivity::Float64 (initially set to 0)
    
    infection.b1 = pathogen.b1_lims[1] + 
                    rand(config.rng_infections) * (pathogen.b1_lims[2] - pathogen.b1_lims[1])#b1::Float64
    
    infection.b2 = pathogen.b2_lims[1] + 
                    (1.0 - infection.q_inc) * (pathogen.b2_lims[2] - pathogen.b2_lims[1])#b2::Float64 ***quantile matching with incubation period. 
    
    infection.b3 = pathogen.b3_lims[1] + 
                    rand(config.rng_infections) * (pathogen.b3_lims[2] - pathogen.b3_lims[1])#b3::Float64
    
    infection.changepoint = infection.t_inc - 
                    (infection.q_inc * (pathogen.c_lims[2] - pathogen.c_lims[1]))#changepoint::Float64
    
end

# infection update function for time-dependent parameters: 
function update_infection!(infection::Infection_T, dt::Float64)::Bool 
    # return value is recovery flag 
    recovered = false

    # increase time since infection
    infection.t_infected += dt

    # check recovery status
    # what to do when an individual recovers? 
    # agent's immunity status should be updated, and the infection should 
    # be removed from memory. 
    if infection.t_infected > (infection.t_inc + infection.t_rec) 
        recovered = true
    end
    
    # check symptom expression
    if (infection.t_infected > infection.t_inc) && infection.symptomatic 
        infection.expressing_symptoms = true
    end 

    # update force of infection 
    # TODO: allow asymptomatic cases to be less infectious 
    compute_FoI!(infection)

    # update detection probability
    #TODO: don't need to do this every update, only if a test is about to be performed. 
    compute_test_sensitivity!(infection)

    return recovered 

end


##initialisation helper functions for computing some derived infection parameters: 
# maps t_plat_fac, t_inc (incubation period of agent i) 
#and Vmax to the growth rate of infectiousness during incubation
function compute_kinc(t_inc::Float64, pathogen::Disease_T) 
    t_plat_inc = pathogen.inc_plat_fac * t_inc
    k_inc = log(pathogen.Vmax)/(t_inc - t_plat_inc)
    return k_inc
end

# maps to the decay rate of infectiousness during 'recovery' i.e.,
# after symptom onset
function compute_krec(t_rec::Float64, pathogen::Disease_T) #t_rec 
    t_plat_rec = pathogen.rec_plat_fac * t_rec
    k_rec = log(1.0/pathogen.Vmax) / (t_rec - t_plat_rec)
    return k_rec
end


## update functions for computing time-dependent infection properties 
# compute FoI 
function compute_FoI!(infection::Infection_T)

    t_pt = 0.0 #local variable used for piecewise timepoint 
    #infection.beta_t = ?  #update force of infection based on t_infected 

    #check if latent: 
    if infection.t_infected <= infection.t_latent #latent 
        t_pt = infection.t_infected
        infection.beta_t = 0.0
    #check if agent is incubating 
    elseif infection.t_infected <= infection.t_inc # incubating
        t_pt = infection.t_infected
        infection.beta_t = 
        infection.beta_min * exp(( infection.k_inc * t_pt) / infection.pathogen.Vmax)
    # check if recovering 
    elseif infection.t_infected <= (infection.t_rec + infection.t_inc) # recovering (i.e., symptomatic period if symptomatic)
        t_pt = infection.t_infected - infection.t_inc 
        infection.beta_t = infection.beta_max * exp(infection.k_rec * t_pt)
    #if none of the above, then recovered. 
    else  
        infection.beta_t = 0.0
    end

    #plateau implemented as cuttoff: 
    if infection.beta_t > infection.beta_max
        infection.beta_t = infection.beta_max
    end

end

# compute test sensitivity 
function compute_test_sensitivity!(infection::Infection_T)

    tau = infection.t_infected - infection.changepoint

    if tau < 0
        infection.test_sensitivity = 
            1.0 / (1.0 + exp(-1.0 * 
            (infection.b1 + infection.b2 * tau )))
    else
        infection.test_sensitivity = 
            1.0 / (1.0 + exp(-1.0 * 
            (infection.b1 + infection.b2 * tau + (-1.0 * infection.b2*infection.b3*tau) )))
    end

end



#***** Define dictionary of pathogens ******
# dictionary of diseases 
#TODO : read these in from parameter files instead of hard-coding. 

mutable struct diseases <: Diseases_T

    dict::Dict{String, Disease_T} #maping from name to disease parameters
    names::Vector{String} # names of the different types of infections
    n::Int64 #number of different types of infections


    diseases() = new() #null constructor

end

function set_disease_dict!(diseases::Diseases_T, config::Setup_RACF.Config_T)

    diseases.dict = Dict{String, Disease_T}() #[disease_name] = Delta_variant 

    #TODO: parameters are still Delta variant- adjust for omicron. 
    diseases.names = ["Default"]#["Delta Variant"]
    diseases.n = length(diseases.names) 

    #NOTE: this is setup for multiple diseases, but only one is implemented. 
    # the parameters below are fixed. For multistrain implementations, 
    # TODO: this would need parameter vectors setup before the loop. 
    for d = 1:diseases.n

        disease_name = diseases.names[d]

        disease_i = disease_params()
        set_disease_params_default!(disease_i, config::Setup_RACF.Config_T)
        # any modifications to the parameters for different
        # pathogens or strains would be made here. 

        diseases.dict[disease_name] = disease_i

    end

end


#***** Define dictionary of immunity types ******

#what defines an immunity type? 
# for our purposes, we've already determined NAT values
# for each individual at the start of an outbreak, so we don't 
# need to worry about vaccination history 
mutable struct Immunity_Profile <:Immunity_Profile_T 

    vaccination_history::Dict{DateTime, String} 

    infection_history::Dict{DateTime, Disease_T}

    NAT_peak::Float64 #scaled peak neuts. 
    dt_peak::Float64 # time since date of peak neuts. (duration of waning) 
    NAT_decay_rate::Float64 #exponential decay rate of NAT levels from peak 
    NAT_t::Float64 #scaled neutralising antibody titers at current time 

    # string key values correspond to names of pathogens 
    protection_Infection::Dict{String, Float64}
    protection_Symptoms::Dict{String, Float64}
    protection_Transmission::Dict{String, Float64}
    protection_Death::Dict{String, Float64} 

    #constructor: 
    Immunity_Profile() = new(
        Dict{DateTime, String}(), #vaccination history 
        Dict{DateTime, Disease_T}(), #infection history
        0.0, # peak NAT (initialised to 0)
        0.0, # delay since peak (initialised to 0)
        0.0, # NAT decay rate (initialised to 0)
        0.0, # NAT level at time t 
        #protection dicts: 
        Dict{String, Float64}(), # Infection
        Dict{String, Float64}(), # Symptoms 
        Dict{String, Float64}(), # Onward transmission
        Dict{String, Float64}()  # Death 
    )

end

#consider these placeholders for now, efficacy (protection) values are imputed directly. 
# in this version but in future versions they may be computed during a run, and can be dynamic. 

# NOTE: these are currently hacked in, the c50 parameters should be contained elsewhere
# either in the disease dictionary or in the immunity provile, this implementation is 
# specific to the scenarios involving SARS-CoV-2 (Omicron)

# copied over from my matlab implementation: 
#c50_hosp = -1.206;
#c50_death = -1.184;
#c50_acquisition = -0.4717;
#c50_transmission = 0.01846;
#c50_symptoms = -0.6349;


#NOTE: these take the raw, relative neut values as input, not the log-10 
function Efficacy_infection_Omicron(neuts::Float64)::Float64
    c50 = -0.4717

    k = exp(1.707)
    log10_neuts = log10(neuts)
    eff = 1.0 / (1.0 + exp(-1.0 * k *(log10_neuts - c50)))
    return eff 
end

function Efficacy_OT_Omicron(neuts::Float64)::Float64
    
    c50 = 0.01846

    k = exp(1.707)
    log10_neuts = log10(neuts)
    eff = 1.0 / (1.0 + exp(-1.0 * k *(log10_neuts - c50)))
    return eff 
end

function Efficacy_Symptoms_Omicron(neuts::Float64)::Float64

    c50 = -0.6349

    k = exp(1.707)
    log10_neuts = log10(neuts)
    eff = 1.0 / (1.0 + exp(-1.0 * k *(log10_neuts - c50)))
    return eff 
end

function Efficacy_Death_Omicron(neuts::Float64)::Float64

    c50 = -1.184

    k = exp(1.707)
    log10_neuts = log10(neuts)
    eff = 1.0 / (1.0 + exp(-1.0 * k *(log10_neuts - c50)))
    return eff 
end



end
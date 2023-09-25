# immutable stuct for general disease parameters, defined 
# for each different pathogen: 
struct disease_params <: Disease_T

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

    disease_params(beta, 
                   p_asymp, 
                   latent_period, 
                   inc_mu, 
                   inc_sig, 
                   rec_min, 
                   rec_max, 
                   b_dispersion, 
                   V_max, 
                   inc_plat_fac, 
                   rec_plat_fac, 
                   name, 
                   c_lims, 
                   b1_lims, 
                   b2_lims, 
                   b3_lims ) = 
            new(
            ## beta, global transmission scaler
            beta,#R0 / 3.94,    #beta::Float64  # 1
    
            # asymptomatic fraction
            p_asymp,    #p_asymp::Float64  

            #latent period 
            latent_period,    #latent_period::Float64 #here initialised as dt (for asynchronous update) 

            ## incubation period, log-normal 
            inc_mu,    #inc_mu::Float64 
            inc_sig,   #inc_sig::Float64 
            LogNormal(inc_mu, inc_sig),    #inc_dist::LogNormal{Float64}

            ## post-incubation (recovery) period 
            rec_min,    #rec_min::Float64
            rec_max,    #rec_max::Float64
            Uniform(rec_min, rec_max),    #rec_dist::Uniform{Float64}
                
            ## distribution of infectiousness 
            b_dispersion,    #b_dispersion::Float64
            beta / b_dispersion,    #b_scale::Float64
            Gamma(b_dispersion, beta/b_dispersion),    #b_dist::Gamma{Float64}

            ## scaling functions for piecewise 'viral load'
            V_max,    #Vmax::Float64 # scaling factor for growth and decline of viral load
            inc_plat_fac,    #inc_plat_fac::Float64 # proportion of incubation period in plateau viral load
            rec_plat_fac,    #rec_plat_fac::Float64 # proportion of recovery period in plateau viral load

            name, #i.e., "Delta variant"

            # ranges of test sensitivity function parameters 
            c_lims, # = [1.0, 5.11]
            b1_lims, # = [0.8, 2.31]
            b2_lims, # = [1.26, 3.47]
            b3_lims # = [1.05, 1.14]

    )

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

    # Initialise new infection
    infection(pathogen::Disease_T, t::Float64, t_inc::Float64, q_inc::Float64, t_rec::Float64, beta_max::Float64) = new(
    
    ## pathogen reference
        pathogen,
    ## name of pathogen 
        pathogen.name, # pathogen::String 
        t, #t_infected::Float64
    ## latent period
        pathogen.latent_period, # TODO: draw from distribution as with other intervals
    ## incubation period (time between exposure and symptom expression)
        t_inc, #rand(rng, pathogen.inc_dist),#t_inc::Float64
        q_inc,#cdf(pathogen.inc_dist, t_inc),#q_inc::Float64 NOTE: not sure if it will let me do this, need to check. *** Nope - need to pass as input to constructor. 
    ## recovery period (time between peak viral load and recovery)
        t_rec,#rand(rng, pathogen.rec_dist),#t_rec::Float64
    ## symptom expression (bool symptomatic or not) "Will they express symptoms?"
        rand(rng_infections) < (1.0 - pathogen.p_asymp),#symptomatic::Bool 
    ## symptom expression (bool expressing symptoms or not) "are they currently expressing symptoms?"
        false,#expressing_symptoms::Bool 
    ## force of infection (i.e., baseline infection transmission rate as a function of time)
        0.0,#beta_t::Float64 
    
    ### parameters determining force of infection
    ##maximum infectiousness (peaks at symptom onset)
    beta_max,#rand(rng, pathogen.b_dist),#beta_max::Float64 
    beta_max / pathogen.Vmax, #beta_min::Float64 #Def: beta_max/Vmax #[NOTE: again, not sure if Julia will let me initialise in this way (does it know the value of beta_max?)] ***

    compute_kinc(t_inc, pathogen),#k_inc::Float64 #rate of exponential increase of infectiousness during incubation *** (may need to initialise these after construction)
    compute_krec(t_rec, pathogen),#k_rec::Float64 #rate of exponential decrease of infectiousness during recovery ***

    ## parameters determining detection probability
    ## RAT detection probability, NOTE: decide whether this is an infection property, or something else... 
    ## for now, I'm going to include these as infection properties, but they could potentially be implemented 
    ## separately... 
    0.0,#test_sensitivity::Float64 (initially set to 0)
    pathogen.b1_lims[1] + rand(rng_infections) * (pathogen.b1_lims[2] - pathogen.b1_lims[1]),#b1::Float64
    pathogen.b2_lims[1] + (1.0 - q_inc) * (pathogen.b2_lims[2] - pathogen.b2_lims[1]),#b2::Float64 ***quantile matching with incubation period. 
    pathogen.b3_lims[1] + rand(rng_infections) * (pathogen.b3_lims[2] - pathogen.b3_lims[1]),#b3::Float64
    t_inc - (  q_inc * (pathogen.c_lims[2] - pathogen.c_lims[1]))#changepoint::Float64
    )

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
diseases = Dict{String, Disease_T}() #[disease_name] = Delta_variant 

#TODO: parameters are still Delta variant- adjust for omicron. 
disease_names = ["Omicron"]#["Delta Variant"]
n_diseases = length(disease_names) 

for d = 1:n_diseases

    disease_name = disease_names[d]

    # initialise disease (Delta variant, SARS-CoV-2) 
    # TODO: update for omicron
    R0 = 6.0
    beta = R0/3.94 
    #note: 3.94 is the scaling factor 
    # determined during the calibration of the hotel quarantine model
    #pre-dating this work 
    # TODO: verify calibration 
    p_asymp = 0.33
    latent_period = dt
    inc_mu = 1.62
    inc_sig = 0.418
    rec_min = 5.0
    rec_max = 10.0
    b_dispersion = 0.15
    V_max = 7.0
    inc_plat_fac = 0.1
    rec_plat_fac = 0.0
    #test sensitivity function parameters: 
    c_lims = [1.0, 5.11]
    b1_lims = [0.8, 2.31]
    b2_lims = [1.26, 3.47]
    b3_lims = [1.05, 1.14]

    disease_i = disease_params(
        beta,
        p_asymp, 
        latent_period,
        inc_mu,
        inc_sig,
        rec_min,
        rec_max,
        b_dispersion,
        V_max,
        inc_plat_fac,
        rec_plat_fac,
        disease_name,
        c_lims,
        b1_lims,
        b2_lims,
        b3_lims )

    diseases[disease_name] = disease_i

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
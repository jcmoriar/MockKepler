module MockKepler
using DataFrames
using MKS
using CurveFit
using JSON
using StatsBase

export load_sims


"""
    kepid_to_combo_index(kepids, combo_stars, reps)

Determine the indices into combo array for each star.

The combo array holds mock observation data for N representative
stars. For each star, the representative is determined (the star
with similar detection characteristics) and then its location in
the combo array.

# Arguments
* `kepids::Array`: kepids of the sample
* `combo_stars::Array`: kepids of the representatives used in the
    combo data array. Note: this must be in the same order as the
    combos data array.
* `reps::Dict`: Mapping of kepid to representative kepid

# Returns
* Array of indices into combo data array for each kepid's
    representative kepid.

"""
function kepid_to_combo_index(kepids, combo_stars, reps)
    indices = Int16[]
    for kepid in kepids
        repid = reps[kepid]
        push!(indices, findfirst(combo_stars, repid))
    end
    return(indices)
end


"""
    sims_from_model(weights, num)

Select random sims with weights according to a model.


# Arguments
* `weights::Array`: The length 10 array giving the weighted
    contributions of each set of intial conditions to the
    model.
* `num::Integer`: The number of samples to take

# Returns
An array of random sim indices distributed according to the
weights.

"""
function sims_from_model(weights, num; simsperset=32)
    all_weights = [1-sum(weights); [weights[div(i,simsperset)+1]/simsperset for i=0:(simsperset*length(weights))-1]]
    return(sample(0:length(weights)*simsperset, WeightVec(all_weights), num))
end

"""
    sims_from_model_subsample(weights, num)

Select random sims with weights according to a model from a resampled
set of the simulations.


# Arguments
* `weights::Array`: The length 10 array giving the weighted
    contributions of each set of intial conditions to the
    model.
* `num::Integer`: The number of samples to take

# Returns
An array of random sim indices distributed according to the
weights.

"""
function sims_from_model_subsample(weights, num; simsperset=32)
    all_weights = [1-sum(weights); [weights[div(i,simsperset)+1]/simsperset for i=0:(simsperset*length(weights))-1]]
    return(sample([0; resampled_sim_indices(length(weights), simsperset)], WeightVec(all_weights), num))
end


"""
    mock_campaign(starnames, stardata, simdata, prob_detect, sims)

    # Arguments
    * `starnames::Array`: kepids of each star
    * `stardata::Array`: Some stellar data (created from function star_array)
    * `simdata::Dict`: Simulation data
    * `prob_detect`: Function to calculate detection probability
    * `sims::Array`: The simulation to use for each star

    # Returns
    * `hosts::Array`: Kepids of the hosting stars in the same order
        as the observed planets.
    * `observed_planets::Array`: Orbital and physical parameters of
        the observed planets. [6, num_planets] where the first index
        corresponds to [ a, e, i, mass, radius, transit duration] in that order.

"""
function mock_campaign2(starnames, stardata, simdata, prob_detect, sims)
    observed_planets = Array(Float64, 6, round(Int, length(sims)/4))
    hosts = Array(Int64, round(Int, length(sims)/4))
    num_planets = 0
    for i=1:length(sims)
        # println(stardata[:,i])
        # println(simdata[sims[i]]-1)
        if sims[i] !=0
            durations = observe_system(stardata[:,i], simdata[sims[i]-1], prob_detect)
            for j=1:length(durations)
                if durations[j] != 0
                    num_planets += 1
                    hosts[num_planets] = starnames[i]

                    observed_planets[:,num_planets] = [transpose(simdata[sims[i]-1][j, [1,2,3,7,8]]);durations[j]]
                end
            end
        end
    end
    return(hosts[1:num_planets], observed_planets[:,1:num_planets])
end

"""
    mock_campaign(combos, stars, star_index, sim_index, planets)

    # Arguments
    * `combos::Array`: 4-D array containing the ocurrances of each
        observed planet combo in the mock observations.
    * `stars::Array`: The kepids of each star in the sample.
    * `star_index::Array`: The indices into combos for each star
        in the sample.
    * `sim_index::Array` The index of the simulation that will be
        used for each star. Must be same length as `star_index`
    * `planets::Array`: The orbital parameters of simulated planets

    # Returns
    * `hosts::Array`: Kepids of the hosting stars in the same order
        as the observed planets.
    * `observed_planets::Array`: Orbital and physical parameters of
        the observed planets. [5, num_planets] where the first index
        corresponds to [ a, e, i, mass, radius] in that order.
"""
function mock_campaign(combos::Array{Int16,4}, stars::Array{Int64,1}, star_index::Array{Int16,1}, sim_index::Array{Int64,1}, planets::Array{Float64,3})
    observed_planets = Array(Float64, 5, round(Int, length(star_index)/4))
    hosts = Array(Int64, round(Int, length(star_index)/4))
    num_planets = 0
    max_combos = size(combos)[2]
    for i=1:length(star_index)
        if sim_index[i] != 0
            stopat = rand(1:10000)
            for j=1:max_combos
                if combos[1,j,sim_index[i], star_index[i]] >= stopat
                    for k=2:7
                        if combos[k,j,sim_index[i], star_index[i]] != 0
                            num_planets += 1
                            hosts[num_planets] = stars[i]
                            observed_planets[:, num_planets] = planets[[1,2,3,7,8], combos[k,j,sim_index[i], star_index[i]], sim_index[i]]
                        end
                    end
                    break
                end
            end
        end
    end
    return(hosts[1:num_planets], observed_planets[:,1:num_planets])
end

""" Calculate the hill separations between all adjacent planets. """
function separations(kepids, planets)
    lastid = 0
    lastsma = 0
    lastmass = 0
    seps = []
    for i=1:length(kepids)
        sma = planets[1,i]
        mass = planets[4,i]
        if kepids[i] == lastid
            rhill = (lastsma+sma)/2*((lastmass+mass)/3)^(1./3)
            push!(seps, (sma-lastsma)/rhill)
        end
        lastid = kepids[i]
        lastsma = sma
        lastmass = mass
    end
    return(seps)
end

""" Calculate the hill separations between all adjacent planets. """
function separations(kois)
    kepids = kois[:kepid]
    planets = kois[[:koi_sma, :koi_prad]]
    lastid = 0
    lastsma = 0
    lastmass = 0
    seps = []
    for i=1:length(kepids)
        sma = planets[1][i]
        mass = planets[2][i].^2.06 * MEARTH/MSUN
        if kepids[i] == lastid
            rhill = (lastsma+sma)/2*((lastmass+mass)/3)^(1./3)
            push!(seps, (sma-lastsma)/rhill)
        end
        lastid = kepids[i]
        lastsma = sma
        lastmass = mass
    end
    return(seps)
end


""" Calculate the period ratios between all adjacent planets. """
function periodratios(kepids, planets)
    lastid = 0
    lastperiod = 0
    seps = []
    for i=1:length(kepids)
        period = planets[1,i]^1.5 * 365.25
        if kepids[i] == lastid
            push!(seps, period/lastperiod)
        end
        lastid = kepids[i]
        lastperiod = period
    end
    return(seps)
end

""" Calculate the period ratios between all adjacent planets. """
function periodratios(kois)
    kepids = kois[:kepid]
    planets = kois[[:koi_period, :koi_prad]]
    lastid = 0
    lastperiod = 0
    seps = []
    for i=1:length(kepids)
        period = planets[1][i]
        if kepids[i] == lastid
            push!(seps, period/lastperiod)
        end
        lastid = kepids[i]
        lastperiod = period
    end
    return(seps)
end

"""
    SimulatedCampaign

    Holds all the data needed for simulating the Kepler campaign.

    # Arguments
    * `combos::Array`: The observed combo data.
    * `stars::Array`: The Kepids of all stars in the sample
    * `star_index::Array`: The indices of the representative star for each star
        in the sample
    * `planets::Array`: The orbital description of each simulated planet
"""
type SimulatedCampaign
    combos::Array{Int16,4}
    stars::Array{Int64,1}
    star_index::Array{Int16,1}
    planets::Array{Float64,3}
end

"""
    mock_campaign(camp, model)

    # Arguments
    * `camp::SimulatedCampaign`: The specifications of the campaign
    * `model::Array`: The model description. (weights of each
        of the 10 possible components)

    # Returns
    * `hosts::Array`: Kepids of the hosting stars in the same order
        as the observed planets.
    * `observed_planets::Array`: Orbital and physical parameters of
        the observed planets. [5, num_planets] where the first index
        corresponds to [ a, e, i, mass, radius] in that order.
"""
function mock_campaign(camp::SimulatedCampaign, model)
    sim_indices = sims_from_model(model, length(camp.star_index))
    return(mock_campaign(camp.combos, camp.stars, camp.star_index, sim_indices, camp.planets))
end





"""
    count_multiplicity(hosts)

Count the number of ocurrances of each host and determine multiplicity dist.

"""
function count_multiplicity(hosts)
    current = hosts[1]
    num = 1
    mult = zeros(Int64, 7)
    for i=2:length(hosts)
        if hosts[i] == current
            num += 1
        else
            mult[num] += 1
            num = 1
            current = hosts[i]
        end
    end
    return(mult)
end

"""
    load_combos(filename)

Reads in mock observation combo data.

# Arguments
* `filename::String`: JSON formatted file where observed planet combos
    data is stored.

# Returns
* `combos::Array`: 4-D array of combo data.
    - Combos[1,combo_index, sim_index, star_index] = Number of times this
    particular combo was observed + counts of previous combos (i.e. the
    cummulative sum of counts).
    - Combos[2:end, combo_index, sim_index, star_index] = The indices of
    the planets observed. (i.e. defines the combo) Zeros are array filler
    which do not correspond to a planet.
* `stars::Array`: The kepids of the stars in the combo array in their
    index order.

"""
function load_combos(filename)
    combos = JSON.parsefile(filename)
    max_combos = 0
    for (i,key) in enumerate(keys(combos))
        for j=0:length(combos[key])-1
            len = length(combos[key][string(j)])
            max_combos = max(max_combos, len)
        end
    end
    stars = []
    combo_array = zeros(Int16, 8, max_combos, length(combos[collect(keys(combos))[1]]), length(combos))
    for (i, key) in enumerate(keys(combos))
        push!(stars, parse(Int,key))
        # println(i)
        for j=0:length(combos[key])-1
            combo_keys = sort(collect(keys(combos[key][string(j)])))[1:end-1]
            if combo_keys[1] == ""
                splice!(combo_keys,1)
                combos[key][string(j)]["none"] += 1
            end
            unshift!(combo_keys, "none")
            count = 0
            for (k, combo) in enumerate(combo_keys)
                count += combos[key][string(j)][combo]
                combo_array[1, k, j+1, i] = count
                if combo == "none"
                    continue
                end
                indices = map(parse, split(combo, ","))
                for (l, ind) in enumerate(indices)
                    if l < 8
                        combo_array[l+1, k, j+1, i] = ind+1
                    else
                        continue
                    end
                end
            end
        end
    end

    return(combo_array, stars)
end


"""
    resampled_sim_indices(num_sets, sets_per_sim)

Create array of sim indices that are resampled within each sim setup

"""
function resampled_sim_indices(num_sets, sims_per_set)
    sims = Int64[]
    for i=1:num_sets
        sims = [sims; rand(1:sims_per_set, sims_per_set)+(i-1)*sims_per_set]
    end
    return(sims)
end


"""
    simulated_multiplicity_probs(combos)

Create array that holds prob of each multiplicity for each star and
    set of initial conditions. Optionally use a fraction or resampled set of simulations.
    arr[mult,sim_conditions,star]

"""
function simulated_multiplicity_probs(combos; model_ind=[], simsperset=32)
    if length(model_ind) == 0
        model_ind = collect(1:size(combos)[3])
    end
    undersamp = size(combos)[3]/length(model_ind)
    num_obs = maximum(combos[1,:,1,1])
    probs = ones(Float64, 8, convert(Int64,size(combos)[3]/simsperset), convert(Int64, size(combos)[4])) #minimum of one count for nonzero Likelihood
    for i=1:size(combos)[4]
        for j=1:length(model_ind)
            for k=1:size(combos)[2]
                if combos[1,k,model_ind[j],i] == 0
                    break
                end
                mult = sum(combos[2:end,k,model_ind[j],i] .!= 0)
                try
                    probs[mult+1, ceil(Int64, model_ind[j]/simsperset), i] += combos[1,k,model_ind[j],i]-combos[1,k-1,model_ind[j],i]
                catch error
                    probs[mult+1, ceil(Int64, model_ind[j]/simsperset), i] += combos[1,k,model_ind[j],i]
                end
            end
        end
    end
    return(probs/simsperset/num_obs*undersamp)
end

"""
    simulated_multiplicity_probs(combos)

Create array that holds prob of each multiplicity for each star and
    set of initial conditions for some number of resamplings of the
    simulations.
    arr[mult,sim_conditions,star,resample_set]

"""
function simulated_multiplicity_probs(combos, num_resamp::Int; simsperset=32)
    num_sets = convert(Int64, size(combos)[3]/simsperset)
    probs = ones(Float64, 8, num_sets, size(combos)[4], num_resamp)
    for i=1:num_resamp
        probs[:,:,:,i] = simulated_multiplicity_probs(combos, model_ind=resampled_sim_indices(num_sets, simsperset), simsperset=simsperset)
    end
    return(probs)
end

"""
    simulated_multiplicity_probs_bysim(combos)

Create array that holds prob of each multiplicity for each star and
    each simulation. arr[mult,sim_conditions,star]

"""
function simulated_multiplicity_probs_bysim(combos; simsperset=32)
    num_obs = maximum(combos[1,:,1,1])
    probs = ones(Float64, 8, convert(Int64,size(combos)[3]), convert(Int64, size(combos)[4])) #minimum of one count for nonzero Likelihood
    for i=1:size(combos)[4]
        for j=1:size(combos)[3]
            for k=1:size(combos)[2]
                if combos[1,k,j,i] == 0
                    break
                end
                mult = sum(combos[2:end,k,j,i] .!= 0)
                try
                    probs[mult+1, j, i] += combos[1,k,j,i]-combos[1,k-1,j,i]
                catch error
                    probs[mult+1, j, i] += combos[1,k,j,i]
                end
            end
        end
    end
    return(probs/num_obs)
end


"""
    observed_multiplicity(kepids, obs)

Create a 1-D array of observed multiplicity for each star.

#Arguments
* `kepids::Array`: Array of kepids for each star in catalog
* `obs::Array`: Dataframe containing koi data
"""
function observed_multiplicity(kepids, obs)
    mult = zeros(Int64, length(kepids))
    grouped = by(obs, :kepid, nrow)
    for i=1:size(grouped)[1]
        mult[findfirst(kepids, grouped[i,1])] = grouped[i,2]
    end
    return(mult)
end


"""
    multiplicity_lnl(model, star_ind, mult, mult_probs)

Calculate the log likelihood of the model.

# Arguments
* `model::Array`: Array giving the weights for each component
    of model.
* `star_ind::Array`: The indices into mult_probs for each star
    in the sample.
* `mult::Array`: The observed multiplicities. Same length as star_ind
* `mult_probs::Array`: The probability of each multiplicity for each
    star and set of initial conditions.

# Returns
Log likelihood of the model

"""
function multiplicity_lnl(model, star_ind, mult, mult_probs::Array{Float64,3})
    noplanets = [1.,0,0,0,0,0,0,0]
    total = sum(model)
    if (total > 1.00000001) | (any(model .< 0))
        return -Inf
    end
    lnl = 0.
    for i=1:length(mult)
        l=0.
        for j=1:length(model)
            l += mult_probs[mult[i]+1,j,star_ind[i]] * model[j]
        end
        l += noplanets[mult[i]+1]*(1-total)
        lnl += log(l)
    end
    return(lnl)
end

"""
    multiplicity_lnl(model, star_ind, mult, mult_probs)

Calculate the log likelihood of the model for a resampled set of sims.

# Arguments
* `model::Array`: Array giving the weights for each component
    of model.
* `star_ind::Array`: The indices into mult_probs for each star
    in the sample.
* `mult::Array`: The observed multiplicities. Same length as star_ind
* `mult_probs::Array`: The probability of each multiplicity for each
    star and set of initial conditions for some number of resamplings
    of the simulations.

# Returns
Log likelihood of the model

"""
function multiplicity_lnl(model, star_ind, mult, mult_probs::Array{Float64,4})
    noplanets = [1.,0,0,0,0,0,0,0]
    total = sum(model)
    if (total > 1.00000001) | (any(model .< 0))
        return -Inf
    end
    lnl = 0.
    for h=1:size(mult_probs)[4]
        for i=1:length(mult)
            l=0.
            for j=1:length(model)
                l += mult_probs[mult[i]+1,j,star_ind[i],h] * model[j]
            end
            l += noplanets[mult[i]+1]*(1-total)
            lnl += log(l)
        end
    end
    return(lnl/10.)
end

"""
    multiplicity_lnl_subsample(model, star_ind, mult, mult_probs)

Calculate the log likelihood of the model.

# Arguments
* `model::Array`: Array giving the weights for each component
    of model.
* `star_ind::Array`: The indices into mult_probs for each star
    in the sample.
* `mult::Array`: The observed multiplicities. Same length as star_ind
* `mult_probs::Array`: The probability of each multiplicity for each
    star and set of initial conditions.
* `simsperset::Int`: The number of simulations for each set of initial conditions

# Returns
Log likelihood of the model

"""
function multiplicity_lnl_resample(model, star_ind, mult, mult_probs_bysim; simsperset=32)
    noplanets = [1.,0,0,0,0,0,0,0]
    total = sum(model)
    if (total > 1.00000001) | (any(model .< 0))
        return -Inf
    end
    lnl = 0.
    r = rand(1:simsperset, simsperset)
    for i=1:length(mult)
        l=0.
        for j=0:length(model)-1
            for k=1:simsperset
                l += mult_probs_bysim[mult[i]+1,j*simsperset+r[k],star_ind[i]] * model[j+1] / simsperset
            end
        end
        l += noplanets[mult[i]+1]*(1-total)
        lnl += log(l)
    end
    return(lnl)
end


"""
    load_sims(filename)

Read planetary system orbital data.

Data is returned as a dictionary with each key-value pair corresponding to
the simID and n x 8 array where n is the number of planets and the columns
are [a, e, i, O, o, M, mass, radius]

"""
function load_sims(filename)
    sims = readtable(filename)
    data = Dict()
    for simID in unique(sims[:simID])
        sim = sims[sims[:simID] .== simID, [:a, :e, :i, :O, :o, :M, :mass, :radius]]
        data[simID] = ones(size(sim))
        for i=1:8
            data[simID][:,i] = convert(Array, sim[:,i])
        end
    end
    return(data)
end

"""
    sims_as_array(sims)

Convert the sims dictionary to a 3-D array for faster calculations.

"""
function sims_as_array(sims)
    data = zeros(Float64, 8, 20, length(sims))
    for key in sort(collect(keys(sims)))
        sim = sims[key]
        data[:,1:size(sim)[1], key+1] = transpose(sim)
    end
    return(data)
end

""" Read in kepler representatives.

Stars are grouped by the detection efficiency characteristics. Each
group is represented by a single member in the group. This allows for
a reduced number of calculations on the representatives only.

"""
function load_representatives(filename)
    f = open(filename)
    lines = readlines(f)
    stars = lines[1]
    reps = lines[2]
    stars = split(stars)
    reps = split(reps)
    stars = map(x->parse(Int64, x), stars)
    reps = map(x->parse(Int64, x), reps)
    return(Dict(zip(stars,reps)))
end


""" Read in the Kepler stellar data and return it as a DataFrame. """
function load_stars(;filename="/Users/johncmoriarty/AeroFS/Software/python/simulatedtransitobs/q1_q17_dr24_stellar.csv")
    readtable(filename)
end


""" Read in the KOI data and return as a DataFrame. """
function load_kois(;filename="/Users/johncmoriarty/AeroFS/Software/python/simulatedtransitobs/q1_q17_dr24_koi.csv")
    readtable(filename)
end

""" Read in the planet injection and recovery data and return as a DataFrame. """
function load_injections(;filename="/Users/johncmoriarty/AeroFS/Software/python/simulatedtransitobs/DR24-Pipeline-Detection-Efficiency-Table.csv")
    readtable(filename)
end

"""
    observe_system!(observed_planets, star, planets, detection_probability_function)

Determine which planets will be observed to transit from a random viewing direction around the given star.

# Arguments
* `observed_planets::DataFrame`: DataFrame to append the observed planets to. Should have columns:
    [kepid, koi_sma, koi_period, koi_duration, koi_prad, koi_eccen, mass]
* `star::DataFrame`: DataFrame row with stellar information
* `planets::Array`: Orbital data of planetary system to observe
* `detection_probability_function`: f(star_data, rp, sma, e) calculates the probability of detecting
    a given transiting plaent.

"""
function observe_system!(observed_planets, star, planets, detection_probability_function)
    println("finish this function. It needs to append observed planets to DataFrame.")
    num_tranets = 0
    sign = rand()-0.5
    sign = sign/abs(sign)
    theta = acos(2*rand()-1)*sign
    phi = 2pi*rand()
    reference = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    for i=1:size(planets)[1]
        a = planets[i,1]
        e = planets[i,2]

        rr = copy(reference)
        rotatez!(rr, -planets[i,4])
        rotatex!(rr, -planets[i,3])
        rotatez!(rr, -planets[i,5])

        nu = atan2(rr[2], rr[1])
        orbital_distance_at_nu = a*(1-e*e)/(1+e*cos(nu))
        impact = abs(orbital_distance_at_nu*sqrt(1-(rr[1]*rr[1]+rr[2]*rr[2]))) * AU / RSUN / star[1]
        if (impact < 1.0) & (detection_probability_function(star, planets[i,8], a, e) > rand())
            num_tranets += 1
        end

    end
    num_tranets
end


"""
    observe_system(star, planets, detection_probability_function)

Determine which planets will be observed to transit from a random viewing direction around the given star.

# Arguments
* `star::Array`: Array of stellar data
* `planets::Array`: Orbital data of planetary system to observe
* `detection_probability_function`: f(star_data, rp, sma, e) calculates the probability of detecting
    a given transiting plaent.

# Returns
* `transit_durations::Array`: Normalized transit durations of observed planets (or zero if not transiting).

"""
function observe_system(star, planets, detection_probability_function)
    sign = rand()-0.5
    sign = sign/abs(sign)
    theta = acos(2*rand()-1)*sign
    phi = 2pi*rand()
    reference = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    durations = zeros(Float64, size(planets)[1])
    for i=1:size(planets)[1]
        a = planets[i,1]
        e = planets[i,2]

        rr = copy(reference)
        rotatez!(rr, -planets[i,4])
        rotatex!(rr, -planets[i,3])
        rotatez!(rr, -planets[i,5])

        nu = atan2(rr[2], rr[1])
        orbital_distance_at_nu = a*(1-e*e)/(1+e*cos(nu))
        impact = abs(orbital_distance_at_nu*sqrt(1-(rr[1]*rr[1]+rr[2]*rr[2]))) * 215.0941768511862 / star[1]
        if (impact < 1.0) & (detection_probability_function(star, planets[i,8], a, e) > rand())
            v_orb = sqrt(887373835.2117007*(2./orbital_distance_at_nu-1./a))
            duration = 2*star[1]*RSUN*sqrt((1+planets[i,8]*REARTH/star[1]/RSUN)^2-impact*impact)/v_orb/(sqrt(a*a*a/star[2])*365.25)
            durations[i] = duration
        end

    end
    return(durations)
end


# """
#     batch_mock_function()

# Create a function that when called returns a batch of simulated obs.


# """

"""
    batch_mock_observation(planets; n=100000, rstar=1., mstar=1.)

Determine which planets will transit from a random viewing direction.

# Arguments
* `planets::Array`: Orbital data of planetary system to observe
* `n::Int`: The number of mock observations
* `rstar::Float`: Radius of the star is solar radii
* `mstar::Float`: Mass of the star is solar masses

# Returns
* `durations::Array`: nplanets rows by n columns array containing 0 if not
    observed, otherwise the normalized transit duration

"""
function batch_mock_observation(planets; n=100000, rstar=1., mstar=1.)
    durations = zeros(size(planets)[1], n)
    for j=1:n
        sign = rand()-0.5
        sign = sign/abs(sign)
        theta = acos(2*rand()-1)*sign
        phi = 2pi*rand()
        reference = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
        for i=1:size(planets)[1]
            a = planets[i,1]
            e = planets[i,2]

            rr = copy(reference)
            rotatez!(rr, -planets[i,4])
            rotatex!(rr, -planets[i,3])
            rotatez!(rr, -planets[i,5])

            nu = atan2(rr[2], rr[1])
            orbital_distance_at_nu = a*(1-e*e)/(1+e*cos(nu))
            impact = abs(orbital_distance_at_nu*sqrt(1-(rr[1]*rr[1]+rr[2]*rr[2]))) * 215.0941768511862 / rstar
            if impact < 1.0
                v_orb = sqrt(887373835.2117007*(2./orbital_distance_at_nu-1./a))
                durations[i,j] = 2*rstar*RSUN*sqrt((1+planets[i,8]*REARTH/rstar/RSUN)^2-impact*impact)/v_orb/(sqrt(a*a*a/mstar)*365.25)
                # push!(durations, duration)
            end

        end
    end
    return(durations)
end


function multiplicity_duration_histograms(sims, stars, detection_probability_function)
    mult = zeros(size(stars)[2], length(keys(sims)), 7)
    for sim in sort(collect(keys(sims)))
        print("x")
        print("\r"*string(sim))
        obs = batch_mock_observation(sims[sim])
        rand_arr = rand(size(obs))
        prob_detect = zeros(size(sims[sim])[1])
        for k=1:size(stars)[2]
            for j=1:size(sims[sim])[1]
                prob_detect[j] = detection_probability_function(stars[:,k], sims[sim][j,8], sims[sim][j,1], sims[sim][j,2])
            end

            detected = (prob_detect .> rand_arr) .* (obs .> 0)
            num_detected = reshape(sum(detected,1), size(detected)[2])
            _, m = hist(num_detected, [0,1,2,3,4,5,6,7]-0.5)
            mult[k, sim+1, :] = m
        end
    end
    return(mult)
end


""" Rotate a length 3 vector about the z axis by the given angle"""
function rotatez!(vec, angle)
    s = sin(angle)
    c = cos(angle)
    temp = vec[1]
    vec[1] = c*vec[1]-s*vec[2]
    vec[2] = s*temp + c*vec[2]
end

""" Rotate a length 3 vector about the x axis by the given angle"""
function rotatex!(vec, angle)
    s = sin(angle)
    c = cos(angle)
    temp = vec[2]
    vec[2] = c*vec[2]-s*vec[3]
    vec[3] = s*temp + c*vec[3]
end

"""
    subsample(;trange=[4200, 6100])

Load the stellar, koi, and injection data and return corresponding subsamples with selection criteria applied.

# Arguments
`trange::Array[2]`: lower and upper stellar temperature limits to select on

# Return
`(stars, kois, inj)`: Correlated subsample of stellar, koi, and injection data stored in DataFrames.

"""
function subsample(;trange=[4200, 6100], min_period=0)
    stars = load_stars()
    kois = load_kois()
    inj = load_injections()

    # apply stellar selection criteria
    stars = stars[stars[:logg] .>= 4,:]
    stars = stars[trange[1] .<= stars[:teff] .< trange[2],:]
    stars = stars[~isna(stars[:dataspan]),:]
    stars = stars[(stars[:dutycycle] .* stars[:dataspan]) .> 2*365.25,:]
    stars = stars[stars[:dutycycle] .> 0.33,:]
    stars = stars[~isna(stars[:mass]),:]

    # correlate koi sample with stellar sample
    kois = join(kois, stars[[:kepid, :teff]], on=:kepid, kind=:inner)
    kois = kois[kois[:koi_pdisposition] .== "CANDIDATE",:]
    kois = kois[~isna(kois[:koi_period]),:]
    kois = kois[kois[:koi_period] .> min_period,:]
    kois = kois[kois[:koi_max_mult_ev] .>= 15.0,:]
    sort!(kois, cols=[:kepid, :koi_sma])

    # correlate injection sample with stellar sample
    inj = join(inj, stars[[:kepid, :teff]], on=:kepid, kind=:inner)

    return(stars, kois, inj)

end

"""
For faster access of stellar parameters. [radius, mass, dataspan, dutycycle, cdpp values]
"""
function star_array(stars)
    cdpp_colnames = filter(function swrrmscdpp(x) startswith(x, "rrmscdpp") end, map(string, stars.colindex.names))
    cdpp_values = [float(replace(replace(name, "rrmscdpp", ""), "p", ".")) for name in cdpp_colnames]
    cdpp_cols = map(symbol, cdpp_colnames)

    arr = transpose([stars[:radius].data stars[:mass].data stars[:dataspan].data stars[:dutycycle].data])
    for col in cdpp_cols
        arr = [arr; transpose(stars[col].data)]
    end
    return arr

end

# This accesses stellar data as a DataFrame which is very slow. The new version uses arrays and is reccomended.
function detect_prob_old(stars, inj)
    cdpp_colnames = filter(function swrrmscdpp(x) startswith(x, "rrmscdpp") end, map(string, stars.colindex.names))
    cdpp_values = [float(replace(replace(name, "rrmscdpp", ""), "p", ".")) for name in cdpp_colnames]
    cdpp_cols = map(symbol, cdpp_colnames)

    # determine the dependance of the detection efficiency on MES
    bin_edges = linspace(15,80,20)
    bins = (bin_edges[2:end]+bin_edges[1:end-1])/2
    total_inj = hist(inj[:expect_mes], bin_edges)[2]
    recovered_inj = hist(inj[(inj[:recovered] .== 1) & (inj[:meas_mes].>15),:][:expect_mes], bin_edges)[2]
    recovered_fraction = float(recovered_inj)./total_inj
    completeness_fit = curve_fit(Poly, collect(bins), recovered_fraction, 1)

    function detection_probability_function(star, rp, sma, e)

        # Calculate period and transit duration
        period = sqrt(sma*sma*sma/star[:mass][1])*365.25
        tau = 0.25 * period * sqrt(1 - e*e) / (215.0941768511862*sma/star[:radius][1]) * 24

        # Get CDPP for duration (change to interpolation depending on speed)
        sigma = star[cdpp_cols[indmin(abs(cdpp_values - tau))]][1]

        # Radius ratio and signal to noise
        k = rp / star[:radius][1] * 0.009170524802300503
        snr = 0.84 * k*k * (1.0874 + 1.0187*k) * 1e6 / sigma

        # scale by the number of transits
        num_transits = star[:dataspan][1] * star[:dutycycle][1] / period
        mes = snr * sqrt(num_transits)

        # detection probability assumes that kois with mes < 15 are tossed out
        # Done this way because of extra dependence of pdet on period? in dr24
        pdet = min(1.,completeness_fit(mes)) * (mes > 15)

        # Calc effect due to window function
        M = star[:dataspan][1] / period
        f = star[:dutycycle][1]
        omf = 1.0 - f
        pw = 1 - omf^M - M*f*omf^(M-1) - 0.5*M*(M-1)*f*f*omf^(M-2)
        msk = (pw >= 0.0) * (M >= 2.0)
        pwin = pw * msk

        return(pwin*pdet)
    end

end

""" Get the cdpp value closest to the transit duration (tau). """
function closest_cdpp(star, cdpp, tau)
    min = 1e10
    best = star[5]
    if tau > cdpp[end]
        return star[end]
    end
    for i=1:length(cdpp)
        if abs(tau-cdpp[i]) < min
            min = abs(tau-cdpp[i])
            best = star[i]+4
        else
            return best
        end
    end
    return best
end

function detect_prob(stars, inj)
    cdpp_colnames = filter(function swrrmscdpp(x) startswith(x, "rrmscdpp") end, map(string, stars.colindex.names))
    cdpp_values = [float(replace(replace(name, "rrmscdpp", ""), "p", ".")) for name in cdpp_colnames]
    cdpp_cols = map(symbol, cdpp_colnames)

    # determine the dependance of the detection efficiency on MES
    bin_edges = linspace(15,80,20)
    bins = (bin_edges[2:end]+bin_edges[1:end-1])/2
    total_inj = hist(inj[:expect_mes], bin_edges)[2]
    recovered_inj = hist(inj[(inj[:recovered] .== 1) & (inj[:meas_mes].>15),:][:expect_mes], bin_edges)[2]
    recovered_fraction = float(recovered_inj)./total_inj
    completeness_fit = curve_fit(Poly, collect(bins), recovered_fraction, 1)

    function detection_probability_function(star, rp, sma, e)

        # Calculate period and transit duration
        period = sqrt(sma*sma*sma/star[2])*365.25
        tau = 0.25 * period * sqrt(1 - e*e) / (215.0941768511862*sma/star[1]) * 24

        # Get CDPP for duration (change to interpolation depending on speed)
        # sigma = star[indmin(abs(cdpp_values - tau))+4]
        sigma = closest_cdpp(star, cdpp_values, tau)

        # Radius ratio and signal to noise
        k = rp / star[1] * 0.009170524802300503
        snr = 0.84 * k*k * (1.0874 + 1.0187*k) * 1e6 / sigma

        # scale by the number of transits
        num_transits = star[3] * star[4] / period
        mes = snr * sqrt(num_transits)

        # detection probability assumes that kois with mes < 15 are tossed out
        # Done this way because of extra dependence of pdet on period? in dr24
        pdet = min(1.,completeness_fit(mes)) * (mes > 15)

        # Calc effect due to window function
        M = star[3] / period
        f = star[4]
        omf = 1.0 - f
        pw = 1 - omf^M - M*f*omf^(M-1) - 0.5*M*(M-1)*f*f*omf^(M-2)
        msk = (pw >= 0.0) * (M >= 2.0)
        pwin = pw * msk

        return(pwin*pdet)
    end

end

end
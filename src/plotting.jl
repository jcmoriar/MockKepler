module plotting
import PyPlot
import MockKepler
using DataFrames
reload("MockKepler")
using MKS

function plot_period_dist(kois, mock_planets)
    bins = linspace(0,3,20)
    bins = 10.^bins
    PyPlot.plt[:hist](kois[:koi_period], bins=bins, alpha=0.5, label="KOIs")
    len = size(mock_planets)[2]
    sma = reshape(mock_planets[1,:],len)
    period = sqrt(sma.*sma.*sma)*365.25
    PyPlot.plt[:hist](period, bins=bins, alpha=0.5, label="Simulated")
    PyPlot.xscale("log")
    PyPlot.xlabel("Period [days]")
    PyPlot.ylabel("# of planets")
    PyPlot.legend()

end

function plot_radius_dist(kois, mock_planets)
    bins = linspace(log10(.5),log10(5),20)
    bins = 10.^bins
    PyPlot.plt[:hist](dropna(kois[:koi_prad]), bins=bins, alpha=0.5, label="KOIs")
    len = size(mock_planets)[2]
    # sma = reshape(mock_planets[1,:],len)
    # period = sqrt(sma.*sma.*sma)*365.25
    radii = reshape(mock_planets[5,:],len)
    PyPlot.plt[:hist](radii, bins=bins, alpha=0.5, label="Simulated")
    PyPlot.xscale("log")
    PyPlot.xlabel("Radius [Earth Radii]")
    PyPlot.ylabel("# of planets")
    PyPlot.legend()
end

function plot_period_radius(kois, mock_planets)
    PyPlot.scatter(kois[:koi_period].data, kois[:koi_prad].data, label="KOIs")
    len = size(mock_planets)[2]
    sma = reshape(mock_planets[1,:],len)
    radii = reshape(mock_planets[5,:],len)
    noise_p = rand(size(radii))*0.4+.8
    noise_r = rand(size(radii))*0.4+.8
    PyPlot.scatter(sqrt(sma.*sma.*sma)*365.25.*noise_p, radii.*noise_r,color="red", alpha=0.5, label="Simulated")
    PyPlot.xscale("log")
    PyPlot.yscale("log")
    PyPlot.ylim([.5,5])
    PyPlot.xlim([1e0, 1e3])
    PyPlot.xlabel("Period [days]")
    PyPlot.ylabel("Radius [\$R_{\\oplus}\$]")
    PyPlot.legend()

end


function plot_sma_mass(planets)
    len = size(planets)[2]
    sma = reshape(planets[1,:],len)
    mass = reshape(planets[4,:],len)*MSUN/MEARTH
    PyPlot.scatter(sma, mass,color="red")
    PyPlot.plot([0.05,0.05],[-10,50])
    PyPlot.plot([0.07,0.07],[-10,50])
    # PyPlot.xscale("log")
    # PyPlot.yscale("log")
    # PyPlot.ylim([1,5])
    PyPlot.xlim([0, 1])
end

function plot_chain(chain)
    skip = round(Int64, size(chain)[2]/500)
    for i=1:size(chain)[1]
        points = chain[i,1:skip:end][:]
        PyPlot.scatter(i*ones(size(points)) + (.8+.4*rand(size(points))), points, alpha=.2)
    end
    # PyPlot.xticks([2, 3, 4, 5], ["10\$M_{\\oplus}\$", "20\$M_{\\oplus}\$", "30\$M_{\\oplus}\$", "40\$M_{\\oplus}\$"])
    PyPlot.xlabel("Disk Component")
    PyPlot.ylabel("Contribution to Model")
    PyPlot.ylim([0,maximum(chain)])
    # PyPlot.savefig("../plots/MCMCchain.pdf")

end

function plot_model_multiplicity(model, camp, actual_mult)
    edges, actual = hist(actual_mult, collect(0:8)-.5)
    maxmult = zeros(Int64, 8)
    minmult = zeros(Int64, 8) + 999999999
    num_observed_systems = 0
    simsper = convert(Int32, size(camp.combos)[3]/length(model))
    for i=1:100

        # mockstars, mockplanets = MockKepler.mock_campaign(camp, model)
        sim_indices = MockKepler.sims_from_model(model, length(camp.star_index), simsperset=simsper)
        # sim_indices = rand(1:1, length(camp.star_index))
        # sim_indices = ones(Int64, length(camp.star_index))+1
        mockstars, mockplanets = MockKepler.mock_campaign(camp.combos, camp.stars, camp.star_index, sim_indices, camp.planets)
        mult = MockKepler.count_multiplicity(mockstars)
        mult = [length(camp.stars)-sum(mult);mult]
        num_observed_systems += sum(mult[2:end])
        for j=1:8
            if mult[j] > maxmult[j]
                maxmult[j] = mult[j]
            end
            if mult[j] < minmult[j]
                minmult[j] = mult[j]
            end
        end
    end
    println("Average number of observed systems: ",num_observed_systems/100., " out of: ", length(camp.stars))
    println()
    # println(edges)
    # println(actual)
    PyPlot.scatter(collect(0:7), actual)
    PyPlot.fill_between(collect(0:7), minmult+(minmult.==0).*.1, maxmult, alpha=.5, color=(0.5,0.5,0.5))
    PyPlot.yscale("log")
    PyPlot.ylim([.1,1.1e5])
    PyPlot.xlim([-0.1,6.1])
    PyPlot.xlabel("Number of planets")
    PyPlot.ylabel("Number of stars")
end

function plot_mass_sma_sims(sims)
    smas = sims[1,:,:][:]
    masses = sims[7,:,:][:]*MSUN/MEARTH
    PyPlot.scatter(smas, masses)
    PyPlot.xlim([0,2])

end


"""
    plot_all(chain, camp, kois; suffix="none")

    Create all the relevant plots:
        Multiplicity Distribution
        Radius Distribution
        Period Distribution
        Transit Duration Ratio Distribution

    # Arguments
    * `chain::Array`: 2-D array of MCMC chain of models to plot [num_model_compenents, num_models]
    * `camp::SimulatedCampaign`: The data necessary to simulate transit observations
    * `kois::DataFrame`: KOI data
    * `dir::String`:  The directory in which to save figures. If "none", the figures
        will not be saved.
    * `suffix::String`: (Optional) The suffix to append to each figure filename.

    # Not yet working
    * `star_data::Array`: (Optional) Include if transit duration ratio plots are desired. This
        will take much longer. Contains stellar data necessary to calculate probability of detection.
    * `sims::Dict`: (Optional) Also required for transit duration ratios
    * `prob_d`: (Optional) Also required for transit duration ratios. Function to calculate detection
        probability for a planet around a given star.

"""


function plot_all(chain, camp, kois; dir="/Users/johncmoriarty/AeroFS/Projects/KeplerDichotomy/ParameterizedInsideOut/plots/", suffix="", star_data=0, sims=0, prob_d=0)
    if dir[end] != "/"
        dir *= "/"
    end

    skip = round(Int64, size(chain)[2]/500)
    models = chain[:,1:skip:end]
    observed_multiplicity = MockKepler.observed_multiplicity(camp.stars, kois)

    edges_mult, actual_mult = hist(observed_multiplicity, collect(0:8)-.5)
    maxmult = zeros(Int64, 8)
    minmult = zeros(Int64, 8) + 999999999
    num_observed = []

    edges_r, actual_r = hist(dropna(kois[:koi_prad]), 10.^linspace(log10(.5),log10(10),15))
    maxr = zeros(Int64, length(edges_r)-1)
    minr = zeros(Int64, length(edges_r)-1) + 999999999
    rhist = Array(Int64, length(edges_r)-1, size(models)[2])

    edges_p, actual_p = hist(dropna(kois[:koi_period]), 10.^linspace(0,3,15))
    maxp = zeros(Int64, length(edges_p)-1)
    minp = zeros(Int64, length(edges_p)-1) + 999999999
    phist = Array(Int64, length(edges_p)-1, size(models)[2])

    edges_sep, actual_sep = hist(MockKepler.separations(kois), linspace(0,80,20))
    maxsep = zeros(Int64, length(edges_sep)-1)
    minsep = zeros(Int64, length(edges_sep)-1) + 999999999
    sephist = Array(Int64, length(edges_sep)-1, size(models)[2])

    edges_perratio, actual_perratio = hist(MockKepler.periodratios(kois), linspace(1,10,20))
    maxperratio = zeros(Int64, length(edges_perratio)-1)
    minperratio = zeros(Int64, length(edges_perratio)-1) + 999999999
    perratiohist = Array(Int64, length(edges_perratio)-1, size(models)[2])

    simsper = convert(Int32, size(camp.combos)[3]/size(models)[1])
    for i=1:size(models)[2]

        sim_indices = MockKepler.sims_from_model_subsample(models[:,i], length(camp.star_index), simsperset=simsper)
        sim_indices = MockKepler.sims_from_model(models[:,i], length(camp.star_index), simsperset=simsper)
        if sims == 0
            mockstars, mockplanets = MockKepler.mock_campaign(camp.combos, camp.stars, camp.star_index, sim_indices, camp.planets)
        else
            mockstars, mockplanets = MockKepler.mock_campaign2(camp.stars, star_data, sims, prob_d, sim_indices)
        end

        len = size(mockplanets)[2]


        # Multiplicity
        mult = MockKepler.count_multiplicity(mockstars)
        mult = [length(camp.stars)-sum(mult);mult]
        for j=1:8
            if mult[j] > maxmult[j]
                maxmult[j] = mult[j]
            end
            if mult[j] < minmult[j]
                minmult[j] = mult[j]
            end
        end

        # Radii
        radii = reshape(mockplanets[5,:],len)
        _, radii = hist(radii, edges_r)
        rhist[:,i] = radii
        for j=1:length(edges_r)-1
            if radii[j] > maxr[j]
                maxr[j] = radii[j]
            end
            if radii[j] < minr[j]
                minr[j] = radii[j]
            end
        end

        # Period
        periods = reshape(mockplanets[1,:],len).^1.5 .* 365.25
        _, periods = hist(periods, edges_p)
        phist[:,i] = periods
        for j=1:length(edges_p)-1
            if periods[j] > maxp[j]
                maxp[j] = periods[j]
            end
            if periods[j] < minp[j]
                minp[j] = periods[j]
            end
        end

        # Hill Separations
        seps = MockKepler.separations(mockstars, mockplanets)
        _, seps = hist(seps, edges_sep)
        sephist[:,i] = seps
        for j=1:length(edges_sep)-1
            if seps[j] > maxsep[j]
                maxsep[j] = seps[j]
            end
            if seps[j] < minsep[j]
                minsep[j] = seps[j]
            end
        end

        # Period Ratios
        perratios = MockKepler.periodratios(mockstars, mockplanets)
        _, perratios = hist(perratios, edges_perratio)
        perratiohist[:,i] = perratios
        for j=1:length(edges_perratio)-1
            if perratios[j] > maxperratio[j]
                maxperratio[j] = perratios[j]
            end
            if perratios[j] < minperratio[j]
                minperratio[j] = perratios[j]
            end
        end

    end

    PyPlot.clf()
    PyPlot.fill_between(collect(0:7), minmult+(minmult.==0).*.1, maxmult, color=color_simulated, label="Simulated")
    PyPlot.plot(collect(0:7), actual_mult, color=color_koi, label="KOIs", marker="o", ms=10, ls="None")
    PyPlot.yscale("log")
    PyPlot.ylim([.1,1.1e5])
    PyPlot.xlim([-0.1,6.1])
    PyPlot.xlabel("Number of planets")
    PyPlot.ylabel("Number of stars")
    PyPlot.legend()
    if dir != "none/"
        PyPlot.savefig(dir*"Multiplicity"*suffix*".pdf")
    end

    PyPlot.clf()
    bins_r = (edges_r[2:end]+edges_r[1:end-1])/2
    PyPlot.errorbar(bins_r, mean(rhist,2), yerr=std(rhist, 2), color=color_simulated, marker="o", ms=10, ls="None", label="Simulated")
    PyPlot.plt[:hist](dropna(kois[:koi_prad]), bins=edges_r, color=color_koi, label="KOIs")
    PyPlot.xscale("log")
    PyPlot.ylim([0,600])
    PyPlot.xlim([.4,10])
    PyPlot.ylabel("Number of planets")
    PyPlot.xlabel("Radius [\$R_{\\oplus}\$]")
    PyPlot.legend()
    if dir != "none/"
        PyPlot.savefig(dir*"Radius"*suffix*".pdf")
    end

    PyPlot.clf()
    bins_p = (edges_p[2:end]+edges_p[1:end-1])/2
    PyPlot.errorbar(bins_p, mean(phist,2), yerr=std(phist, 2), color=color_simulated, marker="o", ms=10, ls="None", label="Simulated")
    PyPlot.plt[:hist](dropna(kois[:koi_period]), bins=edges_p, color=color_koi, label="KOIs")
    PyPlot.xscale("log")
    PyPlot.ylim([0,600])
    PyPlot.xlim([1,1000])
    PyPlot.ylabel("Number of planets")
    PyPlot.xlabel("Orbital Period [Days]")
    PyPlot.legend()
    if dir != "none/"
        PyPlot.savefig(dir*"Period"*suffix*".pdf")
    end

    PyPlot.clf()
    bins_sep = (edges_sep[2:end]+edges_sep[1:end-1])/2
    PyPlot.errorbar(bins_sep, mean(sephist,2), yerr=std(sephist, 2), color=color_simulated, marker="o", ms=10, ls="None", label="Simulated")
    PyPlot.plt[:hist](MockKepler.separations(kois), bins=edges_sep, color=color_koi, label="KOIs")
    PyPlot.ylim([0,150])
    PyPlot.ylabel("Number of Planet Pairs")
    PyPlot.xlabel("Planet Separation [Hill Radii]")
    PyPlot.legend()
    if dir != "none/"
        PyPlot.savefig(dir*"Hill"*suffix*".pdf")
    end

    PyPlot.clf()
    bins_perratio = (edges_perratio[2:end]+edges_perratio[1:end-1])/2
    PyPlot.errorbar(bins_perratio, mean(perratiohist,2), yerr=std(perratiohist, 2), color=color_simulated, marker="o", ms=10, ls="None", label="Simulated")
    PyPlot.plt[:hist](MockKepler.periodratios(kois), bins=edges_perratio, color=color_koi, label="KOIs")
    PyPlot.ylim(ymin=0)
    PyPlot.ylabel("Number of Planet Pairs")
    PyPlot.xlabel("Period Ratio")
    PyPlot.legend()
    if dir != "none/"
        PyPlot.savefig(dir*"PeriodRatio"*suffix*".pdf")
    end

    skip = round(Int64, size(chain)[2]/500)


    PyPlot.clf()
    PyPlot.plt[:hist](sum(chain, 1)[:], bins=linspace(0,1,50), color=color_simulated)
    PyPlot.xlabel("Fraction of stars with planets")
    PyPlot.yticks([])
    if dir != "none/"
        PyPlot.savefig(dir*"OccuranceRate"*suffix*".pdf")
    end

    PyPlot.clf()
    for i=1:size(chain)[1]
        points = chain[i,1:skip:end][:]
        PyPlot.plot(i-1 + zeros(size(points)) + (.8+.4*rand(size(points))), points, alpha=.2, color=color_simulated, marker="o", ms=5, ls="None")
    end
    # PyPlot.xticks([1, 2, 3, 4], ["10\$M_{\\oplus}\$", "20\$M_{\\oplus}\$", "30\$M_{\\oplus}\$", "40\$M_{\\oplus}\$"])
    PyPlot.xticks(collect(1:size(chain)[1]))
    PyPlot.xlabel("Disk Component")
    PyPlot.ylabel("Contribution to Model")
    PyPlot.ylim([0,maximum(chain)])
    PyPlot.xlim([0,size(chain)[1]+1])
    if dir != "none/"
        PyPlot.savefig(dir*"ModelComponents"*suffix*".pdf")
    end

end

color_simulated = [0.492188, 0.492188, 0.984375]
color_koi = [0.578125, 0.746094, 0.507813]


end
module utils
using MKS

"""
    random_stream(;dist="uniform", lowerlim=0, upperlim=1)

Creates a function that when called will provide random numbers according
to the specifications.

"""
function random_stream(;dist="uniform", lowerlim=0, upperlim=1)
    if dist == "uniform"
        function f(n)
            return(rand(n)*(upperlim-lowerlim)+lowerlim)
        end
        return f
    end

end
"""
    basicMCMC(logL, x0, n_samples)

"""
function basicMCMC(logL, x0, n_samples; proposal_width=.1)
    accepted = 0
    npars = length(x0)
    chain = Array(Float64, npars, n_samples)
    lnl = Array(Float64, n_samples)
    chain[:,1] = x0
    chainlength = 1
    currentlnl = logL(x0)
    currentx = copy(x0)
    lnl[1] = currentlnl
    while chainlength < n_samples
        proposed_x = randn(npars)*proposal_width + currentx
        proposed_lnl = logL(proposed_x)
        a = exp(proposed_lnl-currentlnl)
        if a > rand()
            accepted +=1
            chainlength += 1
            chain[:,chainlength] = proposed_x
            lnl[chainlength] = proposed_lnl
            currentx = proposed_x
            currentlnl = proposed_lnl
        else
            chainlength += 1
            chain[:,chainlength] = currentx
            lnl[chainlength] = currentlnl
        end
    end
    println("Acceptance Rate: ",accepted/n_samples)
    return(chain, lnl)
end

"""
    integrate_disk(inner, outer, alpha, sigma0)

Determine mass within limits of powerlaw surface density disk.

# Arguments
* `inner`: Inner integration limit [AU]
* `outer`: Outer integration limit [AU]
* `alpha`: Power law index of disk (include - sign)
* `sigma0`: Surface density at 1 AU [g/cm^2]

# Returns
Mass of disk in Earth Masses
"""
function integrate_disk(inner, outer, alpha, sigma0)
    return(2*pi*AU^(-alpha)*(sigma0*10)/(alpha+2)*((outer*AU)^(alpha+2)-(inner*AU)^(alpha+2))/MEARTH)
end

"""
    determine_sigma0(inner, outer, alpha, mass)

Determine the surface density normalization of the specified disk.

# Arguments
* `inner`: Inner edge of disk [AU]
* `outer`: Outer edge of disk [AU]
* `alpha`: Power law index of disk (include - sign)
* `mass`: Mass of disk in Earth Masses

# Returns
The surface density normalization at 1 AU in g/cm^2.

"""
function determine_sigma0(inner, outer, alpha, mass)
    return(mass*MEARTH / 2. / pi / (AU^(-alpha)) * (alpha+2) / ((outer*AU)^(alpha+2) - (inner*AU)^(alpha + 2))/10)
end

end
# following Joe Fowler's drift_correct from mass

    # """Compute a drift correction that minimizes the spectral entropy of
    # <uncorrected> (a filtered pulse height vector) after correction.
    # The <indicator> vector should be of equal length and should have some 
    # linear correlation with the pulse gain. Generally it will be the
    # pretrigger mean of the pulses, but you can experiment with other choices.
    
    # The entropy will be computed on corrected values only in the range
    # [0, <limit>], so <limit> should be set to a characteristic large
    # value of <uncorrected>. If <limit> is None (the default), then
    # in will be compute as 25% larger than the 99%ile point of
    # <uncorrected>
    
    # The model is that the filtered pulse height PH should be scaled by
    # (1 + a*PTM) where a is an arbitrary parameter computed here, and 
    # PTM is the difference between each record's pretrigger mean and the
    # median value of all pretrigger means. (Or replace "pretrigger mean"
    # with whatever quantity you passed in as <indicator>.)
    # """
using KernelDensity, Distributions, Optim


immutable LinearDriftCorrectParams
    indicator_median::Float64
    param::Float64
end
function drift_correct(indicator, uncorrected, σ_smooth, limit=NaN)
    median_indicator = median(indicator)
    centered_indicator = indicator-median_indicator
    if isnan(limit)
        limit = 1.25*StatsBase.percentile(uncorrected, 99) 
    end   
    edges = pow2edges(σ_smooth*0.4, (0, limit))
    smoother = Normal(0,σ_smooth)
    function tooptimize_entropy(param)
        corrected = uncorrected+centered_indicator*param[1]
        uvh = uvhist(edges, corrected)
        entropy(uvh,σ_smooth)
    end
    optimum = optimize(tooptimize_entropy,-1.0,1.0) # optimize(f(x), lower_bound_x, upper_bound_x) minimizes f with respect to x using brent method
    #check for sucess?
    LinearDriftCorrectParams(median_indicator, optimum.minimum)
end




function pow2edges(stepsize_guess, limits)
    nbins_guess = int(round( 0.5+(last(limits)-first(limits))/stepsize_guess ))
    pow2=1024
    while pow2<nbins_guess
        pow2*=2
    end
    nbins = pow2
    edges = first(limits):(last(limits)-first(limits))/pow2:last(limits)
end
function entropy(h::UnivariateKDE, σ)
    hsmooth=conv(h, Normal(0,σ))
    -sum(StatsBase.xlogx(density(hsmooth)))
end



using Base.Test

@test length(pow2edges(0.1,1:100))==1025
@test isapprox(first(pow2edges(0.1,1:100)),1)
@test isapprox(last(pow2edges(0.1,1:100)),100)

n=10000
indicator = randn(n)*10
uncorrected = 1000+randn(n)+indicator*0.1
drift_correct(indicator, uncorrected, 0.5)
# following Joe Fowler's drift_correct from mass

using KernelDensity, Distributions, Optim

abstract PerPulseTransformation

immutable LinearDriftCorrect <: PerPulseTransformation
    indicator_median::Float64
    param::Float64
end
uvhist(edges, y) = UnivariateKDE(edges, convert(Array{Float64,1},hist(y,edges)[2]))
density(h::UnivariateKDE) = h.density
edges(h::UnivariateKDE) = h.x
Base.midpoints(h::UnivariateKDE) = midpoints(edges(h))
Base.conv(h::Histogram, y) = conv(UnivariateKDE(edges(h), convert(Vector{Float64},counts(h))), y)
"""drift_correct(indicator, uncorrected, σ_smooth, limit=NaN;sloperange=(-0.01, 0.01))
Compute a drift correction that minimizes the spectral entropy of
`uncorrected` (a filtered pulse height vector) after correction.
The `indicator` vector should be of equal length and should have some
linear correlation with the pulse gain. Generally it will be the
pretrigger mean of the pulses, but you can experiment with other choices.

The entropy will be computed on corrected values only in the range
[0, `limit`], so `limit` should be set to a characteristic large
value of `uncorrected`. If `limit` is NaN (the default), then
in will be compute as 25% larger than the 99%ile point of
`uncorrected`

The model is that the filtered pulse height PH should be scaled by
(1 + param*PTM) where `param` is an arbitrary parameter computed here, and
PTM is the difference between each record's pretrigger mean and the
median value of all pretrigger means. (Or replace pretrigger mean
with whatever quantity you passed in as `indicator`.)

`σ_smooth` is used when calculating the entropy.
`first(sloperange)` and `last(sloperange)` are the extreme low and high values
in which an optimum value for `param` is searched for."""
function drift_correct(indicator, uncorrected, σ_smooth, limit=NaN;sloperange=(-0.01, 0.01))
    median_indicator = median(indicator)
    centered_indicator = indicator-median_indicator
    if isnan(limit)
        limit = 1.25*StatsBase.percentile(uncorrected, 99)
    end
    edges = pow2edges(σ_smooth*0.4, (0, limit))
    function tooptimize_entropy(param)
        corrected = uncorrected.*(1+centered_indicator*param[1])
        uvh = uvhist(edges, corrected)
        entropy(uvh,σ_smooth)
    end
    optimum = optimize(tooptimize_entropy,first(sloperange),last(sloperange)) # optimize(f(x), lower_bound_x, upper_bound_x) minimizes f with respect to x using brent method
    #check for sucess?
    LinearDriftCorrect(median_indicator, optimum.minimum)
end



"""pow2edges(stepsize_guess, limits)
Return a range `first(limits):s:last(limits)`. Stepsize `s` is chosen to such that the range has an even number
of midpoints, and `s` is near `stepsize_guess`."""
function pow2edges(stepsize_guess, limits)
    # calculate a range from first(limits) to last(limits), stepsize is chose to be close to guess but given a power of 2 number of midpoints
    nbins_guess = round(Int, 0.5+(last(limits)-first(limits))/stepsize_guess )
    pow2=1024
    while pow2<nbins_guess
        pow2*=2
    end
    nbins = pow2
    first(limits):(last(limits)-first(limits))/pow2:last(limits)
end
"""entropy(h::UnivariateKDE, σ)
Calculate the entropy in spectrum `h` after first smoothing by convolution with a gaussian with std devaition σ."""
function entropy(h::UnivariateKDE, σ)
    hsmooth=conv(h, Normal(0,σ))
    return -sum(StatsBase.xlogx(density(hsmooth)))
end

"""applycorrection(correction::PerPulseTransformation, indicator, uncorrected)
Aplly the correction stored in `correction`, for example a gain correction based on pretrigger median, to
convert `uncorrected` to `corrected` based on the values in `indicator`."""
function applycorrection(correction::LinearDriftCorrect, indicator, uncorrected)
    centered_indicator = indicator-correction.indicator_median
    corrected = uncorrected.*(1+centered_indicator*correction.param)
end

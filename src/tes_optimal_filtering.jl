module TESOptimalFiltering

"Return the autocorrelation of vector data with length `nlags`. Generated by averaging together
the autocorrelations of sequenes of length `seqlen` of points chosen from `data`.
Rejects sequences that contain a point more than `max_excursion` from the mean of `data`.

The preferred way to call is `autocorrelation(data::Vector, nlags::Integer, max_excursion=250)` where `seqlen` is calculated by huersitic."
function autocorrelation(data::Vector, nlags::Integer, seqlen::Integer, max_excursion)
    @assert seqlen >= 2*nlags
    @assert seqlen < length(data)
    nseq = div( length(data), seqlen)
    d = zeros(Float64, seqlen+nlags)
    acorr = zeros(Float64, nlags)
    segments_used = 0

    for i = 1:nseq
        d[1:seqlen] = float(data[(i-1)*seqlen+1:i*seqlen]) - mean(data[(i-1)*seqlen+1:i*seqlen])
        if maximum(d) > max_excursion || minimum(d) < -max_excursion
            continue
        end
        acorr += xcorr(d, d)[seqlen+nlags:seqlen+2*nlags-1]
        segments_used += 1
    end
    if segments_used < max(div(nseq,2), 20)
        error("Rejected too many sequences. Used $segments_used segments out of $nseq available.")
    end
    acorr ./ (segments_used*collect(seqlen:-1:seqlen-nlags+1))
end

"Choose a reasonble sequence length for use in `autocorrelation`."
function autocorr_seqlen_heuristic(nlags::Integer)
    # We'll use sequences at least 7x the # of lags, then
    # round up to one of (1,3,5) times a power of 2
    total =8*nlags
    total = max(total,8)
    lowpow2 = 2^floor(Int, log2(total))
    pow2ratio = total/lowpow2
    if pow2ratio > 1.5
        total = lowpow2*2
    elseif pow2ratio > 1.25
        total = div(lowpow2*3, 2)
    else
        total = div(lowpow2*5, 4)
    end
    @assert total >= 8*nlags
    total-nlags
end

autocorrelation(data::Vector, nlags::Integer, max_excursion=250) = autocorrelation(data, nlags, autocorr_seqlen_heuristic(nlags), max_excursion)

"`calculate_filter(pulse_model, noise_autocorr, f_3db, dt)`
Return filter anf expected variance on pulse height. Filter is rolled off at `f_3db` (hertz) given that the point spacing is `dt` (seconds)"
function calculate_filter(pulse_model, noise_autocorr, f_3db, dt)
    M = hcat(pulse_model, ones(Float64, length(pulse_model))) 
    # M will have 2 columns, or more. 1st column is pulse model. All other columns are
    # the things you want to be insensitive to (like exponentials).
    R = toeplitz(noise_autocorr)
    RinvM = (R \ M)
    A = M' * RinvM # "design matrix"
    variance = inv(A)[1,1] # (dV/V)^2, expected paramter covariance
    filter0 = A \ RinvM'
    filter = smoother(length(pulse_model), f_3db, dt)*vec(filter0[1,:])
    vec(filter), variance
end

"`smoother(n, f_3db, dt)`
Return an (n,n) array `s` such that `s*filter` gives a filter rolled off at `f_3db` (hertz). `dt` is the time (seconds) between samples in `filter`."
function smoother(n, f_3db, dt)
    s = Array(Float64, n, n)
    tau = 1 ./ (2*pi*f_3db*dt)
    for i = 1:n
        for j = 1:n
            s[i,j] = exp(-.5*((i-j)/tau)^2)
        end
        s[i,:] = s[i,:] ./ sum(s[i,:])
    end
    s
end

"Used to efficiently to polyfit in filter5lag."
const SG_filter = [[ -6 24  34 24 -6];
             [-14 -7   0  7 14];
             [ 10 -5 -10 -5 10]]./70.
"`filter5lag(filter, pulses)`
`pulses` can either be one pulse in a `Vector{T<:Number}`, or many pulses in a vector of pulses, or a 2d array where `a[1,:]` is the firse pulse.
Truncate the filter by 4 points, 2 on each end. Dot the trunacated filter with the pulse for 5 different lags. Fit a quadratic to the 5 lags.
Use the peak of the fit as the peak_value and the x center as the phase (aka subsample arrival time)."
function filter5lag{T<:Number}(filter, pulse::Vector{T})
    output = [dot(filter[3:end-2],pulse[i:end-5+i]) for i=1:5]
    a,b,c = SG_filter*output #quadratic fit to 5 lags
    filt_value = a-b^2/c/4 # peak of quadratic
    filt_phase = -b/c/2 # x location of peak
    filt_value, filt_phase
end
function filter5lag{T<:Vector}(filter, pulses::Vector{T})
    filt_values = zeros(Float64, length(pulses))
    filt_phases = zeros(Float64, length(pulses))
    for i in 1:length(pulses)
        filt_values[i], filt_phases[i] = filter5lag(filter, pulses[i])
    end
    filt_values, filt_phases
end
function filter5lag{T<:Number}(filter, pulses::Array{T,2})
    filt_values = zeros(Float64, size(pulses, 2))
    filt_phases = zeros(Float64, size(pulses, 2))
    for j in 1:size(pulses,2)
        output = [dot(filter[3:end-2],pulses[i:end-5+i,j]) for i=1:5]
        pfit = SG_filter*output
        filt_values[j] = pfit[1]-pfit[2]^2/4/pfit[3]
        filt_phases[j] = -pfit[2]/pfit[3]/2
    end
    filt_values, filt_phases
end

"Create and return a toeplitz matrix with c as the first column and
(optional) r as the first row. r[1] will be ignored in
favor of c[1] to determine the diagonal."
function toeplitz{T}(c::Vector{T}, r::Vector{T})
    nr = length(c)
    nc = length(r)
    m = Array(typeof(c[1]), nr, nc)
    for j=1:nc
        for i=1:j-1  # Strict upper triangle
            m[i,j] = r[j-i+1]
        end
        for i=j:nr
            m[i,j] = c[i-j+1]
        end
    end
    m
end
toeplitz{T<:Number}(c::Vector{T}) = toeplitz(c,c)

end # module

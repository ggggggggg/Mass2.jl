using ArrayViews

type PulseSummaries
    pretrig_mean      ::Vector{Float64}
    pretrig_rms       ::Vector{Float64}
    pulse_average     ::Vector{Float64}
    pulse_rms         ::Vector{Float64}
    rise_time         ::Vector{Float64}
    postpeak_deriv    ::Vector{Float64}
    timestamp         ::Vector{Float64}
    peak_index        ::Vector{Int16}
    peak_value        ::Vector{UInt16}
    min_value         ::Vector{UInt16}


    function PulseSummaries(n::Integer)
        pretrig_mean = Array(Float64, n)
        pretrig_rms = Array(Float64, n)
        pulse_average = Array(Float64, n)
        pulse_rms = Array(Float64, n)
        rise_time = Array(Float64, n)
        postpeak_deriv = Array(Float64, n)
        timestamp = Array(Float64, n)
        peak_index = Array(Int16, n)
        peak_value = Array(UInt16, n)
        min_value = Array(UInt16, n)

        new(pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time,
            postpeak_deriv, timestamp, peak_index, peak_value, min_value)
    end
end

function compute_summary(pulseiterator, Npre, frame_time)
    summary = PulseSummaries(length(pulseiterator))
    Nsamp = length(first(pulseiterator))
    Npost = Nsamp-Npre
    post_peak_deriv_vect = zeros(Float64, Npost)
    for (p, data) in enumerate(pulseiterator)
        # Pretrigger computation first
        s = s2 = 0
        min_idx = 0
        peak_idx = 0
        peak_val = UInt16(0)
        min_idx = 0
        min_val = typemax(UInt16)
        for j = 1:Npre
            d=Int(data[j])
            if d > peak_val
                peak_idx, peak_val = j, d
            elseif d < min_val
                min_idx, min_val = j,d
            end
            s+=d
            s2+=d*d
        end
        ptm = s/Npre
        summary.pretrig_mean[p] = ptm
        summary.pretrig_rms[p] = sqrt(abs(s2/Npre - ptm*ptm))

        # Now post-trigger calculations
        s = s2 = 0
        for j = Npre+1:Nsamp
            d=Int(data[j])
            if d > peak_val
                peak_idx, peak_val = j, d
            elseif d < min_val
                min_idx, min_val = j,d
            end
            # d2=d-ptm
            s+=d
            s2+=d*d
        end

        rise_time::Float64 = estimate_rise_time(data, Npre+1:peak_idx,
                                       peak_val, ptm, frame_time)

        postpeak_data = view(data, peak_idx:Nsamp)
        postpeak_deriv = max_timeseries_deriv_simple(postpeak_data)

        # Copy results into the PulseSummaries object
        summary.pulse_average[p] = s/Npost-ptm
        summary.pulse_rms[p] = sqrt(abs(s2/Npost - ptm*(ptm+2*summary.pulse_average[p])))
        summary.pulse_rms[p] = sqrt(abs(s2/Npost - ptm*(ptm+2*summary.pulse_average[p])))
        summary.rise_time[p] = rise_time
        summary.postpeak_deriv[p] = postpeak_deriv
        summary.peak_index[p] = peak_idx
        summary.min_value[p] = min_val
        if ptm < 0
            summary.peak_value[p] = peak_val
        elseif peak_val > ptm
            summary.peak_value[p] = peak_val - round(UInt16,ptm)
        else
            summary.peak_value[p] = UInt16(0)
        end
    end
    (summary.pretrig_mean, summary.pretrig_rms, summary.pulse_average, summary.pulse_rms, summary.rise_time, summary.postpeak_deriv,
	summary.peak_index, summary.peak_value, summary.min_value)
end

"estimate_rise_time(pulserecord, searchrange, peakval, ptm, frametime)
pulserecord = Vector
searchrange = range to look for risetime in, should go from base of pulse to peak, last(searchrange) is peakindex
peakval = maximum(pulserecord[postrig_region], has not had ptm subtracted off)
ptm = pretrigger mean
framtime = time spacing between points"
function estimate_rise_time(pulserecord, searchrange, peakval,ptm,frametime)
    idx10 = first(searchrange)
    peakindex = last(searchrange)
    (peakindex > length(pulserecord) || length(searchrange)==0) && (return length(pulserecord))

    idx90 = peakindex
    thresh10 = 0.1*(peakval-ptm)+ptm
    thresh90 = 0.9*(peakval-ptm)+ptm
    j=0 # to make j exist after the for loop
    for j = 1:peakindex
        if pulserecord[j] > thresh10
            idx10 = j-1
            break
        end
    end
    for j = j+1:peakindex
        if pulserecord[j] > thresh90
            idx90 = j-1
            break
        end
    end
    rise_time = (idx90-idx10)*frametime
    rise_time
end

# Estimate the derivative (units of arbs / sample) for a pulse record or other timeseries.
# This version uses the default kernel of [-2,-1,0,1,2]/10.0
#
max_timeseries_deriv!(deriv, pulserecord, reject_spikes::Bool) =
    max_timeseries_deriv!(deriv, pulserecord, convert(Vector{eltype(deriv)},[.2 : -.1 : -.2]), reject_spikes)


# Post-peak derivative computed using Savitzky-Golay filter of order 3
# and fitting 1 point before...3 points after.
#
max_timeseries_deriv_SG!(deriv, pulserecord, reject_spikes::Bool) =
    max_timeseries_deriv!(deriv, pulserecord, [-0.11905, .30952, .28572, -.02381, -.45238],
                            reject_spikes)

# Estimate the derivative (units of arbs / sample) for a pulse record or other timeseries.
# Caller pre-allocates the full derivative array, which is available as deriv.
# Returns the maximum value of the derivative.
# The kernel should be a short *convolution* (not correlation) kernel to be convolved
# against the input pulserecord.
# If reject_spikes is true, then the max value at sample i is changed to equal the minimum
# of the values at (i-2, i, i+2). Note that this test only makes sense for kernels of length
# 5 (or less), because only there can it be guaranteed insensitive to unit-length spikes of
# arbitrary amplitude.
#
function max_timeseries_deriv!{T}(
        deriv::Vector{T},       # Modified! Pre-allocate an array of sufficient length
        pulserecord, # The pulse record (presumably starting at the pulse peak)
        kernel::Vector{T},      # The convolution kernel that estimates derivatives
        reject_spikes::Bool  # Whether to employ the spike-rejection test
        )
    N = length(pulserecord)
    Nk = length(kernel)
    @assert length(deriv) >= N+1-Nk
    if Nk > N
        return 0.0
    end
    if Nk+4 > N
        reject_spikes = false
    end
    fill!(deriv, zero(eltype(deriv)))
    for i=1:N-Nk+1
        for j=1:Nk
            deriv[i] += pulserecord[i+Nk-j]*kernel[j] #float
        end
    end
    for i=N-Nk+2:length(deriv)
        deriv[i]=deriv[N-Nk+1]
    end
    if reject_spikes
        for i=3:N-Nk-2
            if deriv[i] > deriv[i+2]
                deriv[i] = deriv[i+2]
            end
            if deriv[i] > deriv[i-2]
                deriv[i] = deriv[i-2]
            end
        end
    end
    maximum(deriv)
end

function max_timeseries_deriv_simple(pulserecord)
    max_deriv = realmin(Float64)
    for j = 1:length(pulserecord)-1
        deriv = float(pulserecord[j+1])-float(pulserecord[j])
        deriv > max_deriv && (max_deriv = deriv)
    end
    max_deriv
end

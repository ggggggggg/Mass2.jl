# Glue contains wrapper functions to allow Mass2 steps to call out to other functions with
# convenient syntax. I hope to move functions out of here to better locations. Eventually this file may
# not need to exist.
import TESOptimalFiltering: filter5lag, calculate_filter, autocorrelation
using LsqFit

function compute_whitenoise_filter(pulse, selection_good)
	sumpulse = zeros(Int64, length(pulse[1]))
	n=0
	for i=1:length(selection_good) # take the mean of the pulses in a way that avoids int overflows
		if selection_good[i]
			sumpulse+=pulse[i]
			n+=1
		end
	end
	meanpulse = sumpulse/n
	filter = meanpulse-mean(meanpulse) # 0 average
	normalization = (maximum(filter)-minimum(filter))./dot(filter, filter)
	filter*normalization
end
function filter1lag(pulse, filter)
	out = Array(Float64, length(pulse))
	for i=1:length(pulse)
		out[i] = dot(pulse[i], filter)
	end
	out
end

function selectfromcriteria(x...) # placeholder, kinda ugly to use and probalby a bit slow
	iseven(length(x)) || error("x must be indicator,criteria,indicator,criteria...")
	out = trues(length(x[1]))
	for i = 1:2:length(x)
		low, high = x[i+1]
		out &= low .< x[i] .< high
	end
	out
end

function calc_dc(indicator, uncorrected, selection_good)
	σ_smooth = 1.0
	drift_correct(indicator[:][selection_good], uncorrected[:][selection_good], σ_smooth)
end

function two_exponential_pulse_for_fittng(x,p)
	# check for bad values of rise_point and fall_points
	# return vector record containing points at the quiesence value
	if p[1]<0 || p[2]<0 || p[1]==p[2]
		return zeros(Float64,length(x))
	end
	two_exponential_pulses(length(x), p[1], p[2],p[3], (p[4],), (p[5],))
end
function fit_pulse_two_exponential(pulse, pretrig_nsamples, frametime)
	pguess = zeros(Float64,5)
	pguess[1] = length(pulse)/10 # rise points
	pguess[2] = 3*pguess[1] # fall points
	pguess[3] = pulse[1] # quiescent_value
	pguess[4] = pretrig_nsamples # arrival time (really should be based on ljh)
	pguess[5] = maximum(pulse)-pulse[1] # amplitude
	fit = curve_fit(two_exponential_pulse_for_fittng, collect(1:length(pulse)), pulse, pguess)
	rise_tau_s, fall_tau_s = fit.param[1]*frametime, fit.param[2]*frametime
end


"selectfromcriteria(x...)
selectfromcriteria(indicator1, criteria1, indicator2, criteria2,...) takes 1 or more indicator,criteria inputs and
returns a bitvector that is true if all indictors pass criteria[1] .< indicator .< criteria[2]"
function selectfromcriteria(x...) # placeholder, kinda ugly to use and probalby a bit slow
	iseven(length(x)) || error("x must be indicator,criteria,indicator,criteria...")
	out = trues(length(x[1]))
	for i = 1:2:length(x)
		low, high = x[i+1]
		out &= low .< x[i] .< high
	end
	out
end

"readcontinuous(ljh,maxnpoints=50_000_000)
read from an LJH file, return pulses concatonated together to form a continuour vector
reads either the whole file, or until maxnpoints"
function readcontinuous(ljh,maxnpoints=50_000_000)
	nsamp = record_nsamples(ljh)
	maxnpoints = nsamp*div(maxnpoints, nsamp)
	npoints = min(maxnpoints, nsamp*length(ljh))
	data = zeros(UInt16,npoints)
	for (i,r) in enumerate(ljh[1:div(npoints, nsamp)])
		data[1+(i-1)*nsamp:i*nsamp] = r.data
	end
	data
end

"compute_noise_autocorr(ljh_filename, nsamp)
computer noise autorcorrelation of length namp from an ljh file of name ljh_filename"
function compute_noise_autocorr(ljh_filename, nsamp)
	data = readcontinuous(LJHGroup(ljh_filename))
	autocorrelation(data, nsamp, 1000) # final argument is max_excursion
	# this needs a revamp to reject fewer noise records and also return delta_f
end

"compute_average_pulse(pulse, selection_good)
compute the average (mean) the good pulses in pulse. good determined by bitvector selection_good"
function compute_average_pulse(pulse, selection_good)
	sumpulse = zeros(Int64, length(pulse[1]))
	n=0
	for i=1:length(selection_good) # take the mean of the pulses in a way that avoids int overflows
		if selection_good[i]
			sumpulse+=pulse[i]
			n+=1
		end
	end
	meanpulse = sumpulse/n
end

"compute_filter(average_pulse, noise_autocorr, f_3db, dt)\n
compute filter given average_pulse and noise_autorr of same length, f_3db in hz, dt in seconds
return `filter`, `vdv`"
function compute_filter(average_pulse, noise_autocorr, f_3db, dt)
	filter, variance = calculate_filter(average_pulse, noise_autocorr, f_3db, dt)
	filter*=maximum(average_pulse)-minimum(average_pulse) # like mass, normalize so you get pulseheight similar to raw units
	vdv = 1/sqrt(variance)/2.355
	filter,vdv
end

function estimate_pretrig_rms_and_postpeak_deriv_criteria(fname, pretrig_nsamples)
	# using the noise file we will try to estimate values to use for actual cuts on pulses
	# pretrig_rms and post_peak deriv should be fairly well estiamated by the noise file
	ljh = LJHGroup(fname)

	traces=[r.data for r in ljh[1:end]]
	pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv,
		peak_index, peak_value, min_value=compute_summary(traces,pretrig_nsamples,frametime(ljh))


	# the distribution of pretrig_rms should follow a chisq distributions if the points are independent (IID) and have gaussian noise
	# but the chisq distribution for 256 degrees of freedom looks a lot like a gaussian, so we're just going to stick with the
	# gaussian approximation for now, plus a backup measure that says at most cut 1 % of pulses
	n_std = 10
	pretrig_rms_max = max(median(pretrig_rms)+n_std*std(pretrig_rms), StatsBase.percentile(pretrig_rms,99))
	pretrig_rms_min = max(0, median(pretrig_rms)-n_std*std(pretrig_rms))


	postpeak_deriv_max = max(median(postpeak_deriv)+n_std*std(postpeak_deriv), StatsBase.percentile(postpeak_deriv,99))
	postpeak_deriv_min = median(postpeak_deriv)-n_std*std(postpeak_deriv)
	(pretrig_rms_min, pretrig_rms_max), (postpeak_deriv_min, postpeak_deriv_max)
end

function estimate_peak_index_criteria(peak_index)
	mad = mean(abs(peak_index-median(peak_index)))
	med = median(peak_index)
	(med-10*mad, med+10*mad)
end

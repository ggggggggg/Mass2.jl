
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
	println("inside calc_dc")
	σ_smooth = 1.0
	drift_correct(indicator[:][selection_good], uncorrected[:][selection_good], σ_smooth)
end

# generate a spectra with multiple peaks for calibration testing
dist_energies = [4000,4700,5500,6000,6800,7500,8200,9600,10000]
dist_filt_values = 5*dist_energies.^0.8
dist_use_inds = [1,3,6,7,8,9]
dist_use_filt_values = dist_filt_values[dist_use_inds]
dist_use_energies = dist_energies[dist_use_inds]
filt_value_distribution = MixtureModel([Normal(mu,1) for mu in dist_filt_values])

steps = AbstractStep[]

push!(steps, MockPulsesStep(TupacLikeTwoExponentialPulseGenerator(Int, filt_value_distribution), 5000, [:pulse,:rowstamp]))
push!(steps, PerPulseStep(compute_summary, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :peak_index, :peak_index_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
push!(steps, PerPulseStep(filter1lag, [:pulse, :whitenoise_filter], [:filt_value]))
push!(steps, HistogramStep(update_histogram!, [:filt_value_hist, :selection_good, :filt_value]))
push!(steps, ThresholdStep(calibrate_nofit, [:filt_value_dc_hist,:known_energies],[:calibration],:filt_value_dc_hist, counted, 5000, true))
push!(steps, ThresholdStep(compute_whitenoise_filter, [:pulse, :selection_good], [:whitenoise_filter], :selection_good, sum, 100, true))
push!(steps, ThresholdStep(calc_dc, [:pretrig_mean, :filt_value, :selection_good], [:ptm_dc],:filt_value_hist, counted, 3000, true))
push!(steps, PerPulseStep(applycorrection, [:ptm_dc, :pretrig_mean, :filt_value], [:filt_value_dc]))
push!(steps, HistogramStep(update_histogram!, [:filt_value_dc_hist, :selection_good, :filt_value_dc]))
push!(steps, PerPulseStep(apply_calibration, [:calibration, :filt_value_dc], [:energy]) )
push!(steps, HistogramStep(update_histogram!, [:energy_hist, :selection_good, :energy]))

push!(steps, ToJLDStep([:filt_value,:pretrig_rms, :energy]))
# push!(steps, FreeMemoryStep())

push!(perpulse_symbols, :filt_value, :selection_good, :energy, :pulse, :rowstamp,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good, :filt_value_dc)

c=Channel()
c[:pretrig_mean] = RunningVector(Float64)
c[:pretrig_rms] = RunningVector(Float64)
c[:pulse_average] = RunningVector(Float64)
c[:pulse_rms] = RunningVector(Float64)
c[:rise_time] = RunningVector(Float64)
c[:postpeak_deriv] = RunningVector(Float64)
c[:peak_index] = RunningVector(Uint16)
c[:peak_value] = RunningVector(Uint16)
c[:min_value] = RunningVector(Uint16)
c[:filt_value] = RunningVector(Float64)
c[:filt_value_dc] = RunningVector(Float64)
c[:pretrig_rms] = RunningVector(Float64)
c[:selection_good] = RunningSumBitVector()
c[:energy] = RunningVector(Float64)
c[:filt_value_hist] = Histogram(0:1:20000)
c[:filt_value_dc_hist] = Histogram(0:1:20000)
c[:energy_hist] = Histogram(0:1:20000)
c[:pulse] = RunningVector(Vector{Int})
c[:rowstamp] = RunningVector(Int)
c[:pre_samples] = 100
c[:frame_time] = 1/100000
c[:peak_index_criteria] = (170,200)
c[:pretrig_rms_criteria] = (0.0,10.)
c[:postpeak_deriv_criteria] = (0.0,50.0)
c[:known_energies] = dist_use_energies

g = graph(steps)


workstat(n, s::MockPulsesStep, t) = "MockPulsesStep "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::PerPulseStep, t) = "PerPulse:$(graphlabel(s)) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ToJLDStep, t) = "ToJLDStep $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::HistogramStep, t) = "HistogramStep:$(inputs(s)[1]) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ThresholdStep, t) = "ThresholdStep:$(graphlabel(s)) $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::FreeMemoryStep, t) = "FreeMemoryStep $n executions at "*@sprintf("%0.2f executions/s",n/t)

savegraph("graph",g)
workdone_cumulative = Array(Int, length(steps))
stepelapsed_cumulative = Array(Float64, length(steps))
errors = Any[]
close(jldopen(filename(c),"w")) # wipe the test file
for i = 1:4
	println("** loop iteration $i")
	#do steps
	for (j,s) in enumerate(steps)
		tstart = time()
		try
			wu = workunits(dostep(s,c))
			stepelapsed = time()-tstart
			workdone_cumulative[j] += wu 
			stepelapsed_cumulative[j] += stepelapsed
		catch ex
			println((i,j))
			println(sprint(io->showerror(io, ex)))
			println(sprint(io->Base.show_backtrace(io, catch_backtrace())))
			debug(s,c)
			break
		end

	end

end
[println("step $i: "*workstat(workdone_cumulative[i], steps[i], stepelapsed_cumulative[i])) for i = 1:length(steps)];
[println("$n donethru $(donethru(c[n]))") for n in keys(c)];
jld = jldopen(filename(c), "r+")

using Base.Test

eh = c[:energy_hist]
@test counted(eh) == sum(counts(eh))
#@test counted(eh) == sum(c[:selection_good]) # this test can easily be failed when things are working as intended, just by changing the step order
# but I'd like to enforce the step order for now so that this test makes sense.

# peakinds returns from lowest counts to highest counts
# lets take the largest peaks and match them up to their known energies from distenergies
peakinds = findpeaks(eh, fwhm=15, n=length(dist_energies))
found_energies = sort(bin_centers(eh)[peakinds])
fvh = c[:filt_value_hist]
peakinds_filt_value = findpeaks(fvh, fwhm=15, n=length(dist_energies))
found_filt_values = sort(bin_centers(fvh)[peakinds])
cal = c[:calibration]

# filt value agrees with manual calculation
fvmanual = [dot(p, c[:whitenoise_filter]) for p in c[:pulse][:]]
@test fvmanual == c[:filt_value][:] 

@test all(abs(found_energies-dist_energies).<15) # check that the calibration has the correct peak assignment
close(jld)

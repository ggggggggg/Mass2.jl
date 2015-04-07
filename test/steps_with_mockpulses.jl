using Mass2

function compute_whitenoise_filter(pulse, selection_good) 
	filter = mean(pulse[selection_good])
	filter -= mean(filter) # 0 average
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

# generate a spectra with multiple peaks for calibration testing
dist_energies = [4000,4700,5500,6000,6800,7500,8200,9600,10000]
dist_filt_values = 5*dist_energies.^0.8
dist_use_inds = [2,3,5,9]
dist_use_filt_values = dist_filt_values[dist_use_inds]
dist_use_energies = dist_energies[dist_use_inds]
filt_value_distribution = MixtureModel([Normal(mu,1) for mu in dist_filt_values])

steps = AbstractStep[]

push!(steps, MockPulsesStep(TupacLikeTwoExponentialPulseGenerator(Int, filt_value_distribution), 1000, [:pulse,:rowstamp]))
push!(steps, PerPulseStep(compute_summary, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
push!(steps, PerPulseStep(apply_calibration, [:calibration, :filt_value], [:energy]) )
push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :peak_index, :peak_index_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
push!(steps, ToJLDStep([:filt_value,:pretrig_rms, :energy]))
push!(steps, HistogramStep(update_histogram!, [:filt_value_hist, :selection_good, :filt_value]))
push!(steps, HistogramStep(update_histogram!, [:energy_hist, :selection_good, :energy]))
push!(steps, ThresholdStep(calibrate_nofit, [:filt_value_hist,:known_energies],[:calibration],:filt_value_hist, counted, 5000, true))
push!(steps, ThresholdStep(compute_whitenoise_filter, [:pulse, :selection_good], [:whitenoise_filter], :selection_good, sum, 100, true))
push!(steps, PerPulseStep(filter1lag, [:pulse, :whitenoise_filter], [:filt_value]))
push!(steps, FreeMemoryStep())

push!(perpulse_symbols, :filt_value, :selection_good, :energy, :pulse, :rowstamp,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good)

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
c[:pretrig_rms] = RunningVector(Float64)
c[:selection_good] = RunningSumBitVector()
c[:energy] = RunningVector(Float64)
c[:filt_value_hist] = Histogram(0:1:20000)
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
for i = 1:10
	println("** loop iteration $i")
	#do steps
	for (j,s) in enumerate(steps)
		tstart = time()
		# try
		wu = workunits(dostep(s,c))
		stepelapsed = time()-tstart
		workdone_cumulative[j] += wu 
		stepelapsed_cumulative[j] += stepelapsed
		# catch ex
		# 	showerror(STDOUT, ex, backtrace())
		# 	push!(errors, ex)
		# end
		# workdone_cumulative+=workdone
		# stepelapsed_cumulative+=stepelapsed
	end

end
[println("step $i: "*workstat(workdone_cumulative[i], steps[i], stepelapsed_cumulative[i])) for i = 1:length(steps)];
[println("$n donethru $(donethru(c[n]))") for n in keys(c)];
jld = jldopen(filename(c), "r+")



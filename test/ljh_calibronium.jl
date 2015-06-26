using Mass2
using ReferenceMicrocalFiles



# type GetPulsesStep{T} <: AbstractStep
# 	pulse_source::T
# 	outputs::Vector{Symbol}
# 	previous_pulse_index::Int # used to keep track of location in an ljh file for example
# 	max_pulses_per_step::Int
# end
# inputs(s::GetPulsesStep) = Symbol[]
# function dostep(s::GetPulsesStep{LJHGroup},c::Channel)
# 	r = s.previous_pulse_index+1:min(s.previous_pulse_index+s.max_pulses_per_step, length(s.pulse_source))
# 	length(r)==0 && (return r)
# 	pulses, rowstamps = collect(s.pulse_source[r])
# 	pulses_out, rowstamps_out = perpulse_outputs(s,c)
# 	append!(pulses_out, pulses)
# 	assert(length(pulses_out)==last(r))
# 	append!(rowstamps_out, rowstamps)
# 	assert(length(rowstamps_out)==last(r))
# 	s.previous_pulse_index=last(r)
# 	r
# end
# graphlabel(s::GetPulsesStep) = repr(typeof(s))

# type LJHPulseEmitter
# 	ljh_group::LJHGroup
# 	max_pulses_per_step::Int
# 	donethru::Int
# end
# type LJHPulseEmitterSlice
# 	ljhpulsemitter::LJHPulseEmitter
# 	r::UnitRange{Int64}
# end
# LJHPulseEmitter(filename,max_pulses_per_step=1000) = LJHPulseEmitter(LJHGroup(filename), max_pulses_per_step, max_pulses_per_step)
# donethru(x::LJHPulseEmitter) = x.donethru
# Base.getindex(x::LJHPulseEmitter, r::UnitRange{Int64}) = LJHPulseEmitterSlice(x, r)
# function getpulses!(ljhpulsemitterslight::LJHPulseEmitterSlice, max_pulses_per_step::Int)
# 	pulses = Array(Vector{Uint16},length(ljh_group_slice))
# 	rowstamps = Array(Uint64, length(ljh_group_slice))
# 	ljh_group_slice = ljhpulsemitterslice.ljhpulsemitter.ljh_group[r]
# 	ljhpulsemitterslice.ljhpulsemitter.donethru = last(ljhpulsemitter.r)
# 	for (i,(pulse, rowstamp)) in enumerate(ljh_group_slice)
# 		pulses[i]=pulse
# 		rowstamps[i]=rowstamp
# 	end
# 	pulses, rowstamps
# end

ljhpaths = ReferenceMicrocalFiles.getljhpaths()
ljh = LJHGroup(ljhpaths[1])


steps = AbstractStep[]
# generate a spectra with multiple peaks for calibration testing
dist_energies = [4000,4700,5500,6000,6800,7500,8200,9600,10000]
dist_filt_values = 5*dist_energies.^0.8
dist_use_inds = [1,3,6,7,8,9]
dist_use_filt_values = dist_filt_values[dist_use_inds]
dist_use_energies = dist_energies[dist_use_inds]
filt_value_distribution = MixtureModel([Normal(mu,1) for mu in dist_filt_values])
push!(steps, GetPulsesStep(ljh, [:pulse, :rowstamp], 0,500))
s=steps[1]

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

# push!(steps, ToJLDStep([:filt_value,:pretrig_rms, :energy]))
# push!(steps, FreeMemoryStep())

push!(perpulse_symbols, :filt_value, :selection_good, :energy, :pulse, :rowstamp,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good, :filt_value_dc,:ljh_group)

c=Channel()
# c[:ljh_group] = ljh
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
c[:pulse] = RunningVector(Vector{Uint16})
c[:rowstamp] = RunningVector(Int)
c[:pre_samples] = 100
c[:frame_time] = 1/100000
c[:peak_index_criteria] = (2280,2360)
c[:pretrig_rms_criteria] = (0.0,10.)
c[:postpeak_deriv_criteria] = (0.0,100.0)
c[:known_energies] = [4952, 5898, 6929, 8040]

workstat(n, s::MockPulsesStep, t) = "MockPulsesStep "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::PerPulseStep, t) = "PerPulse:$(graphlabel(s)) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::GetPulsesStep, t) = "GetPulsesStep:$(graphlabel(s)) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ToJLDStep, t) = "ToJLDStep $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::HistogramStep, t) = "HistogramStep:$(inputs(s)[1]) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ThresholdStep, t) = "ThresholdStep:$(graphlabel(s)) $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::FreeMemoryStep, t) = "FreeMemoryStep $n executions at "*@sprintf("%0.2f executions/s",n/t)

workdone_cumulative = zeros(Int, length(steps))
stepelapsed_cumulative = Array(Float64, length(steps))
errors = Any[]
close(jldopen(filename(c),"w")) # wipe the test file
for i = 1:30
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
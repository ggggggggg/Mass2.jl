using Mass2, Base.Test


# generate a spectra with multiple peaks for calibration testing
dist_energies = [4000,4700,5500,6000,6800,7500,8200,9600,10000]
dist_filt_values = 5*dist_energies.^0.8
dist_use_inds = [1,3,6,7,8,9]
dist_use_filt_values = dist_filt_values[dist_use_inds]
dist_use_energies = dist_energies[dist_use_inds]
filt_value_distribution = MixtureModel([Normal(mu,1) for mu in dist_filt_values])

steps = AbstractStep[]

push!(steps, MockPulsesStep(TupacLikeTwoExponentialPulseGenerator(Int, filt_value_distribution), 100, 5000, [:pulse,:rowstamp]))
push!(steps, PerPulseStep(compute_summary, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :peak_index, :peak_index_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
push!(steps, PerPulseStep(filter1lag, [:pulse, :whitenoise_filter], [:filt_value]))
push!(steps, HistogramStep(update_histogram!, [:filt_value_hist, :selection_good, :filt_value]))
push!(steps, ThresholdStep(calibrate_nofit, [:filt_value_dc_hist,:known_energies],[:calibration],:filt_value_dc_hist, counted, 3000, true))
push!(steps, ThresholdStep(compute_whitenoise_filter, [:pulse, :selection_good], [:whitenoise_filter], :selection_good, sum, 100, true))
push!(steps, ThresholdStep(calc_dc, [:pretrig_mean, :filt_value, :selection_good], [:ptm_dc],:filt_value_hist, counted, 3000, true))
push!(steps, PerPulseStep(applycorrection, [:ptm_dc, :pretrig_mean, :filt_value], [:filt_value_dc]))
push!(steps, HistogramStep(update_histogram!, [:filt_value_dc_hist, :selection_good, :filt_value_dc]))
push!(steps, PerPulseStep(apply_calibration, [:calibration, :filt_value_dc], [:energy]) )
push!(steps, HistogramStep(update_histogram!, [:energy_hist, :selection_good, :energy]))

#push!(steps, ToJLDStep([:filt_value,:pretrig_rms, :energy]))
# push!(steps, FreeMemoryStep())

push!(perpulse_symbols, :filt_value, :selection_good, :energy, :pulse, :rowstamp,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv,
	:peak_index, :peak_value, :min_value, :selection_good, :filt_value_dc)

mc=MassChannel()
mc[:pretrig_mean] = RunningVector(Float64)
mc[:pretrig_rms] = RunningVector(Float64)
mc[:pulse_average] = RunningVector(Float64)
mc[:pulse_rms] = RunningVector(Float64)
mc[:rise_time] = RunningVector(Float64)
mc[:postpeak_deriv] = RunningVector(Float64)
mc[:peak_index] = RunningVector(UInt16)
mc[:peak_value] = RunningVector(UInt16)
mc[:min_value] = RunningVector(UInt16)
mc[:filt_value] = RunningVector(Float64)
mc[:filt_value_dc] = RunningVector(Float64)
mc[:pretrig_rms] = RunningVector(Float64)
mc[:selection_good] = RunningSumBitVector()
mc[:energy] = RunningVector(Float64)
mc[:filt_value_hist] = Histogram(0:1:20000)
mc[:filt_value_dc_hist] = Histogram(0:1:20000)
mc[:energy_hist] = Histogram(0:1:20000)
mc[:pulse] = RunningVector(Vector{Int})
mc[:rowstamp] = RunningVector(Int)
mc[:pre_samples] = 100
mc[:frame_time] = 1/100000
mc[:peak_index_criteria] = (170,200)
mc[:pretrig_rms_criteria] = (0.0,10.)
mc[:postpeak_deriv_criteria] = (0.0,50.0)
mc[:known_energies] = dist_use_energies

# metadata
mc[:steps]=steps
mc[:workdone_cumulative] = Dict{AbstractStep, Int64}()
mc[:time_elapsed_cumulative] = Dict{AbstractStep, Float64}()
mc[:workdone_last] = Dict{AbstractStep, Int64}()
mc[:hdf5_filename] = "steps_with_mockpulses.jld"
mc[:oncleanfinish] = (mc)->nothing
make_task(mc)

workstat(n, s::MockPulsesStep, t) = "MockPulsesStep "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::PerPulseStep, t) = "PerPulse:$(graphlabel(s)) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ToJLDStep, t) = "ToJLDStep $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::HistogramStep, t) = "HistogramStep:$(inputs(s)[1]) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ThresholdStep, t) = "ThresholdStep:$(graphlabel(s)) $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::FreeMemoryStep, t) = "FreeMemoryStep $n executions at "*@sprintf("%0.2f executions/s",n/t)

savegraph("graph",graph(steps))


isfile(mc[:hdf5_filename]) && rm(mc[:hdf5_filename])
schedule(mc)
sleep(2) # make sure steps have time do some work before planning to end, test fails without this
# may need better system for this
plan_to_end(mc)
wait(mc[:task])

eh = mc[:energy_hist]
@test counted(eh) == sum(counts(eh))

# peakinds returns from lowest counts to highest counts
# lets take the largest peaks and match them up to their known energies from distenergies
peakinds = findpeaks(eh, fwhm=15, n=length(dist_energies))
found_energies = sort(bin_centers(eh)[peakinds])
fvh = mc[:filt_value_hist]
peakinds_filt_value = findpeaks(fvh, fwhm=15, n=length(dist_energies))
found_filt_values = sort(bin_centers(fvh)[peakinds])
cal = mc[:calibration]

# filt value agrees with manual calculation
fvmanual = [dot(p, mc[:whitenoise_filter]) for p in mc[:pulse][1:end]]
@test fvmanual == mc[:filt_value][1:end]

@test all(abs(found_energies-dist_energies).<25) # check that the calibration has the correct peak assignment

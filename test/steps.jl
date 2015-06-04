steps = AbstractStep[]

f() = nothing
push!(steps, PerPulseStep(f, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
push!(steps, PerPulseStep(f, [:pretrig_rms, :pretrig_rms_criteria, :peak_index, :peak_index_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
push!(steps, PerPulseStep(f, [:pulse, :whitenoise_filter], [:filt_value]))
push!(steps, ThresholdStep(f, [:pretrig_mean, :filt_value, :selection_good], [:ptm_dc],:filt_value_hist, counted, 3000, true))

g = Graphs.graph(steps)
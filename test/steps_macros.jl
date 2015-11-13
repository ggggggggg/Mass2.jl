using Base.Test
cmp{T}(x::T,y::T) = all(Bool[getfield(x,fname) == getfield(y,fname) for fname in fieldnames(x)])

#@perpulse a,b,c,d = do_filter(filter, pulse_record)
#@threshold calibration = calibrate_nofit(filt_value_dc_hist, known_energies) when counted(filt_value_dc_hist)>3000
#macroexpand( quote @histogram update_histogram!(filt_value_hist, selection_good, filt_value) end)
@test cmp(PerPulseStep(sum,[:filter,:pulse_record],[:filt_value]), @perpulse filt_value = sum(filter, pulse_record))
@test cmp(ThresholdStep(calibrate_nofit, [:filt_value_dc_hist,:known_energies],[:calibration],:filt_value_dc_hist, counted, 3000, true),
@threshold calibration = calibrate_nofit(filt_value_dc_hist, known_energies) when counted(filt_value_dc_hist)>3000 )
@test cmp(HistogramStep(update_histogram!, [:selection_good, :filt_value], [:filt_value_hist]),
@histogram update_histogram!(filt_value_hist, selection_good, filt_value))

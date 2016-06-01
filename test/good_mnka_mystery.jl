using Mass2
using ReferenceMicrocalFiles
using ZMQ
using Base.Test



function setup_channel(ljh_filename, noise_filename)
	ljh = LJHGroup(ljh_filename)

	mc=MassChannel()

	mc[:pretrig_mean] = RunningVector(Float32)
	mc[:pretrig_rms] = RunningVector(Float32)
	mc[:pulse_average] = RunningVector(Float32)
	mc[:pulse_rms] = RunningVector(Float32)
	mc[:rise_time] = RunningVector(Float32)
	mc[:postpeak_deriv] = RunningVector(Float32)
	mc[:peak_index] = RunningVector(UInt16)
	mc[:peak_value] = RunningVector(UInt16)
	mc[:min_value] = RunningVector(UInt16)
	mc[:filt_value] = RunningVector(Float32)
	mc[:filt_value_dc] = RunningVector(Float32)
	mc[:filt_phase] = RunningVector(Float32)
	mc[:energy] = RunningVector(Float32)
	mc[:pulse] = RunningVector(Vector{UInt16})
	mc[:rowcount] = RunningVector(Int)
	mc[:timestamp_posix_usec] = RunningVector(Int)
	mc[:selection_good] = RunningSumBitVector()
	mc[:filt_value_hist] = Histogram(0:1:20000)
	mc[:filt_value_dc_hist] = Histogram(0:1:20000)
	mc[:energy_hist] = Histogram(0:1:20000)

	mc[:pretrig_nsamples] = LJH.pretrig_nsamples(ljh)
	mc[:samples_per_record] = LJH.record_nsamples(ljh)
	mc[:frametime] = LJH.frametime(ljh)
	#mc[:rise_time_criteria] = (0,0.0006) # from exafs.basic_cuts, but it is in ms there, s here
	#mc[:pretrig_rms_criteria] = (0.0,30.) # from exafs.basic_cuts
	#mc[:postpeak_deriv_criteria] = (0.0,250.0) # from exafs.basic_cuts
	mc[:f_3db] = 10000
	mc[:known_energies] = [5898.801,6490.59] # FeKAlpha and FeKBeta

	#metadata
	mc[:calibration_nextra] = 1 # when finding peaks, how many peaks other than the largest n to include when assigning peaks to energies

	mc[:noise_filename]=ascii(noise_filename)
	mc[:ljh_filename]=ascii(ljh_filename) # h5py can't read UFT8String written by julia
	mc[:name] = "summarize and filter test"
	mc[:hdf5_filename] = "$(splitext(ljh_filename)[1])_jl.hdf5"
	mc[:oncleanfinish] = markhdf5oncleanfinish
	isfile(mc[:hdf5_filename]) && rm(mc[:hdf5_filename])

	steps = AbstractStep[
	GetPulsesStep(ljh, [:pulse, :rowcount, :timestamp_posix_usec], 0,100)
	@perpulse pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv, peak_index, peak_value, min_value = compute_summary(pulse, pretrig_nsamples, frametime)
	@threshold pretrig_rms_criteria, postpeak_deriv_criteria = estimate_pretrig_rms_and_postpeak_deriv_criteria(noise_filename, pretrig_nsamples) when length(pulse) > 100
	@threshold peak_index_criteria = estimate_peak_index_criteria(peak_index) when length(peak_index)>100
	@perpulse selection_good = selectfromcriteria(pretrig_rms, pretrig_rms_criteria, peak_index, peak_index_criteria, postpeak_deriv, postpeak_deriv_criteria)
	@threshold average_pulse = compute_average_pulse(pulse, selection_good) when sum(selection_good) > 100
	@threshold rise_tau_s, fall_tau_s = fit_pulse_two_exponential(average_pulse, pretrig_nsamples, frametime) when length(filt_value)>10
	@threshold noise_autocorr = compute_noise_autocorr(noise_filename, samples_per_record) when sum(selection_good) > 100
	@threshold filter, vdv = compute_filter(average_pulse, noise_autocorr, f_3db, frametime) when sum(selection_good) > 100
	@perpulse filt_value, filt_phase = filter5lag(filter, pulse)
	@threshold pretrigger_mean_drift_correction_params = calc_dc(pretrig_mean, filt_value, selection_good) when sum(selection_good)>2500
	@perpulse filt_value_dc = applycorrection(pretrigger_mean_drift_correction_params, pretrig_mean, filt_value)
	@threshold calibration = calibrate_nofit(filt_value_dc_hist, known_energies, calibration_nextra) when counted(filt_value_dc_hist)>2500
	@perpulse energy = apply_calibration(calibration, filt_value_dc)
	@histogram update_histogram!(filt_value_hist, selection_good, filt_value)
	@histogram update_histogram!(energy_hist, selection_good, energy)
	@histogram update_histogram!(filt_value_dc_hist, selection_good, filt_value_dc)
	ToJLDStep([:filt_value, :filt_value_dc, :energy, :filt_phase, :pretrig_rms, :postpeak_deriv, :rise_time, :peak_index,
	:pretrig_mean, :pulse_average, :pulse_rms, :peak_value, :min_value, :rowcount, :timestamp_posix_usec],
	Pair[:filter=>"filter/filter", :f_3db=>"filter/f_3db", :frametime=>"filter/frametime", :noise_autocorr=>"filter/noise_autocorr", :average_pulse=>"filter/average_pulse",
	:average_pulse=>"average_pulse",
	:samples_per_record=>"samples_per_record", :frametime=>"frametime", :pretrig_nsamples=>"pretrig_nsamples",
	:ljh_filename=>"ljh_filename", :noise_filename=>"noise_filename",
	:peak_index_criteria=>"selection_criteria/peak_index", :pretrig_rms_criteria=>"selection_criteria/pretrig_rms", :postpeak_deriv_criteria=>"selection_criteria/postpeak_deriv"],
	mc[:hdf5_filename])
	MemoryLimitStep(Int(4e6)) # throw error if mc uses more than 4 MB
	FreeMemoryStep()
	] # end steps
	# write a verification function that makes sure all inputs either exist, or are the output of another step
	setsteps!(mc, steps)
	preptasks!(mc)
	mc
end


masschannels, watcher_exitchannel, watcher_task = schedule_MATTER_watcher(setup_channel,1)

rmf = ReferenceMicrocalFiles.dict["good_mnka_mystery"]

t = @schedule begin
getopenfilelimit() = parse(Int,split(split(readall(`ulimit -a`),"\n")[6])[end])
#if getopenfilelimit()>=1000
	sleep(2)
	LJHUtil.write_sentinel_file(rmf.noise_filename,false)
	sleep(2)
	LJHUtil.write_sentinel_file(rmf.filename,false)
	sleep(2)
	LJHUtil.write_sentinel_file(rmf.noise_filename,false)
	sleep(2)
	put!(watcher_exitchannel,1)
	@schedule begin
		tasks = [mc.task.value for (ch, mc) in masschannels]
		for task in tasks # wait for all tasks to finish
			try wait(task) end # wait rethrows errors from failed tasks, we can always look at them later
		end
		ndone = sum([task.state==:done for task in tasks])
		nfailed = sum([task.state==:failed for task in tasks])
		info("all tasks done!! Yay. $ndone :done and $nfailed :failed")
	end;
#else
#	println("open file limit is too low, try ulimit -n 1000")
#end
end # write to MATTER sentinel file to simulate matter writing various files


wait(t)
mc=masschannels[13];

wait(mc)
@test get(mc.task).state == :done
@test seen(mc[:energy_hist])==length(mc[:selection_good])
@test seen(mc[:filt_value_hist])==length(mc[:selection_good])
@test seen(mc[:filt_value_dc_hist])==length(mc[:selection_good])
@test 0.8*sum(mc[:selection_good])<counted(mc[:energy_hist])<sum(mc[:selection_good])
@test 0.8*sum(mc[:selection_good])<counted(mc[:filt_value_hist])<sum(mc[:selection_good])
@test 0.8*sum(mc[:selection_good])<counted(mc[:filt_value_dc_hist])<sum(mc[:selection_good])

h5open(mc[:hdf5_filename]) do h5
	@test read(h5["selection_criteria/peak_index"]) == mc[:peak_index_criteria]
	@test read(h5["selection_criteria/pretrig_rms"]) == mc[:pretrig_rms_criteria]
	@test read(h5["selection_criteria/postpeak_deriv"]) == mc[:postpeak_deriv_criteria]
	@test length(h5["filt_value"]) == length(mc[:filt_value])
end

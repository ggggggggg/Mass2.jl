using Mass2
import TESOptimalFiltering: filter5lag, calculate_filter, autocorrelation
using ReferenceMicrocalFiles
using ZMQ
using LsqFit
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
	nsamp = LJH.record_nsamples(ljh)
	maxnpoints = nsamp*div(maxnpoints, nsamp)
	npoints = min(maxnpoints, nsamp*length(ljh))
	data = zeros(UInt16,npoints)
	for (i,(p,t)) in enumerate(ljh[1:div(npoints, nsamp)])
		data[1+(i-1)*nsamp:i*nsamp] = p
	end
	data
end

"compute_noise_autocorr(ljh_filename, nsamp)
computer noise autorcorrelation of length namp from an ljh file of name ljh_filename"
function compute_noise_autocorr(ljh_filename, nsamp)
	data = readcontinuous(LJH.LJHGroup(ljh_filename))
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
	ljh = LJH.LJHGroup(fname)
	traces=[p for (p,t) in ljh[1:end]]
	pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv, 
		peak_index, peak_value, min_value=compute_summary(traces,pretrig_nsamples,LJH.frametime(ljh))


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

push!(perpulse_symbols, :filt_value, :selection_good, :pulse, :rowcount,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good, :filt_phase, :energy)

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
	mc[:filt_phase] = RunningVector(Float32)
	mc[:energy] = RunningVector(Float32)
	mc[:selection_good] = RunningSumBitVector()
	mc[:filt_value_hist] = Histogram(0:1:20000)
	mc[:energy_hist] = Histogram(0:1:20000)
	mc[:pulse] = RunningVector(Vector{UInt16})
	mc[:rowcount] = RunningVector(Int)
	mc[:pretrig_nsamples] = LJH.pretrig_nsamples(ljh)
	mc[:samples_per_record] = LJH.record_nsamples(ljh)
	mc[:frametime] = LJH.frametime(ljh)
	#mc[:rise_time_criteria] = (0,0.0006) # from exafs.basic_cuts, but it is in ms there, s here
	#mc[:pretrig_rms_criteria] = (0.0,30.) # from exafs.basic_cuts
	#mc[:postpeak_deriv_criteria] = (0.0,250.0) # from exafs.basic_cuts
	mc[:f_3db] = 10000
	mc[:known_energies] = [6403.84, 7057.98] # FeKAlpha and FeKBeta

	#metadata
	mc[:workdone_cumulative] = Dict{AbstractStep, Int64}()
	mc[:time_elapsed_cumulative] = Dict{AbstractStep, Float64}()
	mc[:workdone_last] = Dict{AbstractStep, Int64}() 
	mc[:calibration_nextra] = 0 # when finding peaks, how many peaks other than the largest n to include when assigning peaks to energies

	mc[:noise_filename]=noise_filename
	mc[:ljh_filename]=ljh_filename
	mc[:name] = "summarize and filter test"
	mc[:hdf5_filename] = "$(splitext(ljh_filename)[1]).hdf5"
	mc[:oncleanfinish] = markhdf5oncleanfinish
	isfile(mc[:hdf5_filename]) && rm(mc[:hdf5_filename])

	steps = AbstractStep[]
	push!(steps, GetPulsesStep(ljh, [:pulse, :rowcount], 0,100))
	push!(steps, PerPulseStep(compute_summary, [:pulse, :pretrig_nsamples, :frametime],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
	push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :peak_index, :peak_index_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
	push!(steps, ThresholdStep(compute_average_pulse, [:pulse, :selection_good], [:average_pulse], :selection_good, sum, 100, true))
	push!(steps, ThresholdStep(compute_filter, [:average_pulse, :noise_autocorr, :f_3db, :frametime], [:filter, :vdv], :selection_good, sum, 100, true))
	push!(steps, ThresholdStep(compute_noise_autocorr,[:noise_filename, :samples_per_record],[:noise_autocorr], :selection_good, sum, 100, true))
	push!(steps, ThresholdStep(estimate_pretrig_rms_and_postpeak_deriv_criteria,[:noise_filename, :pretrig_nsamples],[:pretrig_rms_criteria, :postpeak_deriv_criteria], :pulse, length, 100, true))
	push!(steps, ThresholdStep(estimate_peak_index_criteria,[:peak_index],[:peak_index_criteria], :rise_time, length, 100, true)) 	
	push!(steps, ThresholdStep(fit_pulse_two_exponential,[:average_pulse, :pretrig_nsamples, :frametime], [:rise_tau_s, :fall_tau_s], :filt_value, length, 10, true))
	push!(steps, PerPulseStep(filter5lag, [:filter, :pulse], [:filt_value, :filt_phase]))
	push!(steps, HistogramStep(update_histogram!, [:filt_value_hist, :selection_good, :filt_value]))
	push!(steps, ThresholdStep(calibrate_nofit, [:filt_value_hist,:known_energies, :calibration_nextra],[:calibration],:filt_value_hist, counted, 1000, true))
	push!(steps, PerPulseStep(apply_calibration, [:calibration, :filt_value], [:energy]) )
	push!(steps, HistogramStep(update_histogram!, [:energy_hist, :selection_good, :energy]))
	push!(steps, ToJLDStep([:filt_value, :filt_phase, :pretrig_rms, :postpeak_deriv, :rise_time, :peak_index, 
	:pretrig_mean, :pulse_average, :pulse_rms, :peak_value, :min_value, :rowcount],
	Pair[:filter=>"filter/filter", :f_3db=>"filter/f_3db", :frametime=>"filter/frametime", :noise_autocorr=>"filter/noise_autocorr", :average_pulse=>"filter/average_pulse", 
	:average_pulse=>"average_pulse",
	:samples_per_record=>"samples_per_record", :frametime=>"frametime", :pretrig_nsamples=>"pretrig_nsamples",
	:ljh_filename=>"ljh_filename", :noise_filename=>"noise_filename"],
	mc[:hdf5_filename]))
	#push!(steps, FreeMemoryStep(graph(steps)))
	#push!(steps, MemoryLimitStep(Int(4e6))) # throw error if mc uses more than 4 MB
	# write a verification function that makes sure all inputs either exist, or are the output of another step
	mc[:steps]=steps

	make_task(mc)
	mc
end

function MASS_MATTER_watcher(masschannels, exitchannel)
	ljhname, writingbool = "",false
	last_noise_filename = ""
	analyzing_fname = "__"
	while true
		watch_file(LJHUtil.sentinel_file_path,10) #blocks task until file changes
		oldljhname, oldwritingbool = ljhname, writingbool
		ljhname, writingbool = LJHUtil.matter_writing_status()
		if ljhname != oldljhname
			if oldljhname == analyzing_fname # stop tasks for previous file once they've finished their work
				info("allowing analysis tasks for $oldljhname to end")
				map(plan_to_end, values(masschannels))
				analyzing_fname=""
			end
			if contains(ljhname, "noise") || contains(ljhname,".noi")
				last_noise_filename = ljhname
				info("New noise file $last_noise_filename")
			else
				last_noise_filename == "" && error("must set last_noise_filename before processing data")
				analyzing_fname = ljhname
				info("Starting analysis of $ljhname with noise from $last_noise_filename")
				t0 = time()
				channums = LJHUtil.ljhallchannels(ljhname)
				channums = channums[1:min(240, length(channums))]
				ljh_filenames = [LJHUtil.ljhfnames(ljhname,channum) for channum in channums]
				noise_filenames = [LJHUtil.ljhfnames(last_noise_filename,channum) for channum in channums]
				for i in eachindex(channums)
					masschannels[channums[i]] = setup_channel(ljh_filenames[i], noise_filenames[i])
				end
				tf = time()
				info("setup $(length(channums)) channels in $(tf-t0) seconds")
				tschedule = @elapsed map(schedule, values(masschannels))
				info("launched $(length(channums)) channels in $tschedule seconds")
			end
		end
		isready(exitchannel) && return
	end
end

function schedule_MASS_MATTER_watcher()
	masschannels = Dict()
	exitchannel = Channel{Int}(1)
	task = @schedule MASS_MATTER_watcher(masschannels, exitchannel)
	masschannels, exitchannel, task
end

masschannels, watcher_exitchannel, watcher_task = schedule_MASS_MATTER_watcher()

@schedule begin
getopenfilelimit() = parse(Int,split(split(readall(`ulimit -a`),"\n")[6])[end])
if getopenfilelimit()>=1000
	sleep(2)
	LJHUtil.write_sentinel_file("/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes_noise/20150730_50mM_irontris_100ps_xes_noise.ljh",false)
	sleep(2)
	LJHUtil.write_sentinel_file("/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes/20150730_50mM_irontris_100ps_xes.ljh",false)
	sleep(2)
	LJHUtil.write_sentinel_file("/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes_noise/20150730_50mM_irontris_100ps_xes_noise.ljh",false)
	sleep(2)
	put!(watcher_exitchannel,1) 
	@schedule begin
		tasks = [mc[:task] for (ch, mc) in masschannels]
		for task in tasks # wait for all tasks to finish
			try wait(task) end # wait rethrows errors from failed tasks, we can always look at them later
		end
		ndone = sum([task.state==:done for task in tasks])
		nfailed = sum([task.state==:failed for task in tasks])
		info("all tasks done!! Yay. $ndone :done and $nfailed :failed")
	end;
else
	println("open file limit is too low, try ulimit -n 1000")
end
end # write to MATTER sentinel file to simulate matter writing various files



# publish hists over ZMQ
function make_coadded_hist(masschannels)
	coadded_hist = Histogram(0:1:20000)
	n = 0
	for (ch,mc) in masschannels
		if haskey(mc,:energy_hist)
			coadded_hist+=mc[:energy_hist]
			n+=1
		end
	end
	coadded_hist,n
end

function publish_hists_zmq(masschannels, zmqexitchannel)
	context=Context()
	pub=Socket(context, PUB)
	ZMQ.bind(pub, "tcp://*:5555")
	try
		while !isready(zmqexitchannel)
			coadded_hist,n = make_coadded_hist(masschannels)
			io = IOBuffer()
			write(io,"0")
			serialize(io, (coadded_hist,n))
			ZMQ.send(pub, Message(io))
			sleep(1)
		end
	finally
		ZMQ.close(pub)
		ZMQ.close(context)
	end
end

zmqexitchannel = Channel()
@schedule publish_hists_zmq(masschannels, zmqexitchannel)




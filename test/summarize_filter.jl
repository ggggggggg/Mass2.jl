using Mass2
import TESOptimalFiltering: filter5lag, calculate_filter, autocorrelation
using ReferenceMicrocalFiles

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
	data = zeros(Uint16,npoints)
	for (i,(p,t)) in enumerate(ljh[1:div(npoints, nsamp)])
		data[1+(i-1)*nsamp:i*nsamp] = p
	end
	data
end

"compute_noise_autocorr(ljh_filename, nsamp)
computer noise autorcorrelation of length namp from an ljh file of name ljh_filename"
function compute_noise_autocorr(ljh_filename, nsamp)
	data = readcontinuous(LJH.LJHGroup(ljh_filename))
	autocorrelation(data, nsamp, 600) # final argument is max_excursion
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

push!(perpulse_symbols, :filt_value, :selection_good, :pulse, :rowcount,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good, :filt_phase)

function setup_channel(ljh_filename, noise_filename)
	ljh = LJHGroup(ljh_filename)



	c=MassChannel()
	c[:pretrig_mean] = RunningVector(Float32)
	c[:pretrig_rms] = RunningVector(Float32)
	c[:pulse_average] = RunningVector(Float32)
	c[:pulse_rms] = RunningVector(Float32)
	c[:rise_time] = RunningVector(Float32)
	c[:postpeak_deriv] = RunningVector(Float32)
	c[:peak_index] = RunningVector(Uint16)
	c[:peak_value] = RunningVector(Uint16)
	c[:min_value] = RunningVector(Uint16)
	c[:filt_value] = RunningVector(Float32)
	c[:filt_phase] = RunningVector(Float32)
	c[:selection_good] = RunningSumBitVector()
	c[:filt_value_hist] = Histogram(0:1:20000)
	c[:pulse] = RunningVector(Vector{UInt16})
	c[:rowcount] = RunningVector(Int)
	c[:pre_samples] = LJH.pretrig_nsamples(ljh)
	c[:samples_per_record] = LJH.record_nsamples(ljh)
	c[:frame_time] = LJH.frametime(ljh)
	c[:rise_time_criteria] = (0,0.0006) # from exafs.basic_cuts, but it is in ms there, s here
	c[:pretrig_rms_criteria] = (0.0,30.) # from exafs.basic_cuts
	c[:postpeak_deriv_criteria] = (0.0,250.0) # from exafs.basic_cuts
	c[:f_3db] = 10000

	#metadata
	c[:workdone_cumulative] = Dict{AbstractStep, Int64}()
	c[:stepelapsed_cumulative] = Dict{AbstractStep, Float64}()
	#c[:workdone_last] = Dict{AbstractStep, Int64}([s=>-1 for s in steps]) # start with nonzero values so autoender won't end before starting
	c[:workdone_last] = Dict{AbstractStep, Int64}() 

	c[:noise_filename]=noise_filename
	c[:ljh_filename]=ljh_filename
	c[:name] = "summarize and filter test"
	c[:hdf5_filename] = "$(splitext(ljh_filename)[1]).hdf5"
	c[:trick_mass_str]="this is exists only so that the group filters will exist so mass will skip work in compute_filters"
	isfile(c[:hdf5_filename]) && rm(c[:hdf5_filename])

	steps = AbstractStep[]
	push!(steps, GetPulsesStep(ljh, [:pulse, :rowcount], 0,100))
	push!(steps, PerPulseStep(compute_summary, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
	push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :rise_time, :rise_time_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
	push!(steps, ThresholdStep(compute_average_pulse, [:pulse, :selection_good], [:average_pulse], :selection_good, sum, 100, true))
	push!(steps, ThresholdStep(compute_filter, [:average_pulse, :noise_autocorr, :f_3db, :frame_time], [:filter, :vdv], :selection_good, sum, 100, true))
	push!(steps, ThresholdStep(compute_noise_autocorr,[:noise_filename, :samples_per_record],[:noise_autocorr], :selection_good, sum, 100, true))
	push!(steps, PerPulseStep(filter5lag, [:filter, :pulse], [:filt_value, :filt_phase]))
	push!(steps, HistogramStep(update_histogram!, [:filt_value_hist, :selection_good, :filt_value]))
	push!(steps, ToJLDStep([:filt_value, :filt_phase, :pretrig_rms, :postpeak_deriv, :rise_time, :peak_index, 
	:pretrig_mean, :pulse_average, :pulse_rms, :peak_value, :min_value, :rowcount],
	Pair[:filter=>"filter/filter", :f_3db=>"filter/f_3db", :frame_time=>"filter/frametime", :noise_autocorr=>"filter/noise_autocorr", :average_pulse=>"filter/average_pulse", 
	:average_pulse=>"average_pulse",
	:samples_per_record=>"samples_per_record", :frame_time=>"frametime", :pre_samples=>"pre_samples",
	:ljh_filename=>"ljh_filename", :noise_filename=>"noise_filename"],
	c[:hdf5_filename]))
	push!(steps, FreeMemoryStep(graph(steps)))
	push!(steps, MemoryLimitStep(Int(4e6))) # throw error if c uses more than 4 MB
	# write a verification function that makes sure all inputs either exist, or are the output of another step
	c[:steps]=steps


	c
end

function autoender(c::MassChannel, tasks, exitchannels)
	while true
		if all(collect(values(c[:workdone_last])).==0) || any([istaskdone(t) for t in values(tasks)])
			for exitchannel in values(exitchannels)
				!isready(exitchannel) && put!(exitchannel,1)
			end
			return
		end
	yield()
	end
end


function launch_channel(ljh_filename, noise_filename)
	c = setup_channel(ljh_filename, noise_filename)
	exitchannel = Channel{Int}(1)
	c[:exitchannels] = Dict(1=>exitchannel)
	c[:workdone_last] = Dict()
	t=@task begin   while !isready(exitchannel)
						for s in c[:steps]
							yield()
							time_elapsed = @elapsed workdone = workunits(dostep!(s,c))
							workdone_cumulative = get(c[:workdone_cumulative],s,0)
							c[:workdone_cumulative][s] = workdone_cumulative+workdone
							time_elapsed_cumulative = get(c[:stepelapsed_cumulative],s,0)
							c[:stepelapsed_cumulative][s] = time_elapsed_cumulative+time_elapsed	
							c[:workdone_last][s] = workdone
						end
					end
					h5open(c[:hdf5_filename],"r+") do h5
					h5["clean_exit_posix_s"]=time()
					end # note when/if analysis was finished without error, will be used to tag "good" hdf5 files
			end		
	c[:tasks] = Dict(1=>t) # only one task, it's in a dict for compatability with having a task per step
	c[:wait_tasks] = [@schedule(try wait(t) end) for t in values(c[:tasks])] # supress printing of error messages
	c[:endertask] = @task autoender(c,c[:tasks],c[:exitchannels])
	c
end

function automass(masschannels, exitchannel)
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
				for (channum, mc) in masschannels
					schedule(mc[:endertask]) # when the ljh file name changes we wait until all work is done, then close the tasks processing the old file
				end
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
				channums = channums[1:min(1, length(channums))]
				ljh_filenames = [LJHUtil.ljhfnames(ljhname,channum) for channum in channums]
				noise_filenames = [LJHUtil.ljhfnames(last_noise_filename,channum) for channum in channums]
				for i in eachindex(channums)
					masschannels[channums[i]] = launch_channel(ljh_filenames[i], noise_filenames[i])
				end
				tsetup = time()
				info("setup $(length(channums)) channels in $(tsetup-t0) seconds")
				for masschannel in values(masschannels)
					schedule(masschannel[:tasks][1])
					# dont print or do io in this loop, print gives up control to schedueler
					# basically dont do anything else in this loop
				end
				tf = time()
				info("launched $(length(channums)) channels in $(tf-tsetup) seconds")
			end
		end
		isready(exitchannel) && return
	end
end

function launch_automass()
	masschannels = Dict()
	exitchannel = Channel{Int}(1)
	task = @schedule automass(masschannels, exitchannel)
	masschannels, exitchannel, task
end

# savegraph("graph",graph(c[:steps]))

getopenfilelimit() = parse(Int,split(split(readall(`ulimit -a`),"\n")[6])[end])
if getopenfilelimit()>=1000
	masschannels, exitchannel, task = launch_automass()
	sleep(1)
	LJHUtil.write_sentinel_file("/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes_noise/20150730_50mM_irontris_100ps_xes_noise.ljh",false)
	sleep(1)
	LJHUtil.write_sentinel_file("/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes/20150730_50mM_irontris_100ps_xes.ljh",false)
	sleep(1)
	LJHUtil.write_sentinel_file("/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes_noise/20150730_50mM_irontris_100ps_xes_noise.ljh",false)
else
	println("open file limit is too low, try ulimit -n 1000")
end

@schedule begin
	tasks = [mc[:tasks][1] for (ch, mc) in masschannels]
	for task in tasks # wait for all tasks to finish
		try wait(task) end # wait rethrows errors from failed tasks, we can always look at them later
	end
	ndone = sum([task.state==:done for task in tasks])
	nfailed = sum([task.state==:failed for task in tasks])
	info("all tasks done!! Yay. $ndone :done and $nfailed :failed")
end;
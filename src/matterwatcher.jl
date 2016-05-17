"Use schedule_MATTER_watcher to launch in a task"
function MATTER_watcher(masschannels, exitchannel, setup_channel_func, maxchannels)
	println("starting MATTER_watcher")	
	ljhname, writingbool = "",false
	last_noise_filename = ""
	last_started_analyzing_fname = "" # not neccesary, but I could see using the info in the future
	have_planned_to_end = true
	i=0
	while true
		watch_file(LJHUtil.sentinel_file_path,10) #blocks task until file changes, or timeout given in seconds
		oldljhname, oldwritingbool = ljhname, writingbool
		ljhname, writingbool = LJHUtil.matter_writing_status()
		if (ljhname != oldljhname)
			if contains(ljhname, "noise") || contains(ljhname,".noi")
				last_noise_filename = ljhname
				println("New noise file $last_noise_filename")
			else
				println("new ljh file observed: $ljhname")
				if last_noise_filename == ""
					println("must set last_noise_filename before processing data")
					continue
				end
				last_started_analyzing_fname = ljhname
				println("Starting analysis of $ljhname with noise from $last_noise_filename")
				t0 = time()
				channums = LJHUtil.allchannels(ljhname)
				channums = channums[1:min(maxchannels,length(channums))]
				ljh_filenames = [LJHUtil.fnames(ljhname,channum) for channum in channums]
				noise_filenames = [LJHUtil.fnames(last_noise_filename,channum) for channum in channums]
				for i in eachindex(channums)
					masschannels[channums[i]] = setup_channel_func(ljh_filenames[i], noise_filenames[i])
				end
				tf = time()
				println("setup $(length(channums)) channels in $(tf-t0) seconds")
				tschedule = @elapsed map(schedule, values(masschannels))
				have_planned_to_end = false
				println("launched $(length(channums)) channels in $tschedule seconds")
			end
		end
		if !writingbool && !have_planned_to_end
			# stop tasks for previous file once they've finished their work
			println("allowing analysis tasks for $oldljhname to end")
			map(plantoend, values(masschannels))
			have_planned_to_end=true
		end
		if isready(exitchannel)
			println("exiting MATTER_watcher, due to ready exitchannel")			 
			return
		end
	end
end

"masschannels, watcher_exitchannel, watcher_task = schedule_MATTER_watcher(setup_channel_func, maxchannels=1000)
setup_channel_func(ljh_filename, noise_filename) should be a function that returns a MassChannel designed to analyze those files
maschannels is an `Int` that limits the maximum number of channels to analyze
Launches MATTER_watcher in a task. It will monitor sentinel_file_path for changes. It stores the name 
of files that contain the text `noise` or `.noi` as noise files, then starts analysis on any other files, using
the last seen noise file as the noise. When writing stops, or when a new file starts being written it 
plans to end the analysis steps. This means that when all work has finished, the analysis tasks will end.
masschannels is a dict channel_number->masschannel object for the currently analyzing files
watcher_exitchannel will cause the MATTER_watcher task to end upon put!(watcher_exitchannel,1) (or any other int)
watcher_task is the task executing MATTER_watcher"
function schedule_MATTER_watcher(setup_channel_func, maxchannels=1000)
	masschannels = Dict()
	exitchannel = Channel{Int}(1)
	task = @schedule MATTER_watcher(masschannels, exitchannel, setup_channel_func, maxchannels)
	masschannels, exitchannel, task
end

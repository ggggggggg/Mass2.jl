using Mass2
import TESOptimalFiltering: filter5lag, calculate_filter, autocorrelation
using ReferenceMicrocalFiles


function selectfromcriteria(x...) # placeholder, kinda ugly to use and probalby a bit slow
	iseven(length(x)) || error("x must be indicator,criteria,indicator,criteria...")
	out = trues(length(x[1]))
	for i = 1:2:length(x)
		low, high = x[i+1]
		out &= low .< x[i] .< high
	end
	out
end

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

function compute_noise_autocorr(ljh_filename, nsamp)
	data = readcontinuous(LJH.LJHGroup(ljh_filename))
	autocorrelation(data, nsamp, 600) # final argument is max_excursion
end

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

function compute_filter(average_pulse, noise_autocorr, f_3db, dt)
	filter, variance = calculate_filter(average_pulse, noise_autocorr, f_3db, dt)
	filter*=maximum(average_pulse)-minimum(average_pulse) # like mass, normalize so you get pulseheight similar to raw units
	vdv = 1/sqrt(variance)/2.355
	filter,vdv
end

push!(perpulse_symbols, :filt_value, :selection_good, :pulse, :rowstamp,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good)

function setup_channel(ljh_filename, noise_filename)
	ljh = LJHGroup(ljh_filename)

	steps = AbstractStep[]
	push!(steps, GetPulsesStep(ljh, [:pulse, :rowstamp], 0,500))
	push!(steps, PerPulseStep(compute_summary, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
	push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :rise_time, :rise_time_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
	push!(steps, ThresholdStep(compute_average_pulse, [:pulse, :selection_good], [:average_pulse], :selection_good, sum, 100, true))
	push!(steps, ThresholdStep(compute_filter, [:average_pulse, :noise_autocorr, :f_3db, :frame_time], [:filter, :vdv], :selection_good, sum, 100, true))
	push!(steps, PerPulseStep(filter5lag, [:filter, :pulse], [:filt_value, :filt_phase]))
	push!(steps, HistogramStep(update_histogram!, [:filt_value_hist, :selection_good, :filt_value]))
	push!(steps, ToJLDStep([:filt_value,:pretrig_rms],"$(splitext(basename(ljh_filename))[1]).hdf5"))
	push!(steps, FreeMemoryStep(graph(steps)))

	c=MassChannel()
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
	c[:filt_phase] = RunningVector(Float64)
	c[:selection_good] = RunningSumBitVector()
	c[:filt_value_hist] = Histogram(0:1:20000)
	c[:pulse] = RunningVector(Vector{Int})
	c[:rowstamp] = RunningVector(Int)
	c[:pre_samples] = LJH.pretrig_nsamples(ljh)
	c[:nsamp] = LJH.record_nsamples(ljh)
	c[:frame_time] = LJH.frametime(ljh)
	c[:rise_time_criteria] = (0,0.0006) # from exafs.basic_cuts, but it is in ms there, s here
	c[:pretrig_rms_criteria] = (0.0,30.) # from exafs.basic_cuts
	c[:postpeak_deriv_criteria] = (0.0,250.0) # from exafs.basic_cuts
	c[:noise_autocorr] = compute_noise_autocorr(noise_filename,c[:nsamp])
	c[:f_3db] = 10000

	#metadata
	c[:steps]=steps
	c[:workdone_cumulative] = Dict{AbstractStep, Int64}()
	c[:stepelapsed_cumulative] = Dict{AbstractStep, Float64}()
	c[:workdone_last] = Dict{AbstractStep, Int64}()
	c[:name] = "summarize and filter test"

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

function repeatstep(c::MassChannel, s::AbstractStep, exitchannel::Channel{Int})
	workdone_cumulative=0
	time_elapsed_cumulative=0
	while !isready(exitchannel)
		# do the step, record workdone and time elapsed
		time_elapsed = @elapsed workdone= workunits(dostep!(s,c))
		workdone_cumulative+=workdone
		time_elapsed_cumulative+=time_elapsed
		c[:workdone_cumulative][s] = workdone_cumulative
		c[:stepelapsed_cumulative][s] = time_elapsed_cumulative
		c[:workdone_last][s] = workdone
		yield()
	end
end

function launch_channel(ljh_filename, noise_filename)
	tasks = Dict()
	exitchannels = Dict()
	for s in c[:steps]
		exitchannel = Channel{Int}(1)
		tasks[s] = @schedule repeatstep(c, s, exitchannel)
		exitchannels[s] = exitchannel
	end
	@schedule autoender(c,tasks, exitchannels)
	tasks,exitchannels
end

ljh_filename = "/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes/20150730_50mM_irontris_100ps_xes_chan479.ljh"
noise_filename = "/Volumes/Drobo/exafs_data/20150730_50mM_irontris_100ps_xes_noise/20150730_50mM_irontris_100ps_xes_noise_chan479.ljh"
c = setup_channel(ljh_filename, noise_filename)
savegraph("graph",graph(c[:steps]))

tasks,exitchannels = launch_channel(ljh_filename, noise_filename);
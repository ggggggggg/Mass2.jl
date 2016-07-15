using Mass2
using ReferenceMicrocalFiles
using ZMQ
using Base.Test

"CountRateHistogram contains a RunningVector of Histograms, where hists[i] contains
all counts observed with `first(posix_timestamp_usec_bin_edges)<timestamp<=posix_timestamp_usec_bin_edges[i]`.
Count rates over ROIS can easily be calculated by differences of histograms."
type CountRateHistogram <: AbstractHistogram
  posix_timestamp_usec_bin_edges::Range{Int}
  hists::RunningVector{Histogram}
  x::FloatRange{Float64} # bins of the histograms in hists
  countedname::ASCIIString
end

"`posix_timestamp_usec_bin_edges` should always start at 0 (this should be an inner constructor?)"
function CountRateHistogram(posix_timestamp_usec_bin_step,x, countedname)
  posix_timestamp_usec_bin_edges = 0:posix_timestamp_usec_bin_step:typemax(Int64)
  h = CountRateHistogram(posix_timestamp_usec_bin_edges, RunningVector(Vector{Histogram}(), 0), x, countedname)
end
"update_histogram!(h::CountRateHistogram, selection, events, timestamps) adds the values in `events[selection]`
to the `h`, adding to the appropriate entry in `hists` based on `timestamps`. returns nothing.
timestamps must be sorted!!"
function update_histogram!(h::CountRateHistogram, selection, events, timestamps, last_timestamp_checked)
  @assert(length(selection)==length(events)==length(timestamps))
  if length(selection) == 0 # if the input arrays are empty, return early
    return donethru(h)
  end
  #last_timestamp_checked = max(last_timestamp_checked, last(timestamps)
  last_timestamp_checked = last(timestamps)
  j_hi = 0
  while true
    if length(available(h.hists)) >=1 # figure out what timestamp range we're interested in
      i_lo, i_hi = length(h.hists), length(h.hists)+1 # indexes into h.posix_timestamp_usec_bin_edges
    else
      i_range = searchsorted(h.posix_timestamp_usec_bin_edges, first(timestamps))
      i_lo, i_hi = last(i_range), first(i_range)
      i_lo = i_lo == i_hi ? i_lo-1 : i_lo # if they are equal, it means the fist timestamp is exactly equal to h.posix_timestamp_usec_bin_edges[i_hi]
      freeuntil!(h.hists, i_lo)
    end
    i_lo, i_hi
    t_lo, t_hi = h.posix_timestamp_usec_bin_edges[i_lo], h.posix_timestamp_usec_bin_edges[i_hi]
    if last_timestamp_checked<t_hi
      return donethru(h)
    end # only add histogram for which we have seen all events
    # t_lo and t_hi are the edges of the time bin to be used for hists[i_hi]
    j_lo = j_hi+1
    while j_hi+1 <= length(timestamps) && timestamps[j_hi+1]<=t_hi
      j_hi+=1
    end
    if j_lo<=j_hi
      @assert(timestamps[j_lo]>t_lo)
    else
      return donethru(h)
    end
    newhist = Histogram(h.x)
    update_histogram!(newhist, selection[j_lo:j_hi], events[j_lo:j_hi])
    if i_lo in available(h.hists)
      newhist = newhist+h.hists[i_lo]
    end
    push!(h.hists, newhist)
  end
  return donethru(h)
end

function donethru(h::CountRateHistogram)
  if length(available(h.hists)) == 0
    0
  else
    h.hists[end].seen
  end
end

tstep = 1000000 # usec
Nstep = 10
Npoints = 1000000
tfirst = 8000000000
crh = CountRateHistogram(tstep, 0:1:1000, "test")


events = ones(Npoints);
selection = trues(Npoints);
tlast = tfirst+Nstep*tstep
timestamps = round(Int,linspace(tfirst+1, tlast, Npoints))

update_histogram!(crh, selection, events, timestamps, last(timestamps))
donethru(crh) == Npoints
length(available(crh.hists)) == Nstep
for (i,j) in enumerate(available(crh.hists))
  @assert counts(crh.hists[j])[1]==div(Npoints,Nstep)*i
end

using ZMQ
port = 5555
context = ZMQ.Context()
socket = ZMQ.Socket(context,ZMQ.PUB)
ZMQ.bind(socket, "tcp://*:$port")
function sendtest()
  ZMQ.send(socket,"test")
end
"(channel_number::Int, crh, i) send CountRateHistogram `crh` over a zmq socket with format:
channel number; bin edges; counts; total counts seen; timestamps of start and end of new observations"
function send_crh(channel_number::Integer, crh, i)
  h = crh.hists[i]
  # last argument is SNDMORE, if true, it builds up a multipart message
  ZMQ.send(socket, crh.countedname, true)
  ZMQ.send(socket, "$channel_number", true)
  ZMQ.send(socket, repr(bin_edges(h)), true)
  ZMQ.send(socket, repr(counts(h)), true)
  ZMQ.send(socket, repr(seen(h)), true)
  ZMQ.send(socket, repr(crh.posix_timestamp_usec_bin_edges[[i-1,i]]))
end
function send_and_free_crh(channel_number::Integer, crh, i)
  send_crh(channel_number, crh, i)
  freeuntil!(crh.hists,i)
end
sub_socket = ZMQ.Socket(context, ZMQ.SUB)
ZMQ.subscribe(sub_socket)
ZMQ.connect(sub_socket, "tcp://localhost:5555")


function recv_multipart(socket)
  out = []
  push!(out, recv(socket))
  while ZMQ.ismore(socket)
    push!(out, recv(socket))
  end
  out
end

# @schedule for i = 1:100000 println(i);sleep(0.001); println(map(bytestring,recv_multipart(sub_socket))) end


# sendtest()
# send_and_free_crh(1, crh, 8011)

function update_histogram_send!(h::CountRateHistogram, selection, events, timestamps, last_timestamp_checked, channel_number)
  update_histogram!(h::CountRateHistogram, selection, events, timestamps, last_timestamp_checked)
  if length(available(h.hists))>=2
    i = available(h.hists)[end-1]
    send_and_free_crh(channel_number, h, i)
  end
  return nothing
end

"setup_channel(ljh_filename, noise_filename) should return a `MassChannel`"
function setup_channel(ljh_filename, noise_filename)
  # TODO:
  # pass in HDF5 group, need to close it on exit
  # pass in "init dictionary"?

    mc=MassChannel()
    mc[:ljh] = LJHGroup(ljh_filename)
    mc[:pretrig_nsamples] = LJH.pretrig_nsamples(mc[:ljh])
    mc[:samples_per_record] = LJH.record_nsamples(mc[:ljh])
    mc[:frametime] = LJH.frametime(mc[:ljh])
    mc[:channel_number] = LJH.channel(mc[:ljh])
    mc[:last_timestamp_checked] = 0.0 # this is updated via a total HACK, its hardcoded into GetPulsesStep

    # look for an init_hdf5 file, extract some values into init_dict
    # if init_dict is empty, get(mc.task) should show fail with a somewhat descriptive error
    init_dict = get_init_dict(mc[:channel_number])
    if isempty(init_dict)
	return mc
    end
    for (k,v) in init_dict
	mc[k]=v
    end
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
    mc[:energy_from_pulse_rms] = RunningVector(Float32)
    mc[:pulse] = RunningVector(Vector{UInt16})
    mc[:rowcount] = RunningVector(Int)
    mc[:timestamp_posix_usec] = RunningVector(Int)
    mc[:selection_good] = RunningSumBitVector()
    mc[:filt_value_hist] = CountRateHistogram(10_000_000, 0:0.5:1500,"filt_value")
    mc[:energy_hist] = CountRateHistogram(10_000_000, 0:0.5:1500,"energy_from_filt_value")
    mc[:energy_from_pulse_rms_hist] = CountRateHistogram(10_000_000, 0:0.5:1500,"energy_from_pulse_rms")

    mc[:noise_filename]=ascii(noise_filename)
    mc[:ljh_filename]=ascii(ljh_filename) # h5py can't read UFT8String written by julia
    mc[:name] = "summarize and filter test"
    mc[:hdf5_filename] = "$(splitext(ljh_filename)[1])_jl.hdf5"
    mc[:init_hdf5_filename] = "$(splitext(ljh_filename)[1])_init_jl.hdf5"

    mc[:oncleanfinish] = markhdf5oncleanfinish
    isfile(mc[:hdf5_filename]) && rm(mc[:hdf5_filename])

    steps = AbstractStep[
    GetPulsesStep(mc[:ljh], [:pulse, :rowcount, :timestamp_posix_usec], 0,100, 0)
    @perpulse pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv, peak_index, peak_value, min_value = compute_summary(pulse, pretrig_nsamples, frametime)
    @threshold calibration_pulse_rms = calibrate_xy(calibration_pulse_rms_pulseheights, calibration_pulse_rms_energies) when length(peak_index)>0
    @threshold calibration_filt_value = calibrate_xy(calibration_filt_value_dc_pulseheights, calibration_filt_value_dc_energies) when length(peak_index)>0 # calibration from filt value dc should work fine for filt value, probably still change this in future
    @perpulse selection_good = selectfromcriteria(pretrig_rms, pretrig_rms_criteria, peak_index, peak_index_criteria, postpeak_deriv, postpeak_deriv_criteria)
    @perpulse filt_value, filt_phase = filter5lag(filter, pulse)
    @perpulse energy = apply_calibration(calibration_filt_value, filt_value)
    @perpulse energy_from_pulse_rms = apply_calibration(calibration_pulse_rms, pulse_rms)
    @histogram update_histogram_send!(filt_value_hist, selection_good, filt_value, timestamp_posix_usec, last_timestamp_checked, channel_number)
    @histogram update_histogram_send!(energy_hist, selection_good, energy, timestamp_posix_usec, last_timestamp_checked, channel_number)
    ToJLDStep([:filt_value, :energy, :filt_phase, :pretrig_rms, :postpeak_deriv, :rise_time, :peak_index,
    :pretrig_mean, :pulse_average, :pulse_rms, :peak_value, :min_value, :rowcount],
    Pair[:filter=>"filter/filter", :f_3db=>"filter/f_3db", :frametime=>"filter/frametime", :noise_autocorr=>"filter/noise_autocorr", :average_pulse=>"filter/average_pulse",
    :average_pulse=>"average_pulse",
    :samples_per_record=>"samples_per_record", :frametime=>"frametime", :pretrig_nsamples=>"pretrig_nsamples",
    :ljh_filename=>"ljh_filename", :noise_filename=>"noise_filename",
  	:peak_index_criteria=>"selection_criteria/peak_index", :pretrig_rms_criteria=>"selection_criteria/pretrig_rms", :postpeak_deriv_criteria=>"selection_criteria/postpeak_deriv",
  	:pretrigger_mean_drift_correction_params=>"pretrigger_mean_drift_correction_params"],
    mc[:hdf5_filename])
    MemoryLimitStep(Int(4e6)) # throw error if mc uses more than 4 MB
    FreeMemoryStep()
    ] # end steps
    setsteps!(mc, steps)
    preptasks!(mc)
    mc
end
get_init_hdf5_filename_groupname(channel_number) = "ssrl_realtime_init_20160602.hdf5","chan$channel_number"
function get_init_dict(channel_number::Integer)
	init_hdf5_filename, groupname = get_init_hdf5_filename_groupname(channel_number)
	d = Dict{Symbol, Any}()
	h5open(init_hdf5_filename,"r") do h5
	    if !exists(h5,groupname)
		return nothing
	    end
	    grp = h5[groupname]
	    d[:calibration_pulse_rms_energies] = read(grp["calibrations/calibration_pulse_rms_energies"])
	    d[:calibration_pulse_rms_pulseheights] = read(grp["calibrations/calibration_pulse_rms_pulseheights"])
	    d[:calibration_filt_value_dc_energies] = read(grp["calibrations/calibration_filt_value_dc_energies"])
	    d[:calibration_filt_value_dc_pulseheights] = read(grp["calibrations/calibration_filt_value_dc_pulseheights"])
	    d[:filter] = read(grp["filter/filter"])
	    d[:pretrigger_mean_drift_correction_params] = LinearDriftCorrect(read(grp["pretrigger_mean_drift_correction_params"])...)
	    d[:peak_index_criteria] = read(grp["selection_criteria"]["peak_index"])
	    d[:pretrig_rms_criteria] = read(grp["selection_criteria"]["pretrig_rms"])
	    d[:postpeak_deriv_criteria] = read(grp["selection_criteria"]["postpeak_deriv"])
	    return nothing
	end # do
	return d
end # get_init_dict



if true # change this to false if you want to type in the filename
  ljhname, writingbool = LJHUtil.matter_writing_status()
else
  ljhname = "put filename here"
  writingbool = false
end
println("***********************")
@show ljhname
@show maxchannels = 240
masschannels = Dict()
t0 = time()
channums = LJHUtil.allchannels(ljhname)
channums = channums[1:min(maxchannels,length(channums))]
ljh_filenames = [LJHUtil.fnames(ljhname,channum) for channum in channums]
noise_filenames = ljh_filenames
for (i, channum) in enumerate(channums)
  masschannels[channums[i]] = setup_channel(ljh_filenames[i], noise_filenames[i])
end
tf = time()
println("setup $(length(channums)) channels in $(tf-t0) seconds")
tschedule = @elapsed map(schedule, values(masschannels))
println("launched $(length(channums)) channels in $tschedule seconds")
# stop tasks for previous file once they've finished their work
println("allowing analysis tasks for $oldljhname to end")

while writingbool
  watch(watch_file(LJHUtil.sentinel_file_path,10))
  ljhname, writingbool = LJHUtil.matter_writing_status()
end

map(plantoend, values(masschannels)) # tell each channel it can end

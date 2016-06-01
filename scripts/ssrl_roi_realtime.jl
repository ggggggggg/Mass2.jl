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
end

"`posix_timestamp_usec_bin_edges` should always start at 0 (this should be an inner constructor?)"
function CountRateHistogram(posix_timestamp_usec_bin_step,x)
  posix_timestamp_usec_bin_edges = 0:posix_timestamp_usec_bin_step:typemax(Int64)
  h = CountRateHistogram(posix_timestamp_usec_bin_edges, RunningVector(Vector{Histogram}(), 0), x)
  h
end
"update_histogram!(h::CountRateHistogram, selection, events, timestamps) adds the values in `events[selection]`
to the `h`, adding to the appropriate entry in `hists` based on `timestamps`. returns nothing.
timestamps must be sorted!!"
function update_histogram!(h::CountRateHistogram, selection, events, timestamps, last_timestamp_checked)
  @assert(length(selection)==length(events)==length(timestamps))
  if length(selection) == 0 # if the input arrays are empty, return early
    return donethru(h)
  end
  last_timestamp_checked = max(last_timestamp_checked, last(timestamps))

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
  if length(h.hists) == 0
    0
  else
    h.hists[end].seen
  end
end

tstep = 1000000 # usec
Nstep = 10
Npoints = 1000000
tfirst = 8000000000
crh = CountRateHistogram(tstep, 0:1:1000)


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
    mc[:energy_from_pulse_rms] = RunningVector(Float32)
    mc[:pulse] = RunningVector(Vector{UInt16})
    mc[:rowcount] = RunningVector(Int)
    mc[:timestamp_posix_usec] = RunningVector(Int)
    mc[:selection_good] = RunningSumBitVector()
    mc[:filt_value_hist] = CountRateHistogram(10_000_000, 0:1:2000)
    mc[:energy_hist] = CountRateHistogram(10_000_000, 0:1:2000)
    mc[:energy_from_pulse_rms_hist] = CountRateHistogram(10_000_000, 0:1:2000)

    mc[:pretrig_nsamples] = LJH.pretrig_nsamples(ljh)
    mc[:samples_per_record] = LJH.record_nsamples(ljh)
    mc[:frametime] = LJH.frametime(ljh)
    mc[:channel_number] = LJH.channel(ljh)
    mc[:last_timestamp_checked] = 0.0 # this is updated via a total HACK

    #metadata
    mc[:calibration_nextra] = 1 # when finding peaks, how many peaks other than the largest n to include when assigning peaks to energies

    mc[:noise_filename]=ascii(noise_filename)
    mc[:ljh_filename]=ascii(ljh_filename) # h5py can't read UFT8String written by julia
    mc[:name] = "summarize and filter test"
    mc[:hdf5_filename] = "$(splitext(ljh_filename)[1])_jl.hdf5"
    mc[:init_hdf5_filename] = "$(splitext(ljh_filename)[1])_init_jl.hdf5"

    mc[:oncleanfinish] = markhdf5oncleanfinish
    isfile(mc[:hdf5_filename]) && rm(mc[:hdf5_filename])

    steps = AbstractStep[
    FromHDF5Step([], Pair[:filter=>"filter/filter", :calibration_pulse_rms_energies=>"calibrations/calibration_pulse_rms_energies",
    :calibration_pulse_rms_pulseheights=>"calibrations/calibration_pulse_rms_pulseheights",
     :calibration_filt_value_energies=>"calibrations/calibration_filt_value_energies",
    :calibration_filt_value_pulseheights=>"calibrations/calibration_filt_value_pulseheights"],
    mc[:init_hdf5_filename])
    GetPulsesStep(ljh, [:pulse, :rowcount, :timestamp_posix_usec], 0,100, 0)
    @perpulse pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv, peak_index, peak_value, min_value = compute_summary(pulse, pretrig_nsamples, frametime)
    @threshold pretrig_rms_criteria, postpeak_deriv_criteria = estimate_pretrig_rms_and_postpeak_deriv_criteria(noise_filename, pretrig_nsamples) when length(pulse) > 100
    @threshold peak_index_criteria = estimate_peak_index_criteria(peak_index) when length(peak_index)>100
    @threshold calibration_pulse_rms = calibrate_xy(calibration_pulse_rms_pulseheights, calibration_pulse_rms_energies) when length(peak_index)>0
    @threshold calibration_filt_value = calibrate_xy(calibration_filt_value_pulseheights, calibration_filt_value_energies) when length(peak_index)>0
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
    :ljh_filename=>"ljh_filename", :noise_filename=>"noise_filename"],
    mc[:hdf5_filename])
    MemoryLimitStep(Int(4e6)) # throw error if mc uses more than 4 MB
    FreeMemoryStep()
    ] # end steps
    setsteps!(mc, steps)
    preptasks!(mc)
    mc
end
begin
rmf = ReferenceMicrocalFiles.dict["good_mnka_mystery"]
mc = setup_channel(rmf.filename, rmf.noise_filename)
h5open(mc[:init_hdf5_filename],"w") do h5
  h5["calibrations/calibration_pulse_rms_energies"] = collect(1:4.0)
  h5["calibrations/calibration_pulse_rms_pulseheights"] = collect(1:4.0)*1000
  h5["calibrations/calibration_filt_value_energies"] = collect(1:4.0)
  h5["calibrations/calibration_filt_value_pulseheights"] = collect(1:4.0)*1000
  filter = zeros(mc[:samples_per_record])
  filter[div(end,2)]=1
  h5["filter/filter"] = filter
  end
schedule(mc)
sleep(1)
plantoend(mc)
end

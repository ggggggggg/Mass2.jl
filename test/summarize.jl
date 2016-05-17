using Mass2, Base.Test
let
  # make a single noiseless pulse with known properties
  num_points = 1024
  rise_points = 50
  fall_points = 100
  quiescent_value = 1000
  arrival_point = div(num_points,4)
  amplitude = 10000
  pulse_float = two_exponential_pulses(num_points, rise_points, fall_points, quiescent_value, (arrival_point,), (amplitude,))
  pulse = to_type_and_white_noise(UInt16, 0, pulse_float)
  pulse_int = convert(Vector{Int}, pulse) # for functions where uint causes problems
  frame_time = 5e-6 # time per point, seconds

  # summarize it
  out = compute_summary((pulse,),arrival_point,frame_time)

  (pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv,
  peak_index, peak_value, min_value) = map(first,out)
  npre = arrival_point # number of pretrigger points
  npost = length(pulse)-npre # number of posttrigger points

  @test pretrig_mean ≈ mean(pulse[1:npre]) # average of points in pretrigger
  @test pretrig_rms ≈ std(pulse[1:npre]) # std deviatino of points in pretrigger
  @test pulse_average ≈ mean(pulse[npre+1:end]-quiescent_value)# pulse average deviation from pretrig_mean
  @test pulse_rms ≈ sqrt(mean((pulse_int[npre+1:end]-quiescent_value).^2)) # pulse rms relative to pretrig_rms
  @test rise_time ≈ frame_time*(find(pulse.>10000)[1]-find(pulse.>2000)[1]) # defined as 10%-90% rise time
  # @test postpeak_deriv<minimum(diff(pulse_int))/frame_time # for a noiseless pulse, the postpeak deriv should be negative
  @test peak_index ≈ indmax(pulse) # index of the maximum value
  @test peak_value+pretrig_mean ≈ maximum(pulse) #peak value is defined relative to pretrig_mean
  @test min_value ≈ minimum(pulse) #min_value is the minimum value, not defined relative to anything
end

let
  # make a piled up noiseless pulse with known properties, peak is on 2nd pulse
  num_points = 1024
  rise_points = 50
  fall_points = 100
  quiescent_value = 1000
  arrival_point = div(num_points,4)
  arrival_point2 = div(num_points,2)
  amplitude = 10000
  amplitude2 = 10000
  pulse_float = two_exponential_pulses(num_points, rise_points, fall_points, quiescent_value,
                (arrival_point,arrival_point2), (amplitude, amplitude2))
  pulse = to_type_and_white_noise(UInt16, 0, pulse_float)
  pulse_int = convert(Vector{Int}, pulse) # for functions where uint causes problems
  frame_time = 5e-6 # time per point, seconds

  # summarize it
  out = compute_summary((pulse,),arrival_point,frame_time)

  (pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv,
  peak_index, peak_value, min_value) = map(first,out)
  npre = arrival_point # number of pretrigger points
  npost = length(pulse)-npre # number of posttrigger points

  @test pretrig_mean ≈ mean(pulse[1:npre]) # average of points in pretrigger
  @test pretrig_rms ≈ std(pulse[1:npre]) # std deviatino of points in pretrigger
  @test pulse_average ≈ mean(pulse[npre+1:end]-quiescent_value)# pulse average deviation from pretrig_mean
  @test pulse_rms ≈ sqrt(mean((pulse_int[npre+1:end]-quiescent_value).^2)) # pulse rms relative to pretrig_rms
  @test rise_time > frame_time*(find(pulse.>10000)[1]-find(pulse.>2000)[1]) # defined as 10%-90% rise time, will be confused by pileup
  # @test postpeak_deriv<minimum(diff(pulse_int))/frame_time # for a noiseless pulse, the postpeak deriv should be negative
  @test peak_index ≈ indmax(pulse) # index of the maximum value
  @test peak_value+pretrig_mean ≈ maximum(pulse) #peak value is defined relative to pretrig_mean
  @test min_value ≈ minimum(pulse) #min_value is the minimum value, not defined relative to anything
end

let
  # make a piled up noiseless pulse with known properties, peak is on 1st pulse
  num_points = 1024
  rise_points = 50
  fall_points = 100
  quiescent_value = 1000
  arrival_point = div(num_points,4)
  arrival_point2 = div(num_points,2)
  amplitude = 10000
  amplitude2 = 5000
  pulse_float = two_exponential_pulses(num_points, rise_points, fall_points, quiescent_value,
                (arrival_point,arrival_point2), (amplitude, amplitude2))
  pulse = to_type_and_white_noise(UInt16, 0, pulse_float)
  pulse_int = convert(Vector{Int}, pulse) # for functions where uint causes problems
  frame_time = 5e-6 # time per point, seconds

  # summarize it
  out = compute_summary((pulse,),arrival_point,frame_time)

  (pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv,
  peak_index, peak_value, min_value) = map(first,out)
  npre = arrival_point # number of pretrigger points
  npost = length(pulse)-npre # number of posttrigger points

  @test pretrig_mean ≈ mean(pulse[1:npre]) # average of points in pretrigger
  @test pretrig_rms ≈ std(pulse[1:npre]) # std deviatino of points in pretrigger
  @test pulse_average ≈ mean(pulse[npre+1:end]-quiescent_value)# pulse average deviation from pretrig_mean
  @test pulse_rms ≈ sqrt(mean((pulse_int[npre+1:end]-quiescent_value).^2)) # pulse rms relative to pretrig_rms
  @test rise_time ≈ frame_time*(find(pulse.>10000)[1]-find(pulse.>2000)[1]) # defined as 10%-90% rise time
  @test postpeak_deriv ≈ maximum(diff(pulse_int)[peak_index:end]) # for a noiseless pulse, the postpeak deriv should be negative
  @test peak_index ≈ indmax(pulse) # index of the maximum value
  @test peak_value+pretrig_mean ≈ maximum(pulse) #peak value is defined relative to pretrig_mean
  @test min_value ≈ minimum(pulse) #min_value is the minimum value, not defined relative to anything
end

let
  # make a single noisy pulse with known properties
  num_points = 10240
  rise_points = 50
  fall_points = 100
  quiescent_value = 1000
  arrival_point = div(num_points,4)
  amplitude = 10000
  noise_std = 1
  pulse_float = two_exponential_pulses(num_points, rise_points, fall_points, quiescent_value, (arrival_point,), (amplitude,))
  pulse = to_type_and_white_noise(UInt16, noise_std, pulse_float)
  pulse_int = convert(Vector{Int}, pulse) # for functions where uint causes problems
  frame_time = 5e-6 # time per point, seconds

  # summarize it
  out = compute_summary((pulse,),arrival_point,frame_time)

  (pretrig_mean, pretrig_rms, pulse_average, pulse_rms, rise_time, postpeak_deriv,
  peak_index, peak_value, min_value) = map(first,out)
  npre = arrival_point # number of pretrigger points
  npost = length(pulse)-npre # number of posttrigger points

  @test isapprox(pretrig_mean, mean(pulse[1:npre]) ;rtol=1e-3) # average of points in pretrigger
  @test isapprox(pretrig_rms, std(pulse[1:npre]) ;rtol=1e-3) # std deviatino of points in pretrigger
  @test isapprox(pulse_average, mean(pulse_int[npre+1:end]-quiescent_value) ;rtol=1e-3)# pulse average deviation from pretrig_mean
  @test isapprox(pulse_rms, sqrt(mean((pulse_int[npre+1:end]-quiescent_value).^2));rtol=1e-3) # pulse rms relative to pretrig_rms
  @test rise_time ≈ frame_time*(find(pulse.>10000)[1]-find(pulse.>2000)[1]) # defined as 10%-90% rise time
  @test postpeak_deriv ≈ maximum(diff(pulse_int)[peak_index:end]) # for a noiseless pulse, the postpeak deriv should be negative
  @test peak_index ≈ indmax(pulse) # index of the maximum value
  @test isapprox(peak_value+pretrig_mean, maximum(pulse) ;rtol=1e-3) #peak value is defined relative to pretrig_mean
  @test min_value ≈ minimum(pulse) #min_value is the minimum value, not defined relative to anything
end

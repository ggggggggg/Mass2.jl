abstract PulseGenerator

type TwoExponentialPulseGenerator{T<:Real} <: PulseGenerator
	record_length::Int
    pre_rise_points::Int
	rise_points::Float64
	fall_points::Float64
	noise_fwhm::Float64
	event_rate_hz::Float64
    samples_per_second::Float64
    min_points_since_last_pulse::Float64
    amplitude_distribution::Distribution
    numrows::Int
    rownumber::Int
    quiescent_average::Float64
    quiescent_modulation_period_points::Float64
    quiescent_modulation_amplitude::Float64
    last_rowstamp::Int
    trigger_points_std::Float64 # std of where it triggers relative to actualy arrival time in points (aka frames)
    d_amp_d_quiescent::Float64 # pretrigger mean/amplitude correlation
    function TwoExponentialPulseGenerator(record_length,pre_rise_points,rise_points,fall_points,noise_fwhm,event_rate_hz,samples_per_second,min_points_since_last_pulse,amplitude_distribution,numrows,rownumber,quiescent_average,quiescent_modulation_period_points,quiescent_modulation_amplitude,last_rowstamp,trigger_points_std,d_amp_d_quiescent)
    	rise_points<fall_points || error("rise_points should be less than fall points")
    	new(record_length,pre_rise_points,rise_points,fall_points,noise_fwhm,event_rate_hz,samples_per_second,min_points_since_last_pulse,amplitude_distribution,numrows,rownumber,quiescent_average,quiescent_modulation_period_points,quiescent_modulation_amplitude,last_rowstamp,trigger_points_std,d_amp_d_quiescent)
    end
end
TupacLikeTwoExponentialPulseGenerator{T}(::Type{T},distribution) = TwoExponentialPulseGenerator{T}(520, 100, 50, 200, 5, 13.5,
	104166.6, 520, distribution, 30, 0, 1000, 100_000*3600, 50, 0, 0.5, -0.0001)

function two_exponential_pulses(
	points::Int, rise_points, fall_points,quiescent_value, arrival_point::Tuple, amplitude::Tuple)
	pulse = Array(Float64, points)
	fill!(pulse, quiescent_value)
	@assert length(arrival_point) == length(amplitude)
	k=kfactor(1/rise_points,1/fall_points) 
	for (ap, amp) in zip(arrival_point, amplitude)
		@inbounds for i = 1:points
			t = i-ap
			pulse[i] += t>0.0?amp*k*(exp(-t/rise_points)-exp(-t/fall_points)):0.0
		end
	end
	pulse
end

# factor for calculating mathematical parameters from physical parameters and maintainging
# constant amplitude
# Otto-von-Guericke Univ. Magdeburg, Magdeburg, Germany
# IEEE Transactions on Electromagnetic Compatibility (Impact Factor: 1.35). 12/2010; DOI: 10.1109/TEMC.2010.2052621 
# Estimation of the Mathematical Parameters of Double-Exponential Pulses Using the Nelder–Mead Algorithm - ResearchGate. Available from: http://www.researchgate.net/publication/224153186_Estimation_of_the_Mathematical_Parameters_of_Double-Exponential_Pulses_Using_the_NelderMead_Algorithm [accessed Mar 20, 2015].kfactor(a,b) = 1/(exp(-a*(log(a/b)/(a-b)))-exp(-b*(log(a/b)/(a-b))))
kfactor(a,b) = 1/(exp(-a*(log(a/b)/(a-b)))-exp(-b*(log(a/b)/(a-b)))) # memoize?
const X = [0.842411, 1.660852, 8.887335, 0.624129]
const Y = [2.167370, 0.360349, 0.013289, 1.358701, 0.137057, 122.802861, 1.314850] 
alpha(t_rise, t_fwhm) = (X[1]/((t_fwhm/t_rise)^X[2]-X[3])^X[4])/t_rise
beta(t_rise, t_fwhm) = (Y[1]-Y[2]*exp(-t_fwhm*Y[3]/t_rise)-Y[4]*exp(-t_fwhm*Y[5]/t_rise)-Y[6]*exp(-t_fwhm*Y[7]/t_rise))/t_rise

white_noise(points::Int, noise_fwhm) = rand(points)*noise_fwhm/(2*sqrt(2*log(2))) # fwhm to std
function white_noise!(pulse, noise_fwhm) 
	@inbounds for i = 1:length(pulse) pulse[i]+=randn()*noise_fwhm end
	pulse
end
to_type_and_white_noise{T<:Integer}(::Type{T}, noise_fwhm, pulse) = convert(Vector{T}, round(white_noise!(pulse, noise_fwhm)))


function getcleanpulse{T}(pg::TwoExponentialPulseGenerator{T}, amp)
	to_type_and_white_noise(T, pg.noise_fwhm, 
		two_exponential_pulses(pg.record_length, pg.rise_points, pg.fall_points, pg.quiescent_average, (pg.pre_rise_points+1,), (amp,)))
end
function gettriggeredpulse!{T}(pg::TwoExponentialPulseGenerator{T})
	points_from_last, points_to_next = rand(Exponential(pg.samples_per_second/pg.event_rate_hz),2)
	arrival_points = (pg.pre_rise_points-pg.min_points_since_last_pulse-points_from_last+1
		, pg.pre_rise_points+rand(Normal(0,pg.trigger_points_std))+1, pg.pre_rise_points+points_to_next+1)

	last_rowstamp = pg.last_rowstamp
	rowstamp = round(last_rowstamp+points_from_last*pg.numrows)
	pg.last_rowstamp=rowstamp

	quiescent_value_modulation = pg.quiescent_modulation_amplitude*sinpi(2*rowstamp/(pg.quiescent_modulation_period_points*pg.numrows))
	quiescent_value = pg.quiescent_average+quiescent_value_modulation

	amplitudes = rand(pg.amplitude_distribution,3)

	# introduce pretrigger mean/amplitude correlation
	amplitudes ./= 1+pg.d_amp_d_quiescent*quiescent_value_modulation
	amplitudes_tuple = tuple(amplitudes...)

	pulse = two_exponential_pulses(pg.record_length, pg.rise_points, pg.fall_points, quiescent_value, arrival_points, amplitudes_tuple)
	pulse = to_type_and_white_noise(T, pg.noise_fwhm, pulse)
	last_rowstamp = pg.last_rowstamp
	rowstamp = round(last_rowstamp+points_from_last*pg.numrows)
	pg.last_rowstamp=rowstamp
	pulse, rowstamp
end
function gettriggeredpulse!{T}(pg::TwoExponentialPulseGenerator{T},n)
	pulses = Array(Vector{T},n)
	rowstamps = Array(Int, n)
	for i = 1:n
		pulse, rowstamp = gettriggeredpulse!(pg)
		pulses[i] = pulse
		rowstamps[i] = rowstamp
	end
	pulses, rowstamps
end
times_s(pg) = range(0,1/pg.samples_per_second, pg.record_length)

# Base.start(pg::TwoExponentialPulseGenerator) = 0.0
# function Base.next(pg::TwoExponentialPulseGenerator, last_arrival::Float64)
# 	timestamp = last_arrival+rand(Exponential(pg.samples_per_second/pg.event_rate_hz),2)
# 	pulse = gettriggeredpulse(pg)
# 	return ((timestamp,pulse), timestamp)
# end
# Base.done(pg::TwoExponentialPulseGenerator, last_arrival::Float64) = false

# function getrowstamp(pg::TwoExponentialPulseGenerator,n::Int, last_rowstamp::Int=0)
# 	rowstamp_diffs = pg.numrows*rand(Exponential(pg.samples_per_second/pg.event_rate_hz),n)
# 	last_rowstamp = pg.numrows*div(last_rowstamp,pg.numrows)+pg.rownumber # make sure rowstamp reflect the correct row
# 	rowstamp_diffs[1]+=last_rowstamp
# 	cumsum(convert(Vector{Int}, round(rowstamp_diffs)))
# end


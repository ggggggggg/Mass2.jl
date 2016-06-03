using Mass2

srand(hash("mockpluses.jl"))
pg = TupacLikeTwoExponentialPulseGenerator(Int, Normal(1000,0.5))
getcleanpulse(pg,1000)
gettriggeredpulse!(pg,10)

function compute_whitenoise_filter(pulse, selection_good)
	filter = mean(pulse[selection_good])
	filter -= mean(filter) # 0 average
	normalization = (maximum(filter)-minimum(filter))./dot(filter, filter)
	filter*normalization
end
function filter1lag(pulse, filter)
	out = Array(Float64, length(pulse))
	for i=1:length(pulse)
		out[i] = dot(pulse[i], filter)
	end
	out
end
cleanpulses = [getcleanpulse(pg,1000) for i = 1:10000]
filter = compute_whitenoise_filter(cleanpulses, trues(length(cleanpulses)))
cleanpulses_filt_values = filter1lag(cleanpulses, filter)
@test abs(mean(cleanpulses_filt_values)-1000)<1
@test abs(std(cleanpulses_filt_values)-.655)<.1

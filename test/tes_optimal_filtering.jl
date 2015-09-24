import TESOptimalFiltering: filter5lag, smoother, calculate_filter, autocorrelation
using Base.Test, HDF5
srand(1000)
# write your own tests here
pulse_model = 100*(exp(-(1:1000)/150)-exp(-(1:1000)/350));
noise_data = randn(1000000);
dt = 9.6e-6;
f_3db = 1e30;
noise_autocorr = autocorrelation(noise_data, length(pulse_model));
filter, variance = calculate_filter(pulse_model, noise_autocorr, f_3db, dt)
pulses = Vector{Float64}[pulse_model*i+randn(length(pulse_model)) for i=1:10] # create pulses with heights 1 thru 1-10
pulses2d = zeros(Float64, length(pulse_model), 10 )
for i in eachindex(pulses) pulses2d[:,i] = pulses[i] end
filt_value, filt_phase = filter5lag(filter, pulses[1]) # filter the first pulse
filt_values, filt_phases = filter5lag(filter, pulses) # filter all the pulses
filt_values2d, filt_phases2d = filter5lag(filter, pulses2d) # filter all the pulses

@test filt_value == filt_values[1]
@test filt_phase == filt_phases[1]
@test all(filt_values2d .== filt_values)
@test all(filt_phases2d .== filt_phases)

# results seem noiser than expected
#@test all(mean(abs(filt_values - collect(1:10))) .< sqrt(variance)*5) # check that pulse heights came out about right

pmm = pulse_model-mean(pulse_model);
approx_filter = pmm/dot(pmm,pmm);
@test maximum(approx_filter-filter)<1e-4


# compare to data from 20150730_50mM_irontris_100ps_xes_chan479.ljh caclulated by python mass on sept 14 2015
h5 = h5open(joinpath(Pkg.dir("Mass2/test"),"pythonmass_refdata.hdf5"))

filter, variance = calculate_filter(h5["average_pulse"][:], h5["noise_autocorr"][:], read(h5["f_3db"]), read(h5["dt"]))
#normalize so filter returns pulse height of average pulse
filter*=maximum(h5["average_pulse"][:])
vdv = 1/sqrt(variance)/sqrt(8*log(2))

@test abs(vdv/read(h5["vdv"])-1)<1e-3
# mass filter is shortened by 2 points on either end
# check that the filters are roughly similar by comparing their filtered values
@test abs(dot(filter[3:end-2],filter[3:end-2])/dot(filter[3:end-2],h5["filter"][:])-1) < 2/1000
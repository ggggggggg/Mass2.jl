# test peakassign
known_energies = [3691, 4090, 4510, 4952, 5414, 5898, 6403, 6929, 7477, 8040, 8637] # kalpha energies Ca-Zn

filt_values = known_energies.^0.9 # caluclate filt values with some amount of non-linearity

tested = 0
peakassignwrong=0 

nc = 7
tested+=length(combinations(1:length(known_energies),nc))
for combo in combinations(1:length(known_energies),nc)
	locations_arb = filt_values
	locations_true = known_energies[combo]
	best_combo, best_fom = peakassign(locations_arb, locations_true)
	if !all(best_combo .== filt_values[combo])
		peakassignwrong+=1
	end
end

@test peakassignwrong == 0

# test findpeaks
energy_value_distribution = MixtureModel([Normal(mu,15) for mu in known_energies])

be = 0:1:20000 # bin edges
events = rand(energy_value_distribution, 10000)
be, c = hist(events, be)

# plotpeaks(be, c, 10,10) # dont plot in regular testing
peakinds = findpeaks(midpoints(be), c, fwhm=50,n=length(known_energies))
@test length(peakinds)==length(known_energies)
peakx = sort(midpoints(be)[peakinds])
@test all(abs(peakx-known_energies).<5);


# test findpeaks on Histogram
bin_edges_v = 0:1:20000
h = Histogram(bin_edges_v)
addcounts!(h, events)
peakinds = findpeaks(h,fwhm=50,n=length(known_energies))
@test length(peakinds)==length(known_energies)
peakx = sort(midpoints(be)[peakinds])
@test all(abs(peakx-known_energies).<5);
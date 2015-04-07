known_energies = [3691, 4090, 4510, 4952, 5414, 5898, 6403, 6929, 7477, 8040, 8637] # kalpha energies Ca-Zn

filt_values = known_energies.^0.9 # caluclate filt values with some amount of non-linearity

tested = 0
peakassignwrong=0 

nc = 7
tested+=length(combinations(1:length(known_energies),nc))
for c in combinations(1:length(known_energies),nc)
	locations_arb = filt_values
	locations_true = known_energies[c]
	best_combo, best_fom = peakassign(locations_arb, locations_true)
	if !all(best_combo .== filt_values[c])
		peakassignwrong+=1
	end
end

println("peakassign test got $(peakassignwrong) fraction wrong of $tested")

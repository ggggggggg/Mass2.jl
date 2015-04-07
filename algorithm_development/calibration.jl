function polyfit(x, y, n)
  A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
  A \ y
end

function polyval(p,x::Number)
	y = zero(promote_type(eltype(p), eltype(x)))
	for i=1:length(p)
		y+=p[i]*x^(i-1)
	end
	y
end
function polyval(p,x::AbstractArray)
	y = zeros(promote_type(eltype(p), eltype(x)), length(x))
	for i = 1:length(x)
		y[i] = polyval(p,x[i])
	end
	y
end
function polyresiduals(p,x,y)
	yp = polyval(p,x)
	residuals = y-yp
end


### CALIBRATION
function peakassignlog(locations_arb, locations_true)
# locations_arb is an array of possible peak locations in arbitrary units, assumed to have an unknown, 
# nonlinear (but monotonic and smooth) relationship to locations_true
# locations_arb may contain more peaks than locations_true, including some spurious peaks
# the goal is to find the best possible set of assignments for entries in locations_arb to 
# correspond to locations_true by examining all possible assignments and calcualting some figure of merit
@assert(issorted(locations_arb))
@assert(issorted(locations_true))
@assert(length(locations_arb)>=length(locations_true))
combos = combinations(1:length(locations_arb), length(locations_true))
length(combos)==1 && return (locations_arb[first(combos)], 0.0)

assert(!(0 in locations_true)) # don't pass zero and force zero
assert(!(0 in locations_arb))

log_locations_arb = log(locations_arb)
log_locations_true = log(locations_true)

best_combo = first(combos)
best_fom = realmax(Float64)
for c in combos
	cu = log_locations_arb[c]
	fom = peakassign_fig_of_merit2(cu, log_locations_true)
	if fom < best_fom
		best_fom = fom
		best_combo = c
	end
end
locations_arb[best_combo], best_fom
end

function peakassign2(locations_arb, locations_true)
# locations_arb is an array of possible peak locations in arbitrary units, assumed to have an unknown, 
# nonlinear (but monotonic and smooth) relationship to locations_true
# locations_arb may contain more peaks than locations_true, including some spurious peaks
# the goal is to find the best possible set of assignments for entries in locations_arb to 
# correspond to locations_true by examining all possible assignments and calcualting some figure of merit
@assert(issorted(locations_arb))
@assert(issorted(locations_true))
@assert(length(locations_arb)>=length(locations_true))
combos = combinations(1:length(locations_arb), length(locations_true))
length(combos)==1 && return (locations_arb[first(combos)], 0.0)

assert(!(0 in locations_true)) # don't pass zero and force zero
assert(!(0 in locations_arb))



best_combo = first(combos)
best_fom = realmax(Float64)
for c in combos
	cu = locations_arb[c]
	fom = peakassign_fig_of_merit2(cu, locations_true)
	if fom < best_fom
		best_fom = fom
		best_combo = c
	end
end
locations_arb[best_combo], best_fom
end

function peakassign2zeros(locations_arb, locations_true)
# locations_arb is an array of possible peak locations in arbitrary units, assumed to have an unknown, 
# nonlinear (but monotonic and smooth) relationship to locations_true
# locations_arb may contain more peaks than locations_true, including some spurious peaks
# the goal is to find the best possible set of assignments for entries in locations_arb to 
# correspond to locations_true by examining all possible assignments and calcualting some figure of merit
@assert(issorted(locations_arb))
@assert(issorted(locations_true))
@assert(length(locations_arb)>=length(locations_true))
combos = combinations(1:length(locations_arb), length(locations_true))
length(combos)==1 && return (locations_arb[first(combos)], 0.0)

assert(!(0 in locations_true)) # don't pass zero and force zero
assert(!(0 in locations_arb))
locations_true = unshift!(copy(locations_true), 0)


best_combo = first(combos)
best_fom = realmax(Float64)
for c in combos
	cu = locations_arb[c]
	unshift!(cu,0)
	fom = peakassign_fig_of_merit2(cu, locations_true)
	if fom < best_fom
		best_fom = fom
		best_combo = c
	end
end
locations_arb[best_combo], best_fom
end

function peakassign2(locations_arb, locations_true)
# locations_arb is an array of possible peak locations in arbitrary units, assumed to have an unknown, 
# nonlinear (but monotonic and smooth) relationship to locations_true
# locations_arb may contain more peaks than locations_true, including some spurious peaks
# the goal is to find the best possible set of assignments for entries in locations_arb to 
# correspond to locations_true by examining all possible assignments and calcualting some figure of merit
@assert(issorted(locations_arb))
@assert(issorted(locations_true))
@assert(length(locations_arb)>=length(locations_true))
combos = combinations(1:length(locations_arb), length(locations_true))
length(combos)==1 && return (locations_arb[first(combos)], 0.0)

assert(!(0 in locations_true)) # don't pass zero and force zero
assert(!(0 in locations_arb))



best_combo = first(combos)
best_fom = realmax(Float64)
for c in combos
	cu = locations_arb[c]
	fom = peakassign_fig_of_merit(cu, locations_true)
	if fom < best_fom
		best_fom = fom
		best_combo = c
	end
end
locations_arb[best_combo], best_fom
end

function peakassign_fig_of_merit(x, y)
	# average fractional error trying to predict the y value of the middle of 3 points using the slope
	# from the outer two points and the x value of the middle point
	err = 0.0
	for i = 1:length(y)-2
		yip1_prediction = y[i]+(y[i+2]-y[i])*(x[i+1]-x[i])/(x[i+2]-x[i])
		err += abs(1-y[i+1]/yip1_prediction)
	end
	err/length(y-2)
end

function peakassign_fig_of_merit2(x, y)
	n = length(x) > 2 ? 2 : 1
	p = polyfit(x,y,n)
	std(polyresiduals(p,x,y)./y)
end

# test peakassign

known_energies = [3691, 4090, 4510, 4952, 5414, 5898, 6403, 6929, 7477, 8040, 8637] # ka energies Ca-Zn

filt_values = known_energies.^0.9 # caluclate filt values with some amount of non-linearity

tested = 0
peakassignlogwrong=0 
peakassign2wrong=0 
peakassign2zerowrong=0 
peakassignwrong=0 

nc = 7
	tested+=length(combinations(1:length(known_energies),nc))
	println("testing peakassignlog $nc")
	for c in combinations(1:length(known_energies),nc)
		locations_arb = filt_values
		locations_true = known_energies[c]
		best_combo, best_fom = peakassignlog(locations_arb, locations_true)
		if !all(best_combo .== filt_values[c])
			peakassignlogwrong+=1
			# println(locations_true)
			# println(best_combo)
			# println(filt_values[c])
			# @test best_fom < 1
			# @test all(best_combo .== filt_values[c])
		end
	end

	println("testing peakassign2 $nc")
	for c in combinations(1:length(known_energies),nc)
		locations_arb = filt_values
		locations_true = known_energies[c]
		best_combo, best_fom = peakassign2(locations_arb, locations_true)
		if !all(best_combo .== filt_values[c])
			peakassign2wrong+=1

		end
	end

	println("testing peakassign2zeros $nc")
	for c in combinations(1:length(known_energies),nc)
		locations_arb = filt_values
		locations_true = known_energies[c]
		best_combo, best_fom = peakassign2zeros(locations_arb, locations_true)
		if !all(best_combo .== filt_values[c])
			peakassign2zerowrong+=1

		end
	end
	println("testing peakassign $nc")
	for c in combinations(1:length(known_energies),nc)
		locations_arb = filt_values
		locations_true = known_energies[c]
		best_combo, best_fom = peakassign2(locations_arb, locations_true)
		if !all(best_combo .== filt_values[c])
			peakassignwrong+=1

		end
	end

println("peakassignlog $(peakassignlogwrong/tested)")
println("peakassign2 $(peakassign2wrong/tested)")
println("peakassign2zeros $(peakassign2zerowrong/tested)")
println("peakassign $(peakassignwrong/tested)")

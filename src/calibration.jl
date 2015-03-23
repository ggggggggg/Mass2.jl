using Dierckx, PyCall

abstract AbstractCalibration

type CalibrationSpline1D <: AbstractCalibration
	spline::Spline1D
	known_energies::Vector{Float64}
	indicator::Vector{Float64}
	iscalibrated::Bool
end
function CalibrationSpline1D(indicator, known_energies)
	CalibrationSpline1D(calc_spline(indicator, known_energies), known_energies, indicator, true)
end
iscalibrated(c::AbstractCalibration) = c.iscalibrated
donethru(c::AbstractCalibration) = iscalibrated(c)?DONETHRU_MAX:0
apply_calibration(c::CalibrationSpline1D, estimator) = apply_spline(c.spline, estimator)


### CALIBRATION
function peakassign(locations_arb, locations_true)
# locations_arb is an array of possible peak locations in arbitrary units, assumed to have an unknown, 
# nonlinear (but monotonic and smooth) relationship to locations_true
# locations_arb may contain more peaks than locations_true, including some spurious peaks
# the goal is to find the best possible set of assignments for entries in locations_arb to 
# correspond to locations_true by examining all possible assignments and calcualting some figure of merit
@assert(issorted(locations_arb))
@assert(issorted(locations_true))
@assert(length(locations_arb)>=length(locations_true))
combos = combinations(locations_arb, length(locations_true))
length(combos)==1 && return (first(combos), inf)

best_combo = first(combos)
best_fom = realmax(Float64)
for c in combos
	fom = peakassign_fig_of_merit(c, locations_true)
	if fom < best_fom
		best_fom = fom
		best_combo = c
	end
end
best_combo, best_fom
end

function peakassign_fig_of_merit(x, y)
	# average fractional error trying to predict the y value of the middle of 3 points using the slope
	# from the outer two points and the x value of the middle point
	err = 0.0
	for i = 1:length(y)-2
		yip1_prediction = y[i]+(y[i+2]-y[i])*(x[i+1]-x[i])/(x[i+2]-x[i])
		err += abs(1-y[i+1]/yip1_prediction)
	end
	err /= length(y-2)
end

@pyimport scipy.signal as scipysignal
function findpeaks(y, length_scales; min_snr = 3, min_dist=-1)
	# find peaks using scipy.signal.find_peaks_cwt
	# then make sure any reported peaks are local maxima
	# returnds indicies of peaks orderes from lowest to highest y value
	raw_peakinds = int(scipysignal.find_peaks_cwt(y, length_scales, min_snr=min_snr))
	peakinds = Int[]
	min_dist == -1 && (min_dist = int(maximum(length_scales)/2))
	for p in raw_peakinds
		r = max(1,p-min_dist):min(length(y), p+min_dist)
		local_max_ind = indmax(y[r])+first(r)-2
		local_max_ind in peakinds || push!(peakinds, local_max_ind)
	end
	perm = sortperm(y[peakinds])
	peakinds[perm]
end


function findpeaks(hist::Histogram) 
	peak_inds = findpeaks(counts(hist),[2,4,8,16,32].*binsize(hist))
	peak_x = bin_centers(hist)[peak_inds]
	peak_inds, peak_x
end


function calibrate_nofit(hist::Histogram, known_energies; forcezero=true)
	peak_inds, peak_x = findpeaks(hist)
	println((peak_inds, peak_x))
	if forcezero
		assert(!(0 in known_energies)) # don't pass zero and force zero
		unshift!(peak_x, 0)
		unshift!(known_energies, 0)
	end
	shortened_sorted_peak_x = sort(length(peak_x)>15?peak_x[1:15]:peak_x)
	best_combo, best_fom = peakassign(sort(shortened_sorted_peak_x), known_energies)

	CalibrationSpline1D(best_combo, known_energies)
end

function tempcal(arb, known, forcezero=true)
	if forcezero
		assert(!(0 in known)) # don't pass zero and force zero
		unshift!(arb, 0)
		unshift!(known, 0)
	end
	best_combo, best_fom = peakassign(arb, known)
end

function calc_spline(x,y)
	spl = Spline1D(float(x),float(y),k=1,bc="extrapolate")
end
function apply_spline(spl, energy_estimator)
	energy = evaluate(spl, energy_estimator)
end



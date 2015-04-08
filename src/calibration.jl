
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
function PyPlot.plot(cal::CalibrationSpline1D) 
	figure()
	plot(cal.indicator, apply_spline(cal.spline, cal.indicator))
	plot(cal.indicator, cal.known_energies,"o")
	xlabel("indicator")
	ylabel("energy")
end

function plotlog(cal::CalibrationSpline1D) 
	figure()
	x=log(cal.indicator)
	y=log(cal.known_energies)
	# plot(cal.indicator, apply_spline(cal.spline, cal.indicator))
	plot(x, y,"o")
	p = polyfit(x,y,2)
	plot(x, polyval(p, x))
	xlabel("indicator")
	ylabel("energy")
	title("peakassign_fig_of_merit $(peakassign_fig_of_merit(x,y))")
end


function peakassign(locations_arb, locations_true)
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
	fom = peakassign_fig_of_merit(cu, log_locations_true)
	if fom < best_fom
		best_fom = fom
		best_combo = c
	end
end
locations_arb[best_combo], best_fom
end
function peakassign_fig_of_merit(x, y)
	n = length(x) > 2 ? 2 : 1
	p = polyfit(x,y,n)
	std(polyresiduals(p,x,y)./y)
end
# should remove these and depend on Polynomails.jl or CurveFit or something
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



findpeaks(h::Histogram; fwhm=15,n=typemax(Int)) = findpeaks(bin_edges(h), counts(h), fwhm=fwhm,n=n)
function findpeaks(bin_edges::Range, y::Vector; fwhm=15,n=typemax(Int))
	k = UnivariateKDE(bin_edges, convert(Array{Float64,1},y))
	findpeaks_(k,fwhm,n)
end

function findpeaks_(y::UnivariateKDE, fwhm,n)
	# smooth y with a gaussian of FHWM=fwhm, return indicies of n peaks with highest density in smoothed version
	σ = fwhm/sqrt(8*log(2))
	d = Normal(0,σ)
	smoothed = conv(y,d)
	peakinds = findpeaks_(smoothed.density)
	if length(peakinds)>n
		peakinds = peakinds[end-n+1:end]
	end
	peakinds
end

function findpeaks_{T}(y::Vector{T})
	# find all local maxima, return a vector of indicies sorted by value of y at that indey
	peakinds = Array(Int,0)
	for i = 2:length(y)-1
		if y[i-1]<y[i]>y[i+1]
			push!(peakinds, i)
		end
	end
	p = sortperm(y[peakinds])
	peakinds[p]
end

plotpeaks(h::Histogram; fwhm=15, n=typemax(Int)) = plotpeaks(bin_edges(h), counts(h), fwhm=fwhm, n=n)
function plotpeaks(bin_edges::Range, y::Vector; fwhm=15, n=typemax(Int))
	k = UnivariateKDE(bin_edges, convert(Array{Float64,1},y))
	peakinds = findpeaks_(k, fwhm, n)
	peaky = midpoints(bin_edges)[peakinds]
	figure()
	plot(midpoints(bin_edges), y, label="y")
	plot(midpoints(k.x), conv(k,Normal(0,fwhm/2.355)).density, label="smoothed")
	plot(peaky, y[peakinds],"o",label="raw peaks")
	plot(peaky, conv(k,Normal(0,fwhm/2.355)).density[peakinds],"o",label="smoothed peaks")
	legend()
end


function calibrate_nofit(hist::Histogram, known_energies)
	peak_inds = findpeaks(hist,fwhm=15, n=length(known_energies)+3) # find the largest peaks, at most 3 more than the number of known energies
	# with a 15 unit gaussian applied for smoothing... this needs to be adjustable
	# also we need to check if calibration works, and if not return a failed calibration object with debug info
	peak_x = sort(bin_centers(hist)[peak_inds])
	best_combo, best_fom = peakassign(sort(peak_x), known_energies)
	CalibrationSpline1D(best_combo, known_energies)
end


function calc_spline(x,y)
	spl = Spline1D(float(x),float(y),k=1,bc="extrapolate")
end
function apply_spline(spl, energy_estimator)
	energy = evaluate(spl, energy_estimator)
end



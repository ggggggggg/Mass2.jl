# generate a spectra with multiple peaks for calibration testing
dist_energies = [4000,4700,5500,6000,6800,7500,8200,9600,10000]
energy_value_distribution = MixtureModel([Normal(mu,15) for mu in dist_energies])

be = 0:1:20000

events = rand(energy_value_distribution, 10000)
be, c = hist(events, be)
k=UnivariateKDE(be, float(c))


findpeaks(h::Histogram, fwhm) = findpeaks(bin_edges(h), counts(h), fwhm, typemax(Int))
findpeaks(h::Histogram, fwhm,n) = findpeaks(bin_edges(h), counts(h), fhwm,n)

findpeak(bin_edges::Range, x::Vector, fwhm) = findpeaks(bin_edges, x, fwhm, typemax(Int))
function findpeaks(bin_edges::Range, x::Vector, fwhm,n)
	k = UnivariateKDE(bin_edges, convert(Array{Float64,1},x))
	findpeaks(k,fwhm,n)
end

function findpeaks(x::UnivariateKDE, fwhm,n)
	σ = fwhm/sqrt(8*log(2))
	d = Normal(0,σ)
	smoothed = conv(x,d)
	peakinds = findpeaks_(smoothed.density)
	if length(peakinds)>n
		peakinds = peakinds[end-n+1:end]
	end
	peakinds
end

function findpeaks_{T}(x::Vector{T})
	# find all local maxima, return a vector of indicies sorted by value of x at that index
	peakinds = Array(Int,0)
	for i = 2:length(x)-1
		if x[i-1]<x[i]>x[i+1]
			push!(peakinds, i)
		end
	end
	p = sortperm(x[peakinds])
	peakinds[p]
end

function plotpeaks(bin_edges::Range, x::Vector, fwhm, n)
	k = UnivariateKDE(bin_edges, convert(Array{Float64,1},x))
	peakinds = findpeaks(bin_edges, x, fwhm, n)
	peakx = midpoints(bin_edges)[peakinds]
	figure()
	plot(midpoints(bin_edges), x, label="x")
	plot(midpoints(k.x), conv(k,Normal(0,fwhm/2.355)).density, label="smoothed")
	plot(peakx, x[peakinds],"o",label="raw peaks")
	plot(peakx, conv(k,Normal(0,fwhm/2.355)).density[peakinds],"o",label="smoothed peaks")
	legend()
end


peakinds = findpeaks(k,50,length(dist_energies))
peakx = sort(midpoints(be)[peakinds])

plotpeaks(be, c, 10,10)


figure()
plot(midpoints(be),c)
plot(midpoints(k.x), conv(k,Normal(0,10/2.355)).density)
plot(midpoints(be)[peakinds], c[peakinds],"o")
plot(midpoints(k.x)[peakinds], conv(k,Normal(0,10/2.355)).density[peakinds],"o")


@test all(abs(peakx-dist_energies).<5);
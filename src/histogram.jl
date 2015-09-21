import Base: +

type Histogram
	x::FloatRange{Float64}
	counts::Vector{Int}
	seen::Int # total events seen, including those not in counts
	counted::Int # running sum of counts
end
Histogram(x::Range) = Histogram(convert(FloatRange{Float64},x),zeros(Int, length(x)-1),0,0)
bin_edges(h::Histogram) = h.x
edges(h::Histogram) = h.x
bin_centers(h::Histogram) = midpoints(h.x)
Base.midpoints(h::Histogram) = midpoints(h.x)
binsize(h::Histogram) = step(h.x) 
function addcounts!(h::Histogram, events::Array) 
	newcounts = hist(events,edges(h))[2]
	h.counts+=newcounts
	h.counted+=sum(newcounts)
	h.seen+=length(events)
end
addseen!(h::Histogram,n)=h.seen+=n
counts(h::Histogram) = h.counts
counted(h::Histogram) = h.counted
misses(h::Histogram) = h.seen-counted(h)
donethru(h::Histogram) = h.seen
function update_histogram!(h::Histogram, selection, x)
	addcounts!(h,x[selection])
	addseen!(h,length(selection)-sum(selection))
end
function (+)(a::Histogram, b::Histogram)
	@assert edges(a)==edges(b) "adding histograms only supported if the x is the same, a has $(edges(a)), b has $(edges(b))"
	c = Histogram(edges(a), a.counts + b.counts, a.seen + b.seen, a.counted + b.counted)
end


function PyPlot.plot(h::Histogram)
	figure()
	plot(midpoints(h), counts(h))
end

# interface functions to KernelDensity
uvhist(edges, y) = UnivariateKDE(edges, convert(Array{Float64,1},hist(y,edges)[2]))
density(h::UnivariateKDE) = h.density
edges(h::UnivariateKDE) = h.x
Base.midpoints(h::UnivariateKDE) = midpoints(edges(h))
Base.conv(h::Histogram, y) = conv(UnivariateKDE(edges(h), convert(Vector{Float64},counts(h))), y)

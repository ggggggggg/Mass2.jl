import Base: +

"""`Histogram` is a type representing a histogram. It bundles up both the `x` and `y` information as well a running summs of all
seen event and all counted events. Counted events are those inside one of the bins, while seen events are any event at all
even with value outside any bin. `Histogram(r)` creates a new empty histogram with bin edges given by the range `r`.
`counts(h)[i]` is the number of events observed with value between `bin_edges(h)[i]` and `bin_edges(h)[i+1]`."""
type Histogram
	x::FloatRange{Float64}
	counts::Vector{Int}
	seen::Int # total events seen, including those not in counts
	counted::Int # running sum of counts
end
Histogram(x::Range) = Histogram(convert(FloatRange{Float64},x),zeros(Int, length(x)-1),0,0)
"""Return the bin edges of a Histogram as a FloatRange{Float64}."""
bin_edges(h::Histogram) = h.x
"""Return the bin edges of a Histogram as a FloatRange{Float64}."""
edges(h::Histogram) = h.x
"""Return the bin centers of a Histogram as a FloatRange{Float64}."""
bin_centers(h::Histogram) = midpoints(h.x)
"""Return the bin centers of a Histogram as a FloatRange{Float64}."""
Base.midpoints(h::Histogram) = midpoints(h.x)
"""Return the bin size of a Histogram as a Float64."""
binsize(h::Histogram) = step(h.x)
"""addcounts!(h::Histogram, events::Array)
Add events in the Vector `events` to `h`. If the value of the event is inside one of the bins, `counts` will increase at the corresponding index.
`seen` will increase by `length(events)` and `counted` will increase by the sum of the increase to `counts`. """
function addcounts!(h::Histogram, events::Array)
	newcounts = hist(events,edges(h))[2]
	h.counts+=newcounts
	h.counted+=sum(newcounts)
	h.seen+=length(events)
end
"""addseen!(h::Histogram,n) Add seen events to `h`. Seen events are events whose value is outside of all bins of `h`."""
addseen!(h::Histogram,n)=h.seen+=n
"""counts(h::Histogram) Retrieve the `counts` from `h`."""
counts(h::Histogram) = h.counts
"""counted(h::Histogram) Return the total number of events whose value is inside a bin of `h`."""
counted(h::Histogram) = h.counted
"""misses(h::Histogram) Return the number of events whose value is outside all bins of `h`."""
misses(h::Histogram) = h.seen-counted(h)
"""seen(h::Histogram) Return the number of events seen, counted plus misses."""
seen(h::Histogram) = h.seen
donethru(h::Histogram) = h.seen
"""update_histogram!(h::Histogram, selection, x)
Add events to `h` whose values are in `x`. Only add events for which `selection` (a bit vector or Vector{Bool}) is true. Add
events for which `selection` is false to to `h.seen`."""
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

# interface functions to KernelDensity, mainly to allow convoluation of Histogram with probability distributions like Distributions.Normal(μ=0.0, σ=1.0)
uvhist(edges, y) = UnivariateKDE(edges, convert(Array{Float64,1},hist(y,edges)[2]))
density(h::UnivariateKDE) = h.density
edges(h::UnivariateKDE) = h.x
Base.midpoints(h::UnivariateKDE) = midpoints(edges(h))
Base.conv(h::Histogram, y) = conv(UnivariateKDE(edges(h), convert(Vector{Float64},counts(h))), y)

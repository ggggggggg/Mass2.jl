type Histogram
	x::FloatRange{Float64}
	counts::Vector{Int}
	seen::Int # total events seen, including those not in counts
	counted::Int # running sum of counts
end
Histogram(x::Range) = Histogram(convert(FloatRange{Float64},x),zeros(Int, length(x)-1),0,0)
bin_edges(h::Histogram) = h.x
midpoints(r::Range) = convert(FloatRange{Float64},range(first(r)*step(r)/2,step(r),length(r)-1))
bin_centers(h::Histogram) = midpoints(h.x)
binsize(h::Histogram) = step(h.x)
function addcounts!(h::Histogram, events::Array) 
	newcounts = hist(events,bin_edges(h))[2]
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
bin_edges_v = 0:1:20000
h = Histogram(bin_edges_v)


@test all([isapprox(bin_edges(h)[i], bin_edges_v[i]) for i = 1:length(bin_edges_v)])
bin_centers_v = (bin_edges_v[2:end]+bin_edges_v[1:end-1])*0.5
@test all([isapprox(bin_centers(h)[i], bin_edges_v[i]) for i = 1:length(bin_edges_v)-1])
@test isapprox(binsize(h),1.0)

addcounts!(h, ones(100)*10)
@test counted(h)==100
@test donethru(h) == 100
@test sum(counts(h))==100

addcounts!(h, ones(100)*-10)
@test counted(h)==100
@test misses(h) == 100
@test donethru(h) == 200

s = Bool[1,1,1,1,0,0,0,0]
events = [-1,2,3,4,5,6,7,8]
update_histogram!(h, s, events)
@test counted(h) == 103
@test donethru(h) == 208
@test misses(h) == 105

@test counts(h)[10]==100
@test counts(h)[2]==1

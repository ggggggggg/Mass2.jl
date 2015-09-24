using Base.Test
using LJH
# Write file 1 as an LJH 2.1 file and files 2 as LJH 2.2 with identical data.
name1, f1 = mktemp()
name2, f2 = mktemp()

dt = 9.6e-6
npre = 200
nsamp = 1000

LJH.writeljhheader(f1, dt, npre, nsamp; version="2.1.0")
LJH.writeljhheader(f2, dt, npre, nsamp; version="2.2.0")

rowcount = collect(5000:1000:10000)
timestamps = [Int64(round(r*(dt/30)*1e6)) for r in rowcount]
N = length(rowcount)

data = rand(UInt16, nsamp, N)
data[1,:] = 0xffff

LJH.writeljhdata(f1, data, rowcount)
LJH.writeljhdata(f2, data, rowcount, timestamps)

close(f1)
close(f2)

# Now check that the data are readable
ljh1 = LJH.LJHFile(name1)
ljh2 = LJH.LJHFile(name2)

# Test that the header info is correct
for ljh in (ljh1, ljh2)
    @test ljh.nsamp == nsamp
    @test ljh.npre == npre
    @test ljh.dt == dt
    @test ljh.nrec == N
end

# Test pop!() access
for i = 1:3
    d,r,t = pop!(ljh1)
    @test d==data[:,i]
    # The encoding of rowcount as "time" and subsequent decoding means that
    # we can't insist @test r==rowcount[i] for LJH 2.1 files.
    @test r >= rowcount[i] && r < rowcount[i]+30
    @test t==timestamps[i]

    d,r,t = pop!(ljh2)
    @test d==data[:,i]
    @test r==rowcount[i]
    @test t==timestamps[i]
end

# Test extractljhdata() access and ljhrewind()
rvec = zeros(Int64, 3)
tvec = zeros(Int64, 3)
dmat = zeros(UInt16, nsamp, 3)
LJH.extractljhdata(ljh1, 3, rvec, tvec, dmat)
@test dmat==data[:,4:6]
@test all(rvec .>= rowcount[4:6]) && all(rvec .< rowcount[4:6]+30)
@test tvec==timestamps[4:6]
LJH.extractljhdata(ljh2, 3, rvec, tvec, dmat)
@test dmat==data[:,4:6]
@test rvec==rowcount[4:6]
@test tvec==timestamps[4:6]
LJH.ljhrewind(ljh2)
LJH.extractljhdata(ljh2, 3, rvec, tvec, dmat)
@test dmat==data[:,1:3]
@test rvec==rowcount[1:3]
@test tvec==timestamps[1:3]


# Test data access by ljh[i] going from 1:N, then backwards,
# then skipping data forwards and backwards.
ranges = (1:N, N:-1:1, 1:3:N, N:-2:1)
for thisrange in ranges
    for i=1:N
        d,r,t = ljh1[i]
        @test d==data[:,i]
        # The encoding of rowcount as "time" and subsequent decoding means that
        # we can't insist @test r==rowcount[i] for LJH 2.1 files.
        @test r >= rowcount[i] && r < rowcount[i]+30
        @test t==timestamps[i]

        d,r,t = ljh2[i]
        @test d==data[:,i]
        @test r==rowcount[i]
        @test t==timestamps[i]
    end
end

# Grab a slice from both LJHFile and LJHGroup
@assert N>5
thisrange = 3:5
thisslice = ljh2[thisrange]
for ((d,r,t), i) in zip(thisslice, thisrange)
    @test d==data[:,i]
    @test r >= rowcount[i] && r < rowcount[i]+30
    @test t== timestamps[i]
end

# Be sure to test LJHGroup capability for stringing multiple files together.
# First, groups corresponding to one file
grp1 = LJH.LJHGroup(ljh2)
grp2 = LJH.LJHGroup(name2)
grp3 = LJH.LJHGroup([ljh2])
for g in (grp1, grp2, grp3)
    for i=1:N
        d,r,t = g[i]
        @test d==data[:,i]
        @test r==rowcount[i]
        @test t==timestamps[i]
    end
end

# Now a group corresponding to 3 files (actually, same one 3 times)
grp4 = LJH.LJHGroup([name2, name1, name2])
for j=1:3N
    d,r,t = grp4[j]
    i = mod(j-1, N)+1
    @test d==data[:,i]
    @test r >= rowcount[i] && r < rowcount[i]+30
    @test t==timestamps[i]
end

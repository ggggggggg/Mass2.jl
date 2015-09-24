"LJH is a module for working with LJH files. Unless you have good reason, you should probably use `LJHGroup` most of the time."
module LJH

export LJHGroup, LJHFile, update_num_records, channel, record_nsamples,
        pretrig_nsamples, frametime, filenames, lengths,
        column, row, num_columns, num_rows


"""LJH file header information extracted from the ASCII file header."""
immutable LJHHeader
    filename         ::AbstractString
    version          ::Symbol
    nPresamples      ::Int64
    nSamples         ::Int64
    timebase         ::Float64
    timestampOffset  ::Float64
    date             ::AbstractString
    headerSize       ::Int
    channum          ::Int16
    column           ::Int16
    row              ::Int16
    num_columns      ::Int16
    num_rows         ::Int16
end


"""LJH file abstraction"""
type LJHFile
    name             ::AbstractString        # filename
    str              ::IOStream      # IOStream to read from LJH file
    header           ::LJHHeader     # LJH file header data
    version          ::Symbol
    nrec             ::Int64         # number of (pulse) records in file
    dt               ::Float64       # sample spacing (microseconds)
    npre             ::Int64         # nPresample
    nsamp            ::Int64         # number of sample per record
    reclength        ::Int64         # record length (bytes) including timestamp
    channum          ::Int16         # channel number
    column           ::Int16
    row              ::Int16
    num_columns      ::Int16
    num_rows         ::Int16
    function LJHFile(name::AbstractString)
        hd = readljhheader(name)
        dt = hd.timebase
        pre = hd.nPresamples
        tot = hd.nSamples
        str = open(name)
        datalen = stat(name).size - hd.headerSize
        reclen = 2*tot + 6
        if hd.version == :LJH_22
            reclen = 2*tot + 16
        end
        # assert((datalen%reclen)==0)
        nrec = div(datalen,reclen)
        seek(str,hd.headerSize)
        new(name, str, hd, hd.version, nrec, dt, pre, tot, reclen,
            hd.channum, hd.column, hd.row, hd.num_columns, hd.num_rows)
    end
end

function update_num_records(f::LJHFile)
    datalen = stat(f.name).size - f.header.headerSize
    f.nrec = div(datalen,f.reclength)
end
Base.close(f::LJHFile) = close(f.str)
Base.open(f::LJHFile) = f.str=open(f.name)


"""Represent a slice from an LJH file"""
type LJHSlice{T<:AbstractArray}
    ljhfile::LJHFile
    slice::T
    function LJHSlice(ljhfile, slice)
        maximum(slice)<=ljhfile.nrec || error("$(maximum(slice)) is greater than nrec=$(ljhfile.nrec) in $ljhfile")
        new(ljhfile, slice)
    end
end
LJHSlice{T<:AbstractArray}(ljhfile::LJHFile, slice::T) = LJHSlice{T}(ljhfile, slice)


function Base.show(io::IO, f::LJHFile)
    print(io, "LJHFile channel $(f.channum): $(f.nrec) records with $(f.npre) presamples and $(f.nsamp) samples each\n")
    print(io, f.name*"\n")
end



"""Get pertinent information from LJH file header and return it as an LJHHeader object."""
function readljhheader(filename::AbstractString)
    open(filename) do str # ensures str is closed
    labels=Dict("base"   =>"Timebase:",
            "date"   =>"Date:",
            "date1"  =>"File First Record Time:",
            "end"    =>"#End of Header",
            "offset" =>"Timestamp offset (s):",
            "pre"    =>"Presamples: ",
            "tot"    =>"Total Samples: ",
            "version"=>"Save File Format Version: ",
            "channum"=>"Channel: ",
            "column"=>r"Column number .*: (\d+)",
            "row"=>r"Row number .*: (\d+)",
            "num_columns"=>"Number of columns:",
            "num_rows"=>"Number of rows:")
    nlines=0
    maxnlines=100
    column,row, num_columns, num_rows = -1,-1,-1,-1
    date = "unknown" # If header standard for date labels changes, we don't want a hard error

    # Read channel # from the file name, then update that result from the header, if it exists.
    m = match(r"_chan\d+", filename)
    if m == nothing
        channum = -1
    else
        channum = parse(Int16,m.match[6:end])
    end

    while nlines<maxnlines
        line=readline(str)
        nlines+=1
        if startswith(line,labels["end"])
            headerSize = position(str)
            return(LJHHeader(filename,version,nPresamples,nSamples,
                             timebase,timestampOffset,date,headerSize,channum,column,row,num_columns,num_rows))
        elseif startswith(line,labels["base"])
            timebase = parse(Float64,line[1+length(labels["base"]):end])
        elseif startswith(line,labels["date"]) # Old LJH files
            date = line[7:end-2]
        elseif startswith(line,labels["date1"])# Newer LJH files
            date = line[25:end-2]
        elseif startswith(line,labels["channum"])# Newer LJH files
            channum = Int16(parse(line[10:end]))
        elseif startswith(line,labels["offset"])
            timestampOffset = parse(Float64,line[1+length(labels["offset"]):end])
        elseif startswith(line,labels["pre"])
            nPresamples = parse(Int64,line[1+length(labels["pre"]):end])
        elseif startswith(line,labels["tot"])
            nSamples = parse(Int64,line[1+length(labels["tot"]):end])
        elseif startswith(line, labels["version"])
            version_str = line[1+length(labels["version"]):end]
            if startswith(version_str, "2.1.0")
                version = :LJH_21
            elseif startswith(version_str, "2.2.0")
                version = :LJH_22
            else
                error("read_LJH_header: version '$(version_str)' unknown.")
            end
        elseif ismatch(labels["column"],line)
            m=match(labels["column"],line)
            column = parse(Int64,m.captures[1])
        elseif ismatch(labels["row"],line)
            m=match(labels["row"],line)
            row = parse(Int64,m.captures[1])
        elseif startswith(line, labels["num_columns"])
            num_columns = parse(Int64,line[1+length(labels["num_columns"]):end])
        elseif startswith(line, labels["num_rows"])
            num_rows = parse(Int64,line[1+length(labels["num_rows"]):end])
        end
    end
    error("read_LJH_header: where's '$(labels["end"])' ?")
    end #do
end



"""Read the next nrec records and for each return time and samples
(error if eof occurs or insufficient space in data)."""
function extractljhdata(f::LJHFile, nrec::Integer, rows::Vector{Int64},
                     times::Vector{Int64}, data::Matrix{UInt16})
    #assert(nrec <= min(length(times),size(data,2)) && size(data,1)==f.nsamp)
    if f.version == :LJH_21
        for i=1:nrec
            rows[i], times[i] = record_row_count_v21(read(f.str, UInt8, 6), f.num_rows, f.row, f.dt)
            data[:,i] = read(f.str, UInt16, f.nsamp)
        end
    elseif f.version == :LJH_22
        for i=1:nrec
            rows[i] = read(f.str, Int64)
            times[i] = read(f.str, Int64)
            data[:,i] = read(f.str, UInt16, f.nsamp)
        end
    else
        error("Unknown LJH version number")
    end
end

"""From LJH file, return all data samples as single vector (dropping the
row counts and time stamps)."""
function ljhalldata(filename::AbstractString)
    ljh = LJHFile(filename)
    [d for (d,r,t) in ljh]
    close(ljh.str)
    data
end
ljhalldata(ljh::LJHFile) = [d for (d,r,t) in ljh]

"""Rewind LJHFile so that next record returned by iterating is the 1st record."""
ljhrewind(f::LJHFile) = seek(f.str, f.header.headerSize)


# support for ljhfile[1:7] syntax
seekto(f::LJHFile, i::Int) = seek(f.str,f.header.headerSize+(i-1)*f.reclength)
Base.getindex(f::LJHFile,indexes::AbstractVector)=LJHSlice(f, indexes)
function Base.getindex(f::LJHFile,index::Int)
    seekto(f, index)
    pop!(f)
end
function Base.pop!(f::LJHFile)
    if f.version == :LJH_21
        row_count, time_usec =  record_row_count_v21(read(f.str, UInt8, 6), f.num_rows, f.row, f.dt)
    elseif f.version == :LJH_22
        row_count = read(f.str, Int64)
        time_usec = read(f.str, Int64)
    end
    data = read(f.str, UInt16, f.nsamp)
    data, row_count, time_usec
end

Base.size(f::LJHFile) = (f.nrec,)
Base.length(f::LJHFile) = f.nrec
Base.endof(f::LJHFile) = f.nrec

# access as iterator
Base.start(f::LJHFile) = (seekto(f,1);1)
Base.next(f::LJHFile,j) = pop!(f),j+1
Base.done(f::LJHFile,j) = j==length(f)+1
Base.length(f::LJHFile) = f.nrec


Base.start{T}(f::LJHSlice{T}) = (j=start(f.slice);seekto(f.ljhfile,j);j)
function Base.next{T<:UnitRange}(f::LJHSlice{T}, j)
    n,r=next(f.slice,j)
    pop!(f.ljhfile),r
end
function Base.next{T}(f::LJHSlice{T},j)
    n,r=next(f.slice,j)
    f.ljhfile[n],r
end
Base.done{T}(f::LJHSlice{T}, j) = done(f.slice,j)
Base.length{T}(f::LJHSlice{T}) = length(f.slice)
Base.endof{T}(f::LJHSlice{T}) = length(f.slice)


"""Represent one or more LJHFiles as a seamless sequence that can be addressed
by record number from 1 to the sum of all records in the group."""
type LJHGroup
    ljhfiles::Vector{LJHFile}
    lengths::Vector{Int}
end
LJHGroup(x::Vector{LJHFile}) = LJHGroup(x, Int[length(f) for f in x])
LJHGroup(x) = LJHGroup([LJHFile(f) for f in x])
LJHGroup(x::LJHFile) = LJHGroup([x])
LJHGroup(x::AbstractString) = LJHGroup(LJHFile(x))
Base.length(g::LJHGroup) = sum(g.lengths)
Base.close(g::LJHGroup) = map(close, g.ljhfiles)
Base.open(g::LJHGroup) = map(open, g.ljhfiles)
fieldvalue(g::LJHGroup, s::Symbol) = unique([getfield(f, s) for f in g.ljhfiles])
channel(g::LJHGroup) = (assert(length(fieldvalue(g, :channum))==1);g.ljhfiles[1].channum)
record_nsamples(g::LJHGroup) = (assert(length(fieldvalue(g, :nsamp))==1);g.ljhfiles[1].nsamp)
pretrig_nsamples(g::LJHGroup) = (assert(length(fieldvalue(g, :npre))==1);g.ljhfiles[1].npre)
frametime(g::LJHGroup) = (assert(length(fieldvalue(g, :dt))==1);g.ljhfiles[1].dt)
column(g::LJHGroup) = (assert(length(fieldvalue(g, :column))==1);g.ljhfiles[1].column)
row(g::LJHGroup) = (assert(length(fieldvalue(g, :row))==1);g.ljhfiles[1].row)
num_columns(g::LJHGroup) = (assert(length(fieldvalue(g, :num_columns))==1);g.ljhfiles[1].num_columns)
num_rows(g::LJHGroup) = (assert(length(fieldvalue(g, :num_rows))==1);g.ljhfiles[1].num_rows)
filenames(g::LJHGroup) = convert(Array{ASCIIString},fieldvalue(g, :name))
lengths(g::LJHGroup) = g.lengths
function update_num_records(g::LJHGroup)
    old_lengths = copy(g.lengths)
    update_num_records(last(g.ljhfiles))
    g.lengths = Int[length(f) for f in g.ljhfiles]
    for i = 1:length(g.lengths)-1
        if old_lengths[i]!=g.lengths[i]
            error("a ljh file other than the last file in grew in length $g it was $(g.ljhfiles[i])")
        end
    end
end
function filenum_pulsenum(g::LJHGroup, j::Int)
    for (i,len) in enumerate(g.lengths)
        j <= len ? (return i,j) : (j-=len)
    end
    1,1 # default return value in case of empty range
end
function Base.getindex(g::LJHGroup, i::Int)
    filenum, pulsenum = filenum_pulsenum(g,i)
    g.ljhfiles[filenum][pulsenum]
end
Base.getindex(g::LJHGroup, slice::AbstractArray) = LJHGroupSlice(g, slice)
Base.endof(g::LJHGroup) = length(g)
immutable LJHGroupSlice{T<:AbstractArray}
    g::LJHGroup
    slice::T
    function LJHGroupSlice(ljhgroup, slice)
        isempty(slice) || maximum(slice)<=length(ljhgroup) || error("$(maximum(slice)) is greater than nrec=$(length(ljhgroup)) in $ljhgroup")
        new(ljhgroup, slice)
    end
end
Base.length(g::LJHGroupSlice) = length(g.slice)
Base.endof(g::LJHGroupSlice) = length(g.slice)
LJHGroupSlice{T<:AbstractArray}(ljhfile::LJHGroup, slice::T) = LJHGroupSlice{T}(ljhfile, slice)
function Base.start{T<:UnitRange}(g::LJHGroupSlice{T})
    for f in g.g.ljhfiles seekto(f,1) end
    isempty(g.slice) && return (2,2,1,1) # ensure done condition is immediatley met on empty range
    filenum, pulsenum = filenum_pulsenum(g.g, first(g.slice))
    donefilenum, donepulsenum = filenum_pulsenum(g.g, last(g.slice))
    seekto(g.g.ljhfiles[filenum], pulsenum)
    return (filenum, pulsenum, donefilenum, donepulsenum)
end
function Base.next{T<:UnitRange}(g::LJHGroupSlice{T}, state)
    filenum, pulsenum, donefilenum, donepulsenum = state
    data = pop!(g.g.ljhfiles[filenum])
    pulsenum+=1
    pulsenum > g.g.lengths[filenum] && (pulsenum-=g.g.lengths[filenum];filenum+=1)
    data, (filenum, pulsenum, donefilenum, donepulsenum)
end
function Base.done{T<:UnitRange}(g::LJHGroupSlice{T}, state)
    filenum, pulsenum, donefilenum, donepulsenum = state
    filenum>donefilenum || filenum==donefilenum && pulsenum>donepulsenum
end
function Base.collect{T<:UnitRange}(g::LJHGroupSlice{T})
    pulses = Array(Vector{UInt16}, length(g.slice))
    rowstamps = Array(Int, length(g.slice))
    for (i,(pulse, rowstamp)) in enumerate(g)
        pulses[i]=pulse
        rowstamps[i]=rowstamp
    end
    pulses, rowstamps
end

"""Used only for reading LJH version 2.1.0 files. This parses the ugly
"encoded" version of the frame counter, which is converted into an
approximate time, rounded to 4 microseconds, and divided into the integer
and fractional parts of the millisecond. The latter are stored in bytes
3-6, and the former is in byte 1. Ignore byte 2 (it has some need to be
0 for backward compatibility). Ugly! That's why LJH 2.2.0 does something
totally different.
"""
function record_row_count_v21(header::Vector{UInt8}, num_rows::Integer, row::Integer, frame_time::Float64)
    frac = Int64(header[1])
    ms = UInt64(header[3]) |
         (UInt64(header[4])<<8) |
         (UInt64(header[5])<<16) |
         (UInt64(header[6])<<24)
    count_4usec = Int64(ms*250+frac)
    ns_per_frame = round(Int64,frame_time*1e9)
    ns_per_4usec = Int64(4000)
    count_nsec = count_4usec*ns_per_4usec
    count_frame = cld(count_nsec,ns_per_frame)
    if row == -1 # stupid workaround for the fact that -1 ruins the row calculation
        row = 0
        num_rows = 30
    end
    count_row = count_frame*num_rows+row
    return count_row, 4*count_4usec
end


# writing ljh files

"""Write a header for an LJH file."""
function writeljhheader(filename::AbstractString, dt, npre, nsamp; version="2.2.0")
    open(filename, "w") do f
    writeljhheader(f, dt, npre, nsamp; version=version)
    end #do
end

function writeljhheader(io::IO, dt, npre, nsamp; version="2.2.0")
    write(io,
"#LJH Memorial File Format
Save File Format Version: $(version)
Software Version: Fake LJH file generated by Julia
Software Driver Version: n/a
Date: %(asctime)s GMT
Acquisition Mode: 0
Digitized Word Size in bytes: 2
Location: LANL, presumably
Cryostat: Unknown
Thermometer: Unknown
Temperature (mK): 100.0000
Bridge range: 20000
Magnetic field (mGauss): 100.0000
Detector:
Sample:
Excitation/Source:
Operator: Unknown
SYSTEM DESCRIPTION OF THIS FILE:
USER DESCRIPTION OF THIS FILE:
#End of description
Number of Digitizers: 1
Number of Active Channels: 1
Timestamp offset (s): 1304449182.876200
Digitizer: 1
Description: CS1450-1 1M ver 1.16
Master: Yes
Bits: 16
Effective Bits: 0
Anti-alias low-pass cutoff frequency (Hz): 0.000
Timebase: $(dt)
Number of samples per point: 1
Presamples: $(npre)
Total Samples: $(nsamp)
Trigger (V): 250.000000
Tigger Hysteresis: 0
Trigger Slope: +
Trigger Coupling: DC
Trigger Impedance: 1 MOhm
Trigger Source: CH A
Trigger Mode: 0 Normal
Trigger Time out: 351321
Use discrimination: No
Channel: 1.0
Description: A (Voltage)
Range: 0.500000
Offset: -0.000122
Coupling: DC
Impedance: 1 Ohms
Inverted: No
Preamp gain: 1.000000
Discrimination level (%%): 1.000000
#End of Header\n"
    )
end

"""Write LJH file data to an IO object or a filename given a AbstractString."""
function writeljhdata(filename::AbstractString, a...)
    open(filename, "a") do f
    writeljhdata(f,a...)
    end
end

# Write LJH v2.2+ data, with row # and timestamps
function writeljhdata(io::IO,traces::Array{UInt16,2}, rows::Vector{Int64}, times::Vector{Int64})
    for j = 1:length(times)
        writeljhdata(io, traces[:,j], rows[j], times[j])
    end
end
function writeljhdata(io::IO, trace::Vector{UInt16}, row::Int64, time::Int64)
    write(io, row)
    write(io, time)
    write(io, trace)
end


# Write LJH v2.1 data, with row # but no timestamps
function writeljhdata(io::IO,traces::Array{UInt16,2}, rows::Vector{Int64})
    for j = 1:length(rows)
        writeljhdata(io, traces[:,j], rows[j])
    end
end
function writeljhdata(io::IO, trace::Vector{UInt16}, row::Int64)
    timestamp_us = round(Int32, row*0.32) # made-up line rate of 320 nanoseconds per row.
    timestamp_ms = Int32(div(timestamp_us, 1000))
    subms_part = round(UInt8, mod(div(timestamp_us,4), 250))
    dummy_channum = Int8(0)
    write(io, subms_part)
    write(io, dummy_channum)
    write(io, timestamp_ms)
    write(io, trace)
end

end # end module


"""Functions for working with LJH file names and the MATTER sentinel file. Intended primarily to easily find all LJH files from the same run,
and to enable compatiblity with Python Mass."""
module LJHUtil

function ljhsplit(ljhname::AbstractString)
    if isdir(ljhname)
        dname = dirname(ljhname) # removes trailing /
        bname = last(split(dname,'/'))
        return dname, bname, ".ljh"
    end
    bname,ext = splitext(basename(ljhname))
    ext = isempty(ext) ? ".ljh" : ext
    m = match(r"_chan\d+", bname)
    dirname(ljhname), m == nothing ? bname : bname[1:m.offset-1], ext
end
function channel(ljhname::AbstractString)
    m = match(r"_chan(\d+)", ljhname)
    m == nothing ? -1 : parse(Int,m.captures[1])
end
function fnames(ljhname::AbstractString, chans)
    dname, bname, ext = ljhsplit(ljhname)
    [joinpath(dname, "$(bname)_chan$c$ext") for c in chans]
end
function fnames(ljhname::AbstractString, c::Int)
    dname, bname, ext = ljhsplit(ljhname)
    joinpath(dname, "$(bname)_chan$c$ext")
end
function allchannels(ljhname::AbstractString)
    dname, bname, ext = ljhsplit(ljhname)
    potential_ljh = filter!(s->(startswith(s,bname)), readdir(dname))
    channels = filter!(x->x>=0,[channel(p) for p in filter!(s->startswith(s,bname),readdir(dname))])
    sort!(unique(channels))
end
ljhall(ljhname::AbstractString) = ljhfnames(ljhname, ljhallchannels(ljhname))
function hdf5_name_from_ljh(ljhnames::AbstractString...)
	dname, bname, ext = ljhsplit(ljhnames[1])
	fname = prod([split(f)[2] for f in ljhnames])
	joinpath(dname,hdf5_name_from_ljh(fname))
end
hdf5_name_from_ljh(ljhname::AbstractString) = joinpath(ljhsplit(ljhname)...)*"_mass.hdf5"

# returns (AbstractString, Bool) representing (filename, currently_open?)
const sentinel_file_path = joinpath(expanduser("~"),".daq","latest_ljh_pulse.cur")
function matter_writing_status()
    isfile(sentinel_file_path) || error("$(sentinel_file_path) must be a file")
    open(sentinel_file_path,"r") do sentinel_file
    lines = map(chomp,collect(eachline(sentinel_file)))
    return lines[1], length(lines)==1 # the sentinel file has a second line that says closed when it has closed a file
    end
end

function write_sentinel_file(filename, writingbool)
    open(sentinel_file_path,"w") do sentinel_file
        println(sentinel_file,filename)
        if writingbool == true
            println(sentinel_file, "closed")
        end
    end
end

function change_writing_status()
    fname, writingbool = matter_writing_status()
    write_sentinel_file(fname, !writingbool)
end

end # module

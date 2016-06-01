using RunningVectors
import JLD: JldGroup, JldFile, JldDataset, HDF5Dataset
import Base: keys, length, values, getindex, setindex!, get, haskey
const DONETHRU_MAX = typemax(Int)-1

abstract AbstractStep


type MassChannel
	d::Dict{Symbol,Any} # contains pulse related data
	# workdone_cumulative[i] represents how much work steps[i] has done
	# same idea for time_elapsed_cumulative and workdone_last
	steps::Vector{AbstractStep}
	nextstepind::Int #used for debugging, tells which step it crashed on
	workdone_cumulative::Vector{Int}
	time_elapsed_cumulative::Vector{Float64}
	workdone_last::Vector{Int}
	perpulse_symbols::Set{Symbol}
	task::Nullable{Task}
	waittask::Nullable{Task} # used to supress printing of error message in task
	endertask::Nullable{Task} # schedule this to end task when all(workdone_last.==0)
	exitchannel::Channel{Int}
	oncleanfinish::Function
	graph
end
"Create an empty MassChannel"
MassChannel() = MassChannel(Dict{Symbol,Any}(),Array{AbstractStep}(0),0,Vector{Int}(),Vector{Int}(),Vector{Int}(),
  							Set{Symbol}(), Nullable{Task}(), Nullable{Task}(), Nullable{Task}(),Channel{Int}(1), (mc)->nothing, nothing)
# delegate a few functions to the dict d, to enable mc[:pulse] = blah
for f in (:keys, :length, :values)
	@eval $f(mc::MassChannel) = $f(mc.d)
end
getindex(mc::MassChannel, i) = getindex(mc.d,i)
setindex!(mc::MassChannel, i...) = setindex!(mc.d, i...)
get(mc::MassChannel, i...) = get(mc.d, i...)
haskey(mc::MassChannel,key) = haskey(mc.d,key)


"""PerPulseStep is for things like filtering data, where a perpulse input (say the actual pulse record) is transformed to a perpulse
output (say the filtered pulse height, aka filt_value). For scalars this is easily expressed as
`filt_value = do_filter(filter, pulse_record)`. To express this same thing as a step you would use
`PerPulseStep(do_filter, [:filter, :pulse_record], [:filt_value])` and be sure to add both `pulse_record` and `filt_value` to
`perpulse_symbols`."""
type PerPulseStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
	outputs::Vector{Symbol}
end
type ToJLDStep <:AbstractStep
	inputs::Vector{Symbol}
	one_time_inputs::Vector{Pair{Symbol, ASCIIString}}
	jldfilename::AbstractString
end
type HistogramStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
	outputs::Vector{Symbol} # this contains the Histogram which is also an input, should be redone someday, maybe the histogram can be a member
end
type ThresholdStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
	outputs::Vector{Symbol}
	to_watch::Symbol
	to_watch_func::Function
	threshold::Int
	do_if_able::Bool
end
type MockPulsesStep <: AbstractStep
	pg::TwoExponentialPulseGenerator{Int}
	pulses_per_step::Int
	max_pulses_ever::Int
	outputs::Vector{Symbol}
end
type FromHDF5Step <: AbstractStep # used to popuplate the mc dictionary with items from a jld/hdf5 file
	inputs::Vector{Symbol}
	one_time_inputs::Vector{Pair{Symbol, ASCIIString}}
	hdf5filename::AbstractString
end
type FreeMemoryStep <: AbstractStep
end
isperpulse(q::Symbol, c::MassChannel) = q in c.perpulse_symbols
getfunction(s::AbstractStep) = s.func
graphlabel(s::AbstractStep) = string(getfunction(s))
inputs(s::AbstractStep) = s.inputs
outputs(s::AbstractStep) = s.outputs
inputs(s::AbstractStep, c::MassChannel) = [c[q] for q in inputs(s)]
outputs(s::AbstractStep, c::MassChannel) = [c[q] for q in outputs(s)]
inputs(s::AbstractStep, c::MassChannel, r::UnitRange{Int}) =
[isperpulse(q,c)?c[q][r]:c[q] for q in s.inputs]
perpulse_inputs_key(s::AbstractStep, c::MassChannel) = inputs(s)[[isperpulse(q,c) for q in inputs(s)]]
perpulse_outputs_key(s::AbstractStep, c::MassChannel) = outputs(s)[[isperpulse(q,c) for q in outputs(s)]]
perpulse_inputs(s::AbstractStep, c::MassChannel) = [c[q] for q in perpulse_inputs_key(s,c)]
perpulse_outputs(s::AbstractStep, c::MassChannel) = [c[q] for q in perpulse_outputs_key(s,c)]
other_inputs_key(s::AbstractStep, c::MassChannel) = inputs(s)[[!isperpulse(q,c) for q in inputs(s)]]
other_outputs_key(s::AbstractStep, c::MassChannel) = outputs(s)[[!isperpulse(q,c) for q in outputs(s)]]
other_inputs(s::AbstractStep, c::MassChannel) = [c[q] for q in other_inputs_key(s,c)]
other_outputs(s::AbstractStep, c::MassChannel) = [c[q] for q in other_outputs_key(s,c)]
mindonethru(x) = length(x) == 0 ? DONETHRU_MAX : minimum(map(donethru,x))
Base.range(s::AbstractStep, c::MassChannel) = 1+mindonethru(perpulse_outputs(s,c)):mindonethru(perpulse_inputs(s,c))
other_inputs_exist(s::AbstractStep, c::MassChannel) = all([haskey(c,q) for q in other_inputs_key(s,c)])
function debug(s::AbstractStep, c::MassChannel)
	dump(s)
	println(range(s,c))
	println("perpulse inputs")
	for pi in perpulse_inputs_key(s,c)
		try
			println((pi, typeof(c[pi]), donethru(c[pi])))
		catch
			println("c[:$pi] doesn't exist")
		end
	end
	println("other inputs")
	for oi in other_inputs_key(s,c)
		try
			println((oi, typeof(c[oi]), donethru(c[oi])))
		catch
			println("c[:$oi] doesn't exist")
		end
	end
	println("perpulse outputs")
	for po in perpulse_outputs_key(s,c)
		try
			println((po, typeof(c[po]), donethru(c[po])))
		catch
			println("c[:$po] doesn't exist")
		end
	end
	println("other outputs")
	for oo in other_outputs_key(s,c)
		try
			println((oo, typeof(c[oo]), donethru(c[oo])))
		catch
			println("c[:$oo] doesn't exist")
		end
	end
end

# placeholder versions of exists, probably would use Nullable types here
function dostep!(s::PerPulseStep,c::MassChannel)
	f = getfunction(s)
	r = range(s,c)
	other_inputs_exist(s,c) || (return last(r):-1)
	length(r)==0 && (return r)
	fout = f(inputs(s,c,r)...) # it doesn't know the types here, tho f will, in principle one could write a macro so it knew types
	if length(outputs(s)) == 1 # fout will just be a value
		append!(outputs(s,c)[1],fout)
	elseif length(outputs(s)) > 1 # fout will be a tuple of values, one for each of outputs(s)
		@assert(length(fout) == length(outputs(s)))
		for (fo, co) in zip(fout, outputs(s,c))
			# fo is function out
			# co is Masschannel output eg the thing store in the Masschannel object
			append!(co,fo)
			assert(length(co)==last(r)) # append isn't the right choice if the lengths differ
			# I'm not sure, but there may be no good reason to support differing lengths
		end
	end
	# freeuntil needs information on other dependencies
	r
end

# HistogramStep outputs should contain 1 symbol which points to a Histogram
# the dostep function will use this as both the first input, and the output
Base.range(s::HistogramStep, c::MassChannel) = 1+mindonethru(outputs(s,c)):mindonethru(perpulse_inputs(s,c))
function inputs(s::HistogramStep, c::MassChannel, r::Range)
	out = Any[]
	for inp in inputs(s,c)
		push!(out, typeof(inp) <: AbstractRunningVector ? inp[r] : inp)
	end
	out
end
function dostep!(s::HistogramStep, c::MassChannel)
	f = getfunction(s)
	r = range(s,c)
	other_inputs_exist(s,c) || (return last(r):-1)
	histogram = outputs(s,c)[1]
	inps = inputs(s,c,r)
	finps = (histogram, inps...)
	donethru_before = donethru(histogram)
	fout = f(finps...) # it doesn't know the types here, tho f will, also I could write a macro so it knew types
	workdone = donethru(histogram)-donethru_before
end

function dostep!(s::ThresholdStep, c::MassChannel)
	s.do_if_able || (return false)
	other_inputs_exist(s,c) || (return false)
	n_other = mindonethru(perpulse_outputs(s,c))
	n = s.to_watch_func(c[s.to_watch])
	if n >= s.threshold && n_other >= s.threshold
		f = getfunction(s)
		r=1:s.threshold
		fout = f(inputs(s,c,r)...)
		if length(outputs(s)) == 1 # fout will be a value
			c[outputs(s)[1]] =  fout
		elseif length(outputs(s)) > 1
			@assert(length(fout) == length(outputs(s)))
			for (single_output, output_key) in zip(fout, outputs(s))
				c[output_key]=single_output
			end
		end
		s.do_if_able = false
		return true
	else
		return false
	end
end


graphlabel(s::MockPulsesStep) = "MockPulsesStep"
inputs(s::MockPulsesStep) = []
function dostep!(s::MockPulsesStep, c::MassChannel)
	n_pulses_to_get = min(s.pulses_per_step, s.max_pulses_ever-length(c[s.outputs[1]]))
	n_pulses_to_get<=0 && return 0
	pulses, rowstamps = gettriggeredpulse!(s.pg,n_pulses_to_get)
	append!(c[s.outputs[1]], pulses)
	append!(c[s.outputs[2]], rowstamps)
	s.pulses_per_step
end

graphlabel(s::FreeMemoryStep) = "FreeMemoryStep"
inputs(s::FreeMemoryStep) = Symbol[]
outputs(s::FreeMemoryStep) = Symbol[]
function dostep!(s::FreeMemoryStep, mc::MassChannel)
	n_freed = 0
	for q in mc.perpulse_symbols
		d=mc[q]
		l0 = length(available(d))
		freeuntil!(d,min(earliest_needed_index(mc,q,mc.graph)-1,length(d)))
		n_freed += l0-length(available(d))
	end
	n_freed # steps must return an int that increases with amount of work done, and is zero when no work is done
end

## Get Pulses Step
import LJH.LJHGroup
type GetPulsesStep{T} <: AbstractStep
	pulse_source::T
	outputs::Vector{Symbol}
	previous_pulse_index::Int # used to keep track of location in an ljh file for example
	max_pulses_per_step::Int
	last_timestamp_checked::Int
end
GetPulsesStep(a,b,c,d) = GetPulsesStep(a,b,c,d,0)
inputs(s::GetPulsesStep) = Symbol[]
posix_time() = convert(Int, time()*1000000)
function dostep!(s::GetPulsesStep{LJHGroup},c::MassChannel)
	LJH.update_num_records(s.pulse_source)
	r = s.previous_pulse_index+1:min(s.previous_pulse_index+s.max_pulses_per_step, length(s.pulse_source))
	length(r)==0 && (return r)
	pulses, rowcounts, timestamps = get_data_rowcount_timestamp(s.pulse_source[r])
	if last(r) == length(s.pulse_source) # determine when was the latest time this would have reported a pulse from
		s.last_timestamp_checked = posix_time()
	else
		s.last_timestamp_checked = last(timestamps)
	end
	pulses_out, rowcounts_out, timestamps_out = perpulse_outputs(s,c)
	append!(pulses_out, pulses)
	assert(length(pulses_out)==last(r))
	append!(rowcounts_out, rowcounts)
	assert(length(rowcounts_out)==last(r))
	append!(timestamps_out, timestamps)
	assert(length(timestamps_out)==last(r))
	s.previous_pulse_index=last(r)
	c[:last_timestamp_checked]=s.last_timestamp_checked # THIS IS A HACK, FIX IT, SHOULD NOT DIRECTLY SET VARIABLES
	r
end
graphlabel(s::GetPulsesStep) = repr(typeof(s))


type MemoryLimitStep <: AbstractStep
	maxbytes::Int
end
graphlabel(s::MemoryLimitStep) = "MemoryLimitStep"
inputs(s::MemoryLimitStep) = Symbol[]
outputs(s::MemoryLimitStep) = Symbol[]
function dostep!(s::MemoryLimitStep, c::MassChannel)
	usedbytes = Base.summarysize(c)
	usedbytes > s.maxbytes && error("MemoryLimitStep with maxbytes=$(s.maxbytes) noticed MassChannel using $usedbytes")
	0
end



# JLD/HDF5 helper functions for ToJLDStep
function d_extend(d::HDF5Dataset, value::Vector, range::UnitRange)
	set_dims!(d, (maximum(range),))
	d[range] = value
	d
end
d_extend(d::JldDataset, value::Vector, range::UnitRange) = d_extend(d.plain, value, range)
function d_require(parent::Union{JldFile,JldGroup}, name, elementype::Type ,chunksize = 1000)
	dims = ((1,), (-1,)) # create a minimum size 1d dataset with largest possible maximum dimension, zero is not allowed
	exists(parent,name) ? parent[name] : d_create(parent.plain, name, elementype, dims, "chunk", (chunksize,))
end
function update!(parent::Union{JldFile,JldGroup}, name, value)
	exists(parent, name) && delete!(parent[name])
	parent[name]=value
end
graphlabel(s::ToJLDStep) = "to JLD"
outputs(s::ToJLDStep) = []
function dostep!(s::ToJLDStep,c::MassChannel)
	n_written = 0
	jldopen(s.jldfilename,isfile(s.jldfilename) ? "r+" : "w") do jld
		for (sym,value) in zip(perpulse_inputs_key(s,c),perpulse_inputs(s,c))
			length(value) == 0 && continue #don't create 1 entry datasets with no meaninful data in them
			d=d_require(jld, string(sym), eltype(value))
			# account for the fact that the dataset is created with 1 element, not zero
			start = length(d)==1 ? 1 : length(d)+1
			r = start:length(value)
			n_written += length(r)
			if length(r)>0
				d_extend(d, value[r],r)
			end
		end
		for (sym, path) in s.one_time_inputs
			has(jld.plain, path) && continue # if something is already written to the path, move on
			haskey(c,sym) || continue # if the symbol doesn't exist in the masschannel dict, move on
			jld.plain[path]=c[sym] # otherwise write to file, write to .plain to avoid any julia specific formatting
		end
	end
	# JLD.delete! doesn't work well enough to to this, it can leave spare references around
	# for (sym, value) in  zip(other_inputs_key(s,c),other_inputs(s,c))
	# 	# warn("doesn't work if you do it more than once")
	# 	# update!(jld, string(sym), value)
	# end
	n_written # return workunits
end #dostep!

# FromHDF5Step
graphlabel(s::FromHDF5Step) = "from HDF5"
outputs(s::FromHDF5Step) = [sym for (sym, path) in s.one_time_inputs]
function dostep!(s::FromHDF5Step, c::MassChannel)
	n_gotten = 0
	!isfile(s.hdf5filename) && (return 0)
	h5open(s.hdf5filename, "r") do h5
		for (sym, value) in zip(perpulse_inputs_key(s,c),perpulse_inputs(s,c))
			error("perpulse inputs not yet support in FromHDF5Step")
		end
		for (sym, path) in s.one_time_inputs
			has(h5, path) || continue # if the path doesn't exist in the hdf5 file, move one
			haskey(c,sym) && continue # if the symbol already exists in the masschannel dict, move on
			c[sym] = read(h5[path])
			n_gotten+=1
		end
	end
		n_gotten
end

default_vertex_attrs = Graphs.AttributeDict()
default_vertex_attrs["fontsize"]=18.0
default_vertex_attrs["color"]="red"
default_vertex_attrs["fillcolor"]="lightgrey"
default_vertex_attrs["style"]="filled"
label(v)=v.label
vertex_labels(g) = [label(v) for v in Graphs.vertices(g)]
function vertex(g,label)
	i=first(indexin([label], vertex_labels(g)))
	i==0 && error("no vertex with label $label in $v")
	Graphs.vertices(g)[i]
end
add_exvertex!(g,name::Symbol, attrs_in) = add_exvertex!(g,string(name), attrs_in)
function add_exvertex!(g, name,attrs_in)
    verts = Graphs.vertices(g)
    vert_labels = [label(v) for v in verts]
    if !(name in vert_labels)
        v = Graphs.ExVertex(Graphs.num_vertices(g)+1, name)
        attrs = Graphs.attributes(v,g)
        attrs["label"]=name
        for (key,val) in default_vertex_attrs
        	attrs[key]=val
        end
        for (key,val) in attrs_in
        	attrs[key]=val
        end
        Graphs.add_vertex!(g, v)
        v
    else
        first(verts[indexin([name],vert_labels)])
    end
end
add_perpulse_data!(g,name)=add_exvertex!(g,name,[("shape","doubleoctagon"),("type","perpulse_data")])
add_other_data!(g,name)=add_exvertex!(g,name,[("shape","octagon"),("type","other_data")])
add_func!(g,name)=add_exvertex!(g,name,[("shape","box"),("type","function")])
function Graphs.add_edge!(g,label1, label2)
    verts = Graphs.vertices(g)
    vert_labels = [label(v) for v in verts]
    i1 = first(indexin([label1], vert_labels))
    i2 = first(indexin([label2], vert_labels))
    if i1 < 1 || i2 < 1
        error("error label $label1 or $label2 is not in $g")
    end
    add_edge!(g, verts[i1], verts[i2])
end
function add_step!(g, c::MassChannel,stepind::Int)
		stepind
		s = c.steps[stepind]
		# prepend the step index to each label to make sure different
		# steps have different verticies even if they call the same funcion
    vf = add_func!(g,"$stepind "*graphlabel(s))
    for p in perpulse_inputs_key(s,c)
        v=add_perpulse_data!(g,p)
        Graphs.add_edge!(g,v,vf)
    end
    for o in other_inputs_key(s,c)
        v=add_other_data!(g,o)
        Graphs.add_edge!(g,v,vf)
    end
    for p in perpulse_outputs_key(s,c)
        v=add_perpulse_data!(g,p)
        Graphs.add_edge!(g,vf,v)
    end
    for o in other_outputs_key(s,c)
        v=add_other_data!(g,o)
        Graphs.add_edge!(g,vf,v)
    end
    g
end
"graph(steps::Vector{AbstractStep},c::MassChannel)
Generate a directed graph representing the steps in `steps`. You must pass
`c` to allow it to figure out which things are per_pulse and which aren't."
function graph(steps::Vector{AbstractStep},c::MassChannel)
	oldg=Graphs.inclist(Graphs.ExVertex, is_directed=true)
	for i in eachindex(steps)
    	add_step!(oldg,c,i)
	end
	# graph is a type with a larger and more convenient API
	# I should have constructed it directly, but I didn't know about it, so
	# I'm just going to construct it from oldg for now
	g = Graphs.graph(Graphs.vertices(oldg),vcat(oldg.inclist...),is_directed=true)
	g
end
function checkforunuseddata(g::Graphs.AbstractGraph)
	verts = Graphs.vertices(g)
	for v in verts
		contains(v.attributes["type"],"data") || continue # only check for unused data, not functions or other verticies
		if isempty(Graphs.out_neighbors(v,g)) && isempty(Graphs.in_neighbors(v,g))
			error("vertex $v is a data vertex that is unused")
		end
	end
end
function warn_of_cycles(g)
	verts = Graphs.vertices(g)
	for v in verts
		visited = Graphs.visited_vertices(g,Graphs.BreadthFirst(),v)
		length(visited)<=1 && continue
		v in visited[2:end] && println("cycle containing $v")
	end
end
function savegraph(fname,g)
	dot = Graphs.to_dot(g)
	open("$fname.dot","w") do f
		write(f,dot)
		end
end
# the following code is correct and useful, but it conflicts with @pyimport
# so I'm commenting it until that is resolved
# gviz = GraphViz.Graph(dot)
# GraphViz.layout!(gviz,engine="dot")
# 	open("$fname.svg","w") do f
# 	       GraphViz.writemime(f, MIME"image/svg+xml"(),gviz)
# 	end #do
# 	open("$fname.png","w") do f
# 	       GraphViz.writemime(f, MIME"image/png"(),gviz)
# 	end #do


# figure out what we can free up
donethru(x::AbstractRunningVector) = length(x)
donethru(x::Vector{Float64}) = DONETHRU_MAX
donethru(x)=DONETHRU_MAX
earliest_needed_index(x) = donethru(x)+1

"earliest_needed_index(mc::MassChannel,q::Symbol,p::Symbol)
Internal function, don't call directly. Returns the earliest needed index
for `mc[q]`, or `mc[p]` if `q=:to_disk`."
function earliest_needed_index(mc::MassChannel,q::Symbol,p::Symbol)
	if q==:to_disk
		return donethru_jld(mc::MassChannel,q,p)+1
	elseif q in keys(mc)
			return earliest_needed_index(mc[q])
	else
		return 1 # if the item doesn't exist yet, it needs access to the first pulse
	end
end
function donethru_jld(c::MassChannel,q::Symbol,p::Symbol)
	jld = jldopen(filename(c),"r+")
	l = 0
	if string(p) in names(jld)
		l = length(jld[string(p)])
	end
	close(jld)
	l
end
"earliest_needed_index(mc::MassChannel, parent::Symbol, g::Graphs.AbstractGraph)
Return the earliest needed index for `mc[parent]`. `g` should be the graph
in a `FreeMemoryStep`."
function earliest_needed_index(mc::MassChannel, parent::Symbol, g::Graphs.AbstractGraph)
	v = vertex(g,string(parent))
	next_neighbor_vertices = graph_out_next_neighbors(v,g)
	children_sym = [symbol(label(u)) for u in next_neighbor_vertices]
	eni = [earliest_needed_index(mc,q,parent) for q in children_sym]
	#println([(q,parent,earliest_needed_index(mc,q,parent)) for q in children_sym])
	filter!(x->x != nothing, eni)
	isempty(eni)?length(mc[parent])+1:minimum(eni)
end

"graph_out_next_neighbors(v,g::Graphs.AbstractGraph) Find all the next nearest neighbors of vertex `v` in `g`. Used because
in Mass2 steps graph structure, the nearest neighbors of data are usually functions, and the nearest neighbors
of functions are data. Only data has `earliest_needed_index` defined, so we want to find the nearest
data neighbors of data when calculating what we can free."
function graph_out_next_neighbors(v,g::Graphs.AbstractGraph)
	# use Graphs.visited_vertices(g,Graphs.BreadthFirst(),v) to find all dependents
	children_vertices = Graphs.out_neighbors(v,g)
	next_neighbor_vertices = Any[]
	for vv in children_vertices
		next_neighbor_vertices=vcat(next_neighbor_vertices, collect(Graphs.out_neighbors(vv,g)))
	end
	next_neighbor_vertices
end



# for work reporting
workunits(x::Int) = x
workunits(r::Range) = length(r)
workunits(x::Bool) = convert(Int,x)
# for printing out work reports
workstat(n, s::MockPulsesStep, t) = "MockPulsesStep "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::PerPulseStep, t) = "PerPulse:$(graphlabel(s)) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ToJLDStep, t) = "ToJLDStep $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::HistogramStep, t) = "HistogramStep:$(inputs(s)[1]) "*@sprintf("%0.2f pulses/s",n/t)
workstat(n, s::ThresholdStep, t) = "ThresholdStep:$(graphlabel(s)) $n executions at "*@sprintf("%0.2f executions/s",n/t)
workstat(n, s::FreeMemoryStep, t) = "FreeMemoryStep $n executions at "*@sprintf("%0.2f executions/s",n/t)

#functions for analyings a vector/graph of steps to learn things such as:
# 1) list of perpulse_symbols
# 2) are there any required inputs that do not exist and will never exist? these must be pre-initialized or throw error
"Used to calculate which symbols must be perpulse_outputs. As of the first writing, this is
simply all of the outputs from `PerPulseStep`, `GetPulsesStep`, and `MockPulsesStep`."
perpulse_outputs_for_sure(s::Union{PerPulseStep, GetPulsesStep, MockPulsesStep}) = outputs(s)
perpulse_outputs_for_sure(s::AbstractStep) = AbstractStep[]
function perpulse_outputs_for_sure(steps::Vector{AbstractStep})
		perpulse_symbols = Set{Symbol}()
		for step in steps
			for o in perpulse_outputs_for_sure(step)
				o in perpulse_symbols && error("each perpulse symbol should only be an output of one step")
				push!(perpulse_symbols, o)
			end
		end
		perpulse_symbols
	end


# ## functions for repetetivley running many steps in a MassChannel
""" Write an hdf5 dataset with name `clean_exit_posix_timestamp_s` to mc[:hdf5_filename].
Used to indicate which channels succesfully finished analysis."""
function markhdf5oncleanfinish(mc::MassChannel)
	h5open(mc[:hdf5_filename],"r+") do h5
	h5["clean_exit_posix_timestamp_s"]=time()
	end
end

"""Continually call `dostep!(s,mc)` for `s` in `mc[:steps]`.
For each step measure `time_elapsed_cumulative`, `workdone_cumulative` and `workdone_last`.
Exits when `exitchannel` is ready (aka has a value `put!` into it). This should normally be done
by calling `plan_to_end(mc)`. Just before exiting it calls `mc.oncleanfinish(mc)`."""
function repeatsteps!(mc::MassChannel)
	while !isready(mc.exitchannel)
		for (i,s) in enumerate(mc.steps)
			yield()
			mc.nextstepind=i
			time_elapsed = @elapsed workdone = workunits(dostep!(s,mc))
			mc.workdone_cumulative[i]+=workdone
			mc.time_elapsed_cumulative[i]+=time_elapsed
			mc.workdone_last[i] = workdone
		end
	end
	mc.oncleanfinish(mc)
end

"""Creates prep `mc.task`, which is a task to run the steps in `mc`. Call `schedule(mc)` to schedule the task.
Preps `mc.waittask`, which is used to suppress automatic printing of error message.
Preps `mc.endertask`, which will end `mc.task` when all steps report having done no work on the last attempt.
Use `plan_to_end(mc)` to schedule `mc.endertask`, and to end `task`. You can `put!(mc.exitchannel,1)` to end all tasks
manually if needed, but it is not preffered."""
function preptasks!(mc::MassChannel)
	@assert isnull(mc.task) "trying to prep_tasks on MassChannel that has already been prepped"
	@assert isnull(mc.waittask) "trying to prep_tasks on MassChannel that has already been prepped"
	@assert isnull(mc.endertask) "trying to prep_tasks on MassChannel that has already been prepped"
	mc.task = Nullable(@task repeatsteps!(mc))
	mc.waittask = Nullable(@schedule try wait(mc.task) catch ex isa(ex, InterruptException) && throw(ex) end) # suppress printing of error message
	mc.endertask = Nullable(@task end_when_all_steps_do_no_work(mc.workdone_last, mc.exitchannel, mc.task.value))
end
"""Set `mc.steps=steps` and take care of housecleaning to make sure everything else makes sense. This at least
initializes `workdone_cumulative`, `workdone_last`, `time_elapsed_cumulative` to zero values and
calculates `perpulse_symbols` and `graph`."""
function setsteps!(mc::MassChannel, steps::Vector{AbstractStep})
	@assert length(mc.steps)==0 "using setsteps! on a MassChannel that already has steps seems like a bad idea, try addstep!"
	mc.steps = steps
	mc.perpulse_symbols = perpulse_outputs_for_sure(mc.steps)
	updateandverifygraph(mc)
	mc.workdone_cumulative = zeros(Int, length(steps))
	mc.workdone_last = zeros(Int, length(steps))
	mc.time_elapsed_cumulative = zeros(Float64, length(steps))
end
"""Update `mc.graph` and check for obvious problems. Looks for cycles in the graph, and
looks for totally disconnected data."""
function updateandverifygraph(mc::MassChannel)
	mc.graph = graph(mc.steps,mc)
	if Graphs.test_cyclic_by_dfs(mc.graph)
		warn_of_cycles(mc.graph)
		error("graph contains cycle, it shouldn't")
	end
	checkforunuseddata(mc.graph)
end
"""Add a step to  `mc.steps=steps` and take care of housecleaning to make sure everything else makes sense. This at least
adds new fields with 0 value to `workdone_cumulative`, `workdone_last`, `time_elapsed_cumulative` and
calculates `perpulse_symbols` and `graph`."""
function addstep!(mc::MassChannel, step::AbstractStep)
	push!(mc.steps,step)
	push!(mc.workdone_cumulative,0)
	push!(mc.workdone_last,0)
	push!(mc.time_elapsed_cumulative,0)
	mc.perpulse_symbols = perpulse_outputs_for_sure(mc.steps)
	updateandverifygraph(mc)
	@assert length(mc.steps)==length(mc.workdone_cumulative)==length(mc.workdone_last)==length(mc.time_elapsed_cumulative) "These things always have to have the same length"
end

"""As long as any step has nonzero `workdone_last`, yield. If all steps have zero `workdone_last` put `1` into `exitchannel`
to cause `task` to end. Exits if `istaskdone(task)`, to account for errors in `task`. """
function end_when_all_steps_do_no_work(workdone_last, exitchannel, task)
	while !all(workdone_last.==0) && !istaskdone(task)
		yield()
	end
	!isready(exitchannel) && put!(exitchannel,1)
	return nothing
end

Base.schedule(mc::MassChannel) = schedule(mc.task.value)
plantoend(mc::MassChannel) = schedule(mc.endertask.value) # need the .value since I'm using Nullable tasks
function debug(mc)
	println("failed on step index $(mc.nextstepind)")
	debug(mc.nextstepind, mc)
	println("**** now trying to execute the failed step")
	dostep!(mc.nextstepind, mc)
end
dostep!(i, mc::MassChannel) = dostep!(mc.steps[i], mc)
debug(i, mc::MassChannel) = debug(mc.steps[i], mc)
Base.wait(mc::MassChannel) = wait(get(mc.task))

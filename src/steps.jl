using RunningVectors
import JLD: JldGroup, JldFile, JldDataset, HDF5Dataset

const DONETHRU_MAX = typemax(Int)-1

include("histogram.jl")


typealias Channel Dict{Symbol,Any}
const perpulse_symbols =  Set{Symbol}()
isperpulse(s::Symbol) = s in perpulse_symbols
filename(c::Channel) = joinpath(pwd(), "flowsimple_test.jld")






abstract AbstractStep

type PerPulseStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
	outputs::Vector{Symbol}
end
type ToJLDStep <:AbstractStep
	inputs::Vector{Symbol}
end	
type HistogramStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
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
	outputs::Vector{Symbol}
end
type FreeMemoryStep <: AbstractStep
end
getfunction(s::AbstractStep) = s.func
graphlabel(s::AbstractStep) = string(getfunction(s))
inputs(s::AbstractStep) = s.inputs
outputs(s::AbstractStep) = s.outputs
inputs(s::AbstractStep, c::Channel) = [c[q] for q in inputs(s)]
outputs(s::AbstractStep, c::Channel) = [c[q] for q in outputs(s)]
inputs(s::AbstractStep, c::Channel, r::UnitRange{Int}) = 
[isperpulse(q)?c[q][r]:c[q] for q in s.inputs]
perpulse_inputs(s::AbstractStep) = inputs(s)[[isperpulse(q) for q in inputs(s)]]
perpulse_outputs(s::AbstractStep) = outputs(s)[[isperpulse(q) for q in outputs(s)]]
perpulse_inputs(s::AbstractStep, c::Channel) = [c[q] for q in perpulse_inputs(s)]
perpulse_outputs(s::AbstractStep, c::Channel) = [c[q] for q in perpulse_outputs(s)]
other_inputs(s::AbstractStep) = inputs(s)[[!isperpulse(q) for q in inputs(s)]]
other_outputs(s::AbstractStep) = outputs(s)[[!isperpulse(q) for q in outputs(s)]]
other_inputs(s::AbstractStep, c::Channel) = [c[q] for q in other_inputs(s)]
other_outputs(s::AbstractStep, c::Channel) = [c[q] for q in other_outputs(s)]
mindonethru(x) = length(x) == 0 ? DONETHRU_MAX : minimum(map(donethru,x))
Base.range(s::AbstractStep, c::Channel) = 1+mindonethru(perpulse_outputs(s,c)):mindonethru(perpulse_inputs(s,c))
other_inputs_exist(s::AbstractStep, c::Channel) = all([haskey(c,q) for q in other_inputs(s)])
function debug(s::AbstractStep, c::Channel)
	dump(s)
	println(range(s,c))
	println("perpulse inputs")
	for pi in perpulse_inputs(s)
		try
			println((pi, typeof(c[pi]), donethru(c[pi])))
		catch
			println("c[:$pi] doesn't exist")
		end
	end
	println("other inputs")
	for oi in other_inputs(s)
		try
			println((oi, typeof(c[oi]), donethru(c[oi])))
		catch
			println("c[:$oi] doesn't exist")
		end
	end
	println("perpulse outputs")
	for po in perpulse_outputs(s)
		try 
			println((po, typeof(c[po]), donethru(c[po])))
		catch
			println("c[:$po] doesn't exist")
		end
	end
	println("other outputs")
	for oo in other_outputs(s)
		try 
			println((oo, typeof(c[oo]), donethru(c[oo])))
		catch
			println("c[:$oo] doesn't exist")
		end
	end
end


# placeholder versions of exists, probably would use Nullable types here
function dostep!(s::PerPulseStep,c::Channel)
	f = getfunction(s)
	r = range(s,c)
	other_inputs_exist(s,c) || (return last(r):-1)
	length(r)==0 && (return r)
	fout = f(inputs(s,c,r)...) # it doesn't know the types here, tho f will, also I could write a macro so it knew types
	if length(outputs(s)) == 1 # fout will just be a value
		append!(outputs(s,c)[1],fout)
	elseif length(outputs(s)) > 1 # fout will be a tuple of values, one for each of outputs(s)
		@assert(length(fout) == length(outputs(s)))
		for (fo, co) in zip(fout, outputs(s,c))
			# fo is function out
			# co is channel output eg the thing store in the channel object
			append!(co,fo)
			assert(length(co)==last(r)) # append isn't the right choice if the lengths differ
			# I'm not sure, but there may be no good reason to support differing lengths
		end
	end
	# freeuntil needs information on other dependencies
	r
end

outputs(s::HistogramStep) = other_inputs(s) # for now I'm conflating the idea of perpulse and inputs vs inputs/outputs
perpulse_outputs(s::HistogramStep) = Symbol[]
other_outputs(s::HistogramStep) = outputs(s)
function inputs(s::HistogramStep, c::Channel, r::Range)
	out = Any[]
	for inp in inputs(s,c)
		push!(out, isa(inp,Histogram)?inp:inp[r])
	end
	out
end
function Base.range(s::HistogramStep, c::Channel)
	# for now I'm conflating the idea of per_pulse and histogram on inputs
	start = mindonethru(other_inputs(s,c))+1
	stop = mindonethru(perpulse_inputs(s,c))
	start:stop
end
function dostep!(s::HistogramStep, c::Channel)
	f = getfunction(s)
	r = range(s,c)
	other_inputs_exist(s,c) || (return last(r):-1)
	fout = f(inputs(s,c,r)...) # it doesn't know the types here, tho f will, also I could write a macro so it knew types
	r
end

function dostep!(s::ThresholdStep, c::Channel)
	s.do_if_able || (return false)
	other_inputs_exist(s,c) || (return false)
	n_other = mindonethru(perpulse_outputs(s,c))
	n = s.to_watch_func(c[s.to_watch])
	if n >= s.threshold && n_other >= s.threshold
		s.do_if_able = false
		f = getfunction(s)
		r=1:s.threshold
		fout = f(inputs(s,c,r)...)
		if length(outputs(s)) == 1 # fout will be a value
			c[outputs(s)[1]] =  fout
		elseif length(outputs(s)) > 1
			@assert(length(fout) == length(outputs(s)))
			for (fo, co) in zip(fout, outputs(s))
				c[fo]=co
			end
		end
		return true
	else
		return false
	end
end

graphlabel(s::MockPulsesStep) = "MockPulsesStep"
inputs(s::MockPulsesStep) = []
function dostep!(s::MockPulsesStep, c::Channel)
	pulses, rowstamps = gettriggeredpulse!(s.pg,s.pulses_per_step)
	append!(c[s.outputs[1]], pulses)
	append!(c[s.outputs[2]], rowstamps)
	s.pulses_per_step
end

graphlabel(s::FreeMemoryStep) = "FreeMemoryStep"
inputs(s::FreeMemoryStep) = Symbol[]
outputs(s::FreeMemoryStep) = Symbol[]
function dostep!(s::FreeMemoryStep, c::Channel)
	for q in perpulse_symbols
		d=c[q]
		freeuntil!(d,min(earliest_needed_index(q,c,g)-1,length(d)))
	end
	1 # return 1 unit of work done
end

## Get Pulses Step
import LJH.LJHGroup
type GetPulsesStep{T} <: AbstractStep
	pulse_source::T
	outputs::Vector{Symbol}
	previous_pulse_index::Int # used to keep track of location in an ljh file for example
	max_pulses_per_step::Int
end
inputs(s::GetPulsesStep) = Symbol[]
# placeholder versions of exists, probably would use Nullable types here
function dostep!(s::GetPulsesStep{LJHGroup},c::Channel)
	r = s.previous_pulse_index+1:min(s.previous_pulse_index+s.max_pulses_per_step, length(s.pulse_source))
	length(r)==0 && (return r)
	pulses, rowstamps = collect(s.pulse_source[r])
	pulses_out, rowstamps_out = perpulse_outputs(s,c)
	append!(pulses_out, pulses)
	assert(length(pulses_out)==last(r))
	append!(rowstamps_out, rowstamps)
	assert(length(rowstamps_out)==last(r))
	s.previous_pulse_index=last(r)
	r
end
graphlabel(s::GetPulsesStep) = repr(typeof(s))


# JLD/HDF5 helper functions for ToJLDStep
function d_extend(d::HDF5Dataset, value::Vector, range::UnitRange)
	set_dims!(d, (maximum(range),))
	d[range] = value
	d
end
d_extend(d::JldDataset, value::Vector, range::UnitRange) = d_extend(d.plain, value, range)
function d_require(parent::Union(JldFile,JldGroup), name, elementype::Type ,chunksize = 1000)
	dims = ((1,), (-1,)) # create a minimum size 1d dataset with largest possible maximum dimension, zero is not allowed
	exists(parent,name) ? parent[name] : d_create(parent.plain, name, elementype, dims, "chunk", (chunksize,))
end
function update!(parent::Union(JldFile,JldGroup), name, value)
	exists(parent, name) && delete!(parent[name])
	parent[name]=value
end
graphlabel(s::ToJLDStep) = "to JLD"
outputs(s::ToJLDStep) = []
function dostep!(s::ToJLDStep,c::Channel)
	jld = jldopen(filename(c),"r+")
	for (sym,value) in zip(perpulse_inputs(s),perpulse_inputs(s,c))
		d=d_require(jld, string(sym), eltype(value)) 
		# account for the fact that the dataset is created with 1 element, not zero 			
		start = length(d)==1 ? 1 : length(d)+1
		r = start:length(value)
		if length(r)>0
			println("writing $sym[$r] to $jld")
			d_extend(d, value[r],r)
		end
	end
	# JLD.delete! doesn't work well enough to to this, it can leave spare references around
	# for (sym, value) in  zip(other_inputs(s),other_inputs(s,c))
	# 	# warn("doesn't work if you do it more than once")
	# 	# update!(jld, string(sym), value)
	# end
	close(jld)
	1 # return workunits info
end #dostep!

default_vertex_attrs = AttributeDict()
default_vertex_attrs["fontsize"]=18.0
default_vertex_attrs["color"]="red"
default_vertex_attrs["fillcolor"]="lightgrey"
default_vertex_attrs["style"]="filled"
label(v)=v.label
vertex_labels(g) = [label(v) for v in vertices(g)]
function vertex(v,label)
	i=first(indexin([label], vertex_labels(g)))
	i==0 && error("no vertex with label $label in $v")
	vertices(g)[i]
end
add_exvertex!(g,name::Symbol, attrs_in) = add_exvertex!(g,string(name), attrs_in)
function add_exvertex!(g, name,attrs_in)
    verts = vertices(g)
    vert_labels = [label(v) for v in verts]
    if !(name in vert_labels)
        v = ExVertex(num_vertices(g)+1, name)
        attrs = attributes(v,g)
        attrs["label"]=name
        for (key,val) in default_vertex_attrs
        	attrs[key]=val
        end
        for (key,val) in attrs_in
        	attrs[key]=val
        end
        add_vertex!(g, v)
        v
    else
        first(verts[indexin([name],vert_labels)])
    end
end
add_perpulse_data!(g,name)=add_exvertex!(g,name,[("shape","doubleoctagon"),("type","perpulse_data")])
add_other_data!(g,name)=add_exvertex!(g,name,[("shape","octagon"),("type","other_data")])
add_func!(g,name)=add_exvertex!(g,name,[("shape","box"),("type","function")])
function Graphs.add_edge!(g,label1, label2)
    verts = vertices(g)
    vert_labels = [label(v) for v in verts]
    i1 = first(indexin([label1], vert_labels))
    i2 = first(indexin([label2], vert_labels))
    if i1 < 1 || i2 < 1
        error("error label $label1 or $label2 is not in $g")
    end
    add_edge!(g, verts[i1], verts[i2])
end
function add_step!(g, s)
    vf = add_func!(g,graphlabel(s))
    for p in perpulse_inputs(s)
        v=add_perpulse_data!(g,p)
        add_edge!(g,v,vf)
    end
    for o in other_inputs(s)
        v=add_other_data!(g,o)
        add_edge!(g,v,vf)
    end
    for p in perpulse_outputs(s)
        v=add_perpulse_data!(g,p)
        add_edge!(g,vf,v)
    end
    for o in other_outputs(s)
        v=add_other_data!(g,o)
        add_edge!(g,vf,v)
    end
    g
end
function Graphs.graph(steps::Vector{AbstractStep})
	g=Graphs.inclist(ExVertex, is_directed=true)
	for s in steps
    	add_step!(g,s)
	end
	g
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
function earliest_needed_index(c::Channel, q::Symbol,p::Symbol)
	if q in keys(c)
		return earliest_needed_index(c[q])
	elseif q==:to_disk
		return donethru_jld(c::Channel,q,p)+1
	end
end
function donethru_jld(c::Channel,q::Symbol,p::Symbol)
	jld = jldopen(filename(c),"r+")
	l = 0
	if string(p) in names(jld)
		l = length(jld[string(p)])
	end
	close(jld)
	l
end
function earliest_needed_index(parent::Symbol, c::Channel, g::AbstractGraph) 
	v = vertex(g,string(parent))
	children_vertices = visited_vertices(g,BreadthFirst(),v) # consider memoizing this
	shift!(children_vertices) # remove the first element, which is always v
	children_sym = [symbol(label(u)) for u in children_vertices]
	eni = [earliest_needed_index(c,q,parent) for q in children_sym]
	# println([(q,earliest_needed_index(c,q,parent)) for q in children_sym])
	filter!(x->x != nothing, eni)
	isempty(eni)?length(c[parent])+1:minimum(eni)
end


# for work reporting
workunits(x::Int) = x
workunits(r::Range) = length(r)
workunits(x::Bool) = convert(Int,x)
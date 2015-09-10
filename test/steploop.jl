# psuedo code for running steps on multiple channels in a single thread in fault tolerant way

function stepsloop!(channels::Vector{Channel}, steps::Vector{AbstractStep})
	while true
		try
			dosteps!(channels, steps)
		catch
			if isa(ex, InterruptException) # break on Ctrl-C for clean exit
				println("stepsloop ended by Ctrl-C")
				break
			end
		end
		if all([isdone(c,steps) for c in channels]) 
			println("stepsloop ended because all channels isdone=true")
			break
		end
	end
end

function dosteps!(channels::Vector{Channel})
	for c in channels
		try
			dosteps!(c)
		catch
			if isa(ex, InterruptException) # rethrow Ctrl-C
				rethrow(ex)
			end
		end
	end
end

function dosteps!(c::Channel, steps::Vector{AbstractStep})
	haskey(c, :exception_info)  && return # stop doing steps on this channel afte error
	workdone_cumulative = c[:workdone_cumulative]
	stepelapsed_cumulative = c[:stepelapsed_cumulative]
	for (i,s) in enumerate(steps)
		try
			tstart = time()
			wu = workunits(dostep!(c,s))
			stepelapsed = time()-tstart
			workdone_cumulative[j] += wu 
			stepelapsed_cumulative[j] += stepelapsed
		catch ex
			if isa(ex, InterruptException) # rethrow Ctrl-C
				rethrow(ex)
			else # store a backtrace of the exception
				ex_str = sprint(io->showerror(io, ex))) do io
					showerror(io, ex)
					print(io,"\n")
					Base.show_backtrace(io, catch_backtrace())
					end
				c[:exception_info] = ex_str
				return # stop doing steps on this channel afte error
			end
		end
	end
end
dosteps!(c::Channel) = dosteps(c, c[:steps])


#### TODO ####
 # switch channel to Dict{AbstractString, Any}
 # put steps, workdone, stepelased, exception_info into channel
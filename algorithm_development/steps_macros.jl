"""For use in a macro accepting arguments like
@m a,b,c = f(d,e)
returns a Vector{Symbol} containing [:a,:b,:c]"""
function lhs(expr::Expr)
    @assert expr.head == :(=)
    a = expr.args[1]
    isa(a,Symbol) ? Symbol[a] : convert(Vector{Symbol},a.args)
end
"""For use in a macro accepting arguments like
@m a,b,c = f(d,e)
returns a tuple (:f, [:d,:e])"""
function rhs(expr::Expr)
    @assert expr.head == :(=) 
    @assert expr.args[2].head == :call
    expr.args[2].args[1], convert(Vector{Symbol},expr.args[2].args[2:end])
end
"""Accepts something like 
@perpulse `filt_value = do_filter(filter, pulse_record)`
Converts it to `PerPulseStep(do_filter, [:filter, :pulse_record], [:filt_value])`."""
macro perpulse(exp)
f,args = rhs(exp)
outs = lhs(exp)
	quote 
		PerPulseStep($f, $args, $outs)
	end
end

using Base.Test
@test PerPulseStep(do_filter,[:filter,:pulse_record],[:a,:b,:c,:d]) == @perpulse filt_value = sum(filter, pulse_record)
@perpulse a,b,c,d = do_filter(filter, pulse_record)

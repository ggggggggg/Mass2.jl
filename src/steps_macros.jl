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
"""For use in a macro accepting arguments like
f(a)>b
returns a tuple(:f, :a, :b)"""
function conditionsyms(expr::Expr)
    @assert expr.head == :comparison
    @assert expr.args[2] == :>
    expr.args[1].args[1], expr.args[1].args[2], expr.args[3]
end
"""Accepts something like
@perpulse `filt_value = do_filter(filter, pulse_record)`
Converts it to `PerPulseStep(do_filter, [:filter, :pulse_record], [:filt_value])`."""
macro perpulse(expr)
f,args = rhs(expr)
outs = lhs(expr)
	quote
		PerPulseStep($f, $args, $outs)
	end
end
""" Accepts something like
@threshold calibration = calibrate_nofit(filt_value_dc_hist, known_energies) when counted(filt_value_dc_hist)>3000
ThresholdStep(calibrate_nofit, [:filt_value_dc_hist,:known_energies],[:calibration],:filt_value_dc_hist, counted, 3000, true))
"""
macro threshold(expr, when, condition)
    @assert when == :when "2nd argument to @threshold must be :when, this is just for readability"
    f,args = rhs(expr)
    outs = lhs(expr)
    cf, a, b = conditionsyms(condition)
    aa=string(a)
    quote
        ThresholdStep($f, $args, $outs, Symbol($aa), $cf, $b, true)
    end
end
""" Accepts something like
@histogram update_histogram!(filt_value_hist, selection_good, filt_value)
Converts it to HistogramStep(update_histogram!, [:selection_good, :filt_value], [:filt_value_hist])
"""
macro histogram(expr)
    f = expr.args[1]
    args = expr.args[2:end]
    quote
        HistogramStep($f, $(args[2:end]), $(args[1:1]))
    end
end

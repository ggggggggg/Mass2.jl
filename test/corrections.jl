using Base.Test

srand(87687634)
@test length(pow2edges(0.1,1:100))==1025
@test isapprox(first(pow2edges(0.1,1:100)),1)
@test isapprox(last(pow2edges(0.1,1:100)),100)

known_dc = LinearDriftCorrect(0,0.0001)
n=10000
indicator = randn(n)*10+known_dc.indicator_median
ideal_corrected = 1000+randn(n)
uncorrected = ideal_corrected./(1+indicator*known_dc.param)
found_dc = drift_correct(indicator, uncorrected, 0.5)
@test abs(found_dc.indicator_median-known_dc.indicator_median)<1.0
@test abs(1-found_dc.param/known_dc.param)<0.05

corrected = applycorrection(found_dc, indicator, uncorrected)
known_dc = LinearDriftCorrect(0,0.0001)
known_corrected = applycorrection(known_dc, indicator, uncorrected)

@test all(known_corrected-ideal_corrected.<1e-8)

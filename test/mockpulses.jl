using Mass2

pg = TupacLikeTwoExponentialPulseGenerator(Int, Normal(1000,0.5))
getcleanpulse(pg,1000) 
gettriggeredpulse!(pg,10)


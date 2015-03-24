using Mass2

pg = TwoExponentialPulseGenerator{Int}(520, 100, 50, 200,10, 13.3, 100000,520,Normal(1000,1),30,0,1000.,0,0.5)
getcleanpulse(pg,1000) 
gettriggeredpulse!(pg,10)


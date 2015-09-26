using ZMQ, Mass2
context=Context()
sub=Socket(context, SUB)
ZMQ.connect(sub, "tcp://localhost:5555")
ZMQ.set_subscribe(sub, "0")
t0 = time()
function plot_hist(coadded_hist,n)
	fig = figure(1)
	fig[:clf]()
	ax = gca()
	ax[:plot](midpoints(coadded_hist), counts(coadded_hist),label="coadded $n ch")
	xlabel("energy (eV)")
	ylabel("counts per $(binsize(coadded_hist)) bin")
	title("using $n channels, $(counted(coadded_hist)) total counts")
	legend()
	show()
end

try
	while true
		msg = ZMQ.recv(sub)
		println("recieved message $(time()-t0) after starting")
		out=seek(convert(IOStream, msg),1)
		coadded_hist,n = Base.deserialize(out)
		plot_hist(coadded_hist,n)
	end
finally
	close(sub)
	close(context)
end
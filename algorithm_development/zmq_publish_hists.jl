# publish hists over ZMQ
function make_coadded_hist(masschannels)
	coadded_hist = Histogram(0:1:20000)
	n = 0
	for (ch,mc) in masschannels
		if haskey(mc,:energy_hist)
			coadded_hist+=mc[:energy_hist]
			n+=1
		end
	end
	coadded_hist,n
end

function publish_hists_zmq(masschannels, zmqexitchannel)
	context=Context()
	pub=Socket(context, PUB)
	ZMQ.bind(pub, "tcp://*:5555")
	try
		while !isready(zmqexitchannel)
			coadded_hist,n = make_coadded_hist(masschannels)
			io = IOBuffer()
			write(io,"0")
			serialize(io, (coadded_hist,n))
			ZMQ.send(pub, Message(io))
			sleep(1)
		end
	catch ex
		isa(ex, InterruptException) && return
		println(ex)
	finally
		ZMQ.close(pub)
		ZMQ.close(context)
	end
end

zmqexitchannel = Channel()
@schedule publish_hists_zmq(masschannels, zmqexitchannel)
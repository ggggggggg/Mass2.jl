# this is a hack to allow Mass2 to depend on unregistered Julia Packages
# it should be removed if/when there is an official way of doing this other
# than maintaining and using my own copy of METADATA
if !isdir(Pkg.dir("ReferenceMicrocalFiles"))
	Pkg.clone("https://github.com/ggggggggg/ReferenceMicrocalFiles.jl")
end
isdir(expanduser("~/.daq")) || mkdir(expanduser("~/.daq"))
isfile(expanduser("~/.daq/latest_ljh_pulse.cur")) || touch(expanduser("~/.daq/latest_ljh_pulse.cur"))

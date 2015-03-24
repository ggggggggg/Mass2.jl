using GraphViz
using Graphs
using PyCall
# this file shuld be removed eventually, but it make GraphViz work
# without it the @pyimport call in calibration.jl breaks GraphViz somehow
g = simple_wheel_graph(6)

dot = to_dot(g)
gviz = GraphViz.Graph(dot)
GraphViz.layout!(gviz,engine="dot")
open("workaround.svg","w") do f 
       GraphViz.writemime(f, MIME"image/svg+xml"(),gviz)
end #do



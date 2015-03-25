using GraphViz
fname = "graph"
dot = open(readall,"$fname.dot")
gviz = Graph(dot)
GraphViz.layout!(gviz,engine="dot")
open("$fname.svg","w") do f 
       GraphViz.writemime(f, MIME"image/svg+xml"(),gviz)
end #do
open("$fname.png","w") do f 
       GraphViz.writemime(f, MIME"image/png"(),gviz)
end #do
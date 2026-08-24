module QuiverToolsGraphVizExt

# Renders a quiver as SVG by handing its DOT description (see `to_dot` in
# src/Quivers.jl) to Graphviz. Loaded automatically when a user runs
# `using GraphViz` alongside QuiverTools; from then on a `Quiver` shows as a
# drawing in notebooks (Pluto, IJulia) and VS Code, and `draw` opens it in the
# system viewer from the REPL.

using QuiverTools: Quiver, to_dot
using GraphViz: GraphViz
import QuiverTools: draw

Base.show(io::IO, mime::MIME"image/svg+xml", Q::Quiver) =
  show(io, mime, GraphViz.Graph(to_dot(Q)))

function _write_svg(Q::Quiver, path::AbstractString=tempname() * ".svg")
  open(path, "w") do io
    show(io, MIME("image/svg+xml"), Q)
  end
  return path
end

function _open_in_default_viewer(path::AbstractString)
  if Sys.isapple()
    run(`open $path`)
  elseif Sys.iswindows()
    run(`cmd /c start "" $path`)
  else
    run(`xdg-open $path`)
  end
  return nothing
end

function draw(Q::Quiver)
  path = _write_svg(Q)
  _open_in_default_viewer(path)
  return path
end

end

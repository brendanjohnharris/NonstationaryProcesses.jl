using Plots
using Colors
import Plots.PlotThemes._themes
const cornflowerblue = colorant"cornflowerblue"; export cornflowerblue
const crimson = colorant"crimson"; export crimson
const cucumber = colorant"#77ab58"; export cucumber
const california = colorant"#EF9901"; export california
const juliapurple = colorant"#9558b2"; export juliapurple
const keppel = colorant"#46AF98"; export keppel
const darkbg = colorant"#282C34"; export darkbg

function torgba(c::RGB, a::Real=1)
    Colors.RGBA(c.r, c.g, c.b, a)
end

fourseas_palette = torgba.([
    cornflowerblue,
    crimson,
    cucumber,
    california,
    juliapurple
], (0.7,))


_themes[:fourseas] = Plots.PlotThemes.PlotTheme(
    foreground_color_text = :black,
    fgguide = :black,
    fglegend = :black,
    legendfontcolor = :black,
    legendtitlefontcolor = :black,
    titlefontcolor = :black,
    linewidth = 2.5,
    palette = fourseas_palette,
    colorgradient = :viridis,
    framestyle = :grid,
    grid = true,
    minorgrid = true,
    minorgridalpha = 1.0,
    foreground_color_minor_grid = :gray91,
    foreground_color_grid = :gray88,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    gridlinewidth = 1.5,
    gridalpha = 1.0,
    titlefontsize = 12,
    tickfontsize = 10,
    legend = nothing,
    legendfontsize = 10,
    legendtitlefontsize = 10,
    fontfamily = "Computer Modern",
    minorticks = 2,
); Plots.showtheme(:fourseas)

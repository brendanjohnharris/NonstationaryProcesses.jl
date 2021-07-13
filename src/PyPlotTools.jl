import .PyPlot

function set_pane_color(color=(0, 0, 0), ax=PyPlot.gca())
    PyPlot.svg(true)
    ax.xaxis.set_pane_color(color)
    ax.yaxis.set_pane_color(color)
    ax.zaxis.set_pane_color(color)
    f = PyPlot.gcf()
end
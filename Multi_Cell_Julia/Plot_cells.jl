using Plots



# create an ellipse in the sky
pts = Plots.partialcircle(0, 2π, 100, 1.0)
x, y = Plots.unzip(pts)

# spotlight
plot!(Shape(x, y), c = :blue, fillalpha=0.5)

mkdir("Plots")
savefig(plt, "Plots/test.png")

function plot_ellipses(cells::Array{Float64}, cparams::Array{Float64}, 
    Adh::Array{Float64}, Adh0::Array{Float64}, t::Int64)

    # background and limits
    plt = plot(
        bg = :white,
        size = (800, 800),
        legend = false,
        xtickfont = font(10, "match"),
        ytickfont = font(10, "match")
    )

    # For each cell create outer and inner ellipses 
    # Along with that plot all the adhesions too 
    for i in 1:size(cells)[1]÷3
        # create an ellipse
        pts = Plots.partialcircle(0, 2π, 100, 1.0)
        x, y = Plots.unzip(pts)
        x = @. cparams[4*i+1]*x 
        y = @. cparams[4*i+2]*y

        # outer ellipse
        x_tr = @. x*cos(cparams[4*i+3]) - y*sin(cparams[4*i+3]) + cells[3*i+1]
        y_tr = @. x*sin(cparams[4*i+3]) + y*cos(cparams[4*i+3]) + cells[3*i+2]
        plot!(Shape(x_tr, y_tr), c = :blue, fillalpha=0.5)

        # inner ellipse
        x_tr = @. x*cos(cparams[4*i+3])/2. - y*sin(cparams[4*i+3])/2. + cells[3*i+1]
        y_tr = @. x*sin(cparams[4*i+3])/2. + y*cos(cparams[4*i+3])/2. + cells[3*i+2]
        plot!(Shape(x_tr, y_tr), c = :red, fillalpha=1)

    end    

end    


using Random
using Distributions 
include("EE_Ovlap_Area.jl")

# Function to check if new cell collides with any other cell present
function noCollision(cells, cparams, cell_index, a, b, x, 
                    y, theta)

    # No need to check for initial cell                    
    if cell_index == 1
        return 1
    end               
    
    x1 = x
    y1 = y 
    a1 = a
    b1 = b 
    t1 = theta 

    for i in 1:cell_index
        x2 = cells[3*i+1]
        y2 = cells[3*i+2]
        a2 = cparams[4*i+1]
        b2 = cparams[4*i+2]
        t2 = cparams[4*i+3]

        #intersection area 
        intersec = ellipse_ellipse_overlap(a1, b1, x1, y1, t1, 
        a2, b2, x2, y2, t2)

        if intersec > 0.0
            return 0
        end
    end

    # cell doesn't overlaps with existing cells
    return 1
end                    

# A two cell preset
function two_cell_preset(a, b, Num)

    # x, y and dtheta D.O.F 
    cells = zeros(Float64, 3*Num)
    # a b theta and phase
    cparams = zeros(Float64, 4*Num)
    # Overlap array
    Ovlaps = zeros(Int64, Num, Num)

    # Cell 1
    ind = 1
    x = 20.
    y = 20.
    theta = 0.0
    cells[3*(ind-1)+1:(3*(ind-1)+3)] = Float64[x, y, 0.0]
    cparams[4*(ind-1)+1:(4*(ind-1)+4)] = Float64[a, b, theta, 1]

    # Cell 2
    ind = 2
    x = 18.
    y = 13.
    theta = pi/4
    cells[3*(ind-1)+1:(3*(ind-1)+3)] = Float64[x, y, 0.0]
    cparams[4*(ind-1)+1:(4*(ind-1)+4)] = Float64[a, b, theta, 1]

    return cells, cparams, Ovlaps
end

# Linear preset 
function linear_preset(a, b, Num)

    # x, y and dtheta D.O.F 
    cells = zeros(Float64, 3*Num)
    # a b theta and phase
    cparams = zeros(Float64, 4*Num)
    # Overlap array
    Ovlaps = zeros(Int64, Num, Num)
    
    # Linear in X - direction
    for ind in 1:Num
        cells[3*(ind-1)+1:(3*(ind-1)+3)] = Float64[20+2*(a-1)*ind, 10, 0.0]
        cparams[4*(ind-1)+1:(4*(ind-1)+4)] = Float64[a, b, 0.0, 1] 
    end    

    return cells, cparams, Ovlaps
end    

# Random configuration of non-overlapping cells 
function random_cells(L, a, b, Num)

    # x, y and dtheta D.O.F 
    cells = zeros(Float64, 3*Num)
    # a b theta and phase
    cparams = zeros(Float64, 4*Num)
    # Overlap array
    Ovlaps = zeros(Int64, Num, Num)

    # cell index 
    cell_ind = 0

    while (cell_ind <= Num)
        x = rand(Uniform(0.0, L))
        y = rand(Uniform(0.0, L))
        theta = rand(Uniform(0.0, 2.0*pi))

        # check if the cell is overlapping with others 
        if (noCollision(cells, cparams, cell_ind, a, b, x,
                        y, theta) == 1)
            cells[3*(ind-1)+1:(3*(ind-1)+3)] = Float64[x, y, 0.0]
            cparams[4*(ind-1)+1:(4*(ind-1)+4)] = Float64[a, b, theta, 1]  
            cell_ind = cell_ind + 1  
        end                
    end    

    return cells, cparams, Ovlaps
end    

# Get the overlap indices of all the cells
function find_overlap!(cells, cparams, Ovlaps)

    # iterate through all the cells 
    for i in 1:(size(cells)[1]÷3)
        x1 = cells[3*(i-1)+1]
        y1 = cells[3*(i-1)+2]
        a1 = cparams[4*(i-1)+1]
        b1 = cparams[4*(i-1)+2]
        t1 = cparams[4*(i-1)+3]
        for j in (i+1):(size(cells)[1]÷3)
            x2 = cells[3*(j-1)+1]
            y2 = cells[3*(j-1)+2]
            a2 = cparams[4*(j-1)+1]
            b2 = cparams[4*(j-1)+2]
            t2 = cparams[4*(j-1)+3] 

            # Find the intersection area
            if ellipse_ellipse_overlap(a1, b1, x1, y1, t1, 
                a2, b2, x2, y2, t2) > 0.0
                Ovlaps[i, j] = 1
                Ovlaps[j, i] = 1
            else
                Ovlaps[i, j] = 0
                Ovlaps[j, i] = 0     
            end    
        end
    end    

    return Ovlaps
end    
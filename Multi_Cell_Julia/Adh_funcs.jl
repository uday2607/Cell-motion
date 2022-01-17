include("Cell_funcs.jl")
include("Energy_funcs.jl")
using NLsolve

# Spawn random adhesions 
function random_adhesions(cells, cparams, Nadh, k_plus, dt)

    # Adhesions arrays 
    Adh = -10^8 .* ones((size(cells)[1]÷3, Nadh, 2))
    Adh0 = -10^8 .* ones((size(cells)[1]÷3, Nadh, 2))
    cAdh0 = -10^8 .* ones((size(cells)[1]÷3, Nadh, 2))
    cp0 = -10^8 .* ones((size(cells)[1]÷3, Nadh, 2))

    # For every cell create adhesions
    for ind in 1:(size(cells)[1]÷3)
        
        # Random Numbers 
        phi = rand(Normal(0, pi), Nadh)
        rho = rand(Uniform(0.0, 1.0), Nadh)

        # Check if the adhesion can form 
        for i in 1:Nadh
            
            # In the system of reference
            x = cparams[4*(ind-1)+1]*sqrt(rho[i])*cos(phi[i])
            y = cparams[4*(ind-1)+2]*sqrt(rho[i])*sin(phi[i])

            # Transform to ellipse frame of reference
            xp = (x*cos(cparams[4*(ind-1)+3]) - 
                    y*sin(cparams[4*(ind-1)+3]))
            yp = (x*sin(cparams[4*(ind-1)+3]) + 
                    y*cos(cparams[4*(ind-1)+3]))

            # Probability of adhesion formation
            if rand() < k_plus*dt
                
                # from ground frame of reference
                Adh0[ind, i, 1] = xp + cells[3*(ind-1)+1]
                Adh0[ind, i, 2] = yp + cells[3*(ind-1)+2]

                # from the center of the ellipse 
                Adh[ind, i, 1] = xp
                Adh[ind, i, 2] = yp 
                cAdh0[ind, i, 1] = xp
                cAdh0[ind, i, 2] = yp 

                # cell coordinates
                cp0[ind, i, 1] = cparams[4*(ind-1)+1]
                cp0[ind, i, 2] = cparams[4*(ind-1)+3]
            end    
        end
    end    

    return Adh0, Adh, cAdh0, cp0
end    

# create random adhesions for one cell 
function one_cell_random_adh!(cells, ind, cparams, Adh, Adh0, cAdh0,
                            cp0, Nadh, k_plus, dt, flag)

    # Random Numbers 
    if flag == 1
        phi = rand(Normal(0, pi/4), Nadh)
        rho = rand(Uniform(0.5, 1.0), Nadh)
    else 
        phi = rand(Normal(0, pi), Nadh)
        rho = rand(Uniform(0.0, 1.0), Nadh)
    end

    # Check if the adhesion can form 
    for i in 1:Nadh
        
        # In the system of reference
        x = cparams[4*(ind-1)+1]*sqrt(rho[i])*cos(phi[i])
        y = cparams[4*(ind-1)+2]*sqrt(rho[i])*sin(phi[i])

        # Transform to ellipse frame of reference
        xp = (x*cos(cparams[4*(ind-1)+3]) - 
                y*sin(cparams[4*(ind-1)+3]))
        yp = (x*sin(cparams[4*(ind-1)+3]) + 
                y*cos(cparams[4*(ind-1)+3]))

        # Probability of adhesion formation
        if rand() < k_plus*dt
            
            # from ground frame of reference
            Adh0[ind, i, 1] = xp + cells[3*(ind-1)+1]
            Adh0[ind, i, 2] = yp + cells[3*(ind-1)+2]

            # from the center of the ellipse 
            Adh[ind, i, 1] = xp
            Adh[ind, i, 2] = yp 
            cAdh0[ind, i, 1] = xp
            cAdh0[ind, i, 2] = yp

            # cell coordinates
            cp0[ind, i, 1] = cparams[4*(ind-1)+1]
            cp0[ind, i, 2] = cparams[4*(ind-1)+3]
        end    
    end
end                            

# rotate and translate the adhesions 
function rotation_and_shift!(cells, cparams, Adh)

    # translate and shift all the adhesions of the cells 
    for ind in 1:(size(cells)[1]÷3)
        
        # iterate through all the adhesions 
        for i in LinearIndices(Adh[ind, :, 1])
            
            # If adhesion exists 
            if Adh[ind, i, 1] != -10^8 && Adh[ind, i, 2] != -10^8
                
                # Adhesion coordinates
                x = Adh[ind, i, 1]
                y = Adh[ind, i, 2]
                dtheta = cells[3*(ind-1)+3]   

                # rotate the adhesions
                Adh[ind, i, 1] = cos(dtheta)*x - sin(dtheta)*y 
                Adh[ind, i, 2] = sin(dtheta)*x + cos(dtheta)*y
            end    
        end    

        # rotate the angle of polarity vector 
        cparams[4*(ind-1)+3] += cells[3*(ind-1)+3]
        cells[3*(ind-1)+3] = 0
    end
end    

# Contraction of cells -> Contraction of Adhesion sites 
function contraction!(cells, num, cparams, Ovlap_indices, Adh, Adh0,
                    lamda, tau, dt, T_S, k_oo, k_io, k_ii,
                    k_s, a_min)

    # compute the tension of a cell
    T = compute_tension(cells, num, cparams, Ovlap_indices, Adh, 
                        Adh0, k_oo, k_io, k_ii, k_s)
                        
    # Different cases depending on the Tension 
    # Tension is +ve and less than critical tension 
    # So there is a stalling tension which lowers contraction
    # speed
    if T < T_S && T > 0.0
        c = lamda*dt*(1-T/T_S)/tau

        # Update cell length and phase 
        cparams[4*(num-1)+1] -= c*cparams[4*(num-1)+1]
        cparams[4*(num-1)+4] += 1
        phi = cparams[4*(num-1)+3]

        # Iterate through all the adhesions 
        for i in LinearIndices(Adh[num, :, 1])
            # If adhesion exists 
            if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
                x = Adh[num, i, 1]
                y = Adh[num, i, 2]

                Adh[num, i, 1] = (x - c*(cos(phi)*cos(phi)*x+
                                sin(phi)*cos(phi)*y))
                Adh[num, i, 2] = (y - c*(sin(phi)*cos(phi)*x+
                                sin(phi)*sin(phi)*y))
            end    
        end
    # Tension is -ve. Cell contraction speed is normal    
    elseif T < 0.0    
        c = lamda*dt/tau

        # Update cell length and phase 
        cparams[4*(num-1)+1] -= c*cparams[4*(num-1)+1]
        cparams[4*(num-1)+4] += 1
        phi = cparams[4*(num-1)+3]

        # Iterate through all the adhesions 
        for i in LinearIndices(Adh[num, :, 1])
            # If adhesion exists 
            if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
                x = Adh[num, i, 1]
                y = Adh[num, i, 2]

                Adh[num, i, 1] = (x - c*(cos(phi)*cos(phi)*x+
                                sin(phi)*cos(phi)*y))
                Adh[num, i, 2] = (y - c*(sin(phi)*cos(phi)*x+
                                sin(phi)*sin(phi)*y))
            end    
        end
    # Tension is more than Stalling Tension
    # Cell is stalled
    # No translation of adhesions 
    else
        # Update phase 
        cparams[4*(num-1)+4] += 1    
    end     
    
    # Change the phase if the contraction phase reaches it's
    # end (when a < a_min)
    if cparams[4*(num-1)+4] > 0 && cparams[4*(num-1)+1] < a_min
        cparams[4*(num-1)+4] = -1 # -ve times -> Protrusion
    end
end                    

# Function to find the center of the cell which makes the 
# ellipse pass through the rear-most adhesion site
function solve_center!(F, vals, a, b, phi, h, k, xc, yc)

    # variables
    x = vals[1]
    y = vals[2]

    # angle in the range [-pi, pi]
    phi -= 2.0*pi*floor((phi + pi)*(1.0/(2.0*pi)))

    # Evaluate the function 
    F[1] = (((h+xc-x)*cos(phi)+(k+yc-y)*sin(phi))^2.0/a^2.0 + 
            ((h+xc-x)*sin(phi)-(k+yc-y)*cos(phi))^2.0/b^2.0 - 1)
    F[2] = atan(y-yc, x-xc) - phi
end    

# Function to find the shortest distance from a line 
function shortest_distance(points, a, b, c)

    x1 = points[1]
    y1 = points[2]

    d = abs(a*x1+b*y1+c)/sqrt(a*a+b*b)

    return d
end    

# Protrusion of cells -> Protrusion of Adhesion sites
function protrusion!(cells, num, Adh, Adh0, cparams, Ovlap_indices, T_S,
                    lamda, k_s, k_oo, k_io, k_ii, dt, tau, a)

    # compute the tension of a cell
    T = -1.0*compute_tension(cells, num, cparams, Ovlap_indices, Adh, 
                        Adh0, k_oo, k_io, k_ii, k_s)  
    # (-ve sign for Tension to account for opposite direction of gradient)     
    
    # Adhesion speed 
    vf = 1.0 

    # Different cases depending on the Tension 
    # Tension is +ve and less than critical tension 
    # So there is a stalling tension which lowers contraction
    # speed
    if (T < T_S && T > 0.0)

        c = lamda*dt*(1 - T/T_S)/(tau÷5)

        # Update 'a' and phase 
        cparams[4*(num-1)+1] += c*cparams[4*(num-1)+1]
        cparams[4*(num-1)+4] -= 1

        # Find the perpenicular line to semi-major axis 
        a1 = -1/tan(cparams[4*(num-1)+3] + 10^-8)
        b1 = -1 
        c1 = (-a1*(cparams[4*(num-1)+1]*cos(cparams[4*(num-1)+3]) + 
                 cparams[4*(num-1)+1]*sin(cparams[4*(num-1)+3])))

        # Find the rear most adhesion   
        ad = -1
        max = 0.0
        for i in LinearIndices(Adh[num, :, 1])
            if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
                dis = shortest_distance(Adh[num, i, :], a1, b1, c1)
                if dis > max
                   max = dis 
                   ad = i 
                end  
            end      
        end    
        
        # If there are adhesions, find the new center of the cell 
        if ad != -1
            theta = cparams[4*(num-1)+3]

            # Solving center
            xsol = nlsolve((F,x) -> solve_center!(F, x, cparams[4*(num-1)+1],
            cparams[4*(num-1)+2], theta, Adh[num, ad, 1], 
            Adh[num, ad, 2], cells[3*(num-1)+1], cells[3*(num-1)+2]), 
            cells[3*(num-1)+1:3*(num-1)+2] .+ Float64[cos(theta), sin(theta)])
            x_c, y_c = xsol.zero

            dis = sqrt((cells[3*(num-1)+1] - x_c)^2.0 +
                        (cells[3*(num-1)+2] - y_c)^2.0)            
            xn = cells[3*(num-1)+1] + ((vf*5*dis)/tau)*cos(theta)
            yn = cells[3*(num-1)+2] + ((vf*5*dis)/tau)*sin(theta)     
            
            # Iterate through all the adhesions 
            for i in LinearIndices(Adh[num, :, 1])     

                # If adhesion exists 
                if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
                    x = Adh[num, i, 1]
                    y = Adh[num, i, 2]   

                    Adh[num, i, 1] = (x + c*(cos(theta)*cos(theta)*x+
                                    sin(theta)*cos(theta)*y) + cells[3*(num-1)+1] - xn)
                    Adh[num, i, 2] = (y + c*(sin(theta)*cos(theta)*x+
                                    sin(theta)*sin(theta)*y) + cells[3*(num-1)+2] - yn)
                end    
            end

            # Shift the center of the cells 
            cells[3*(num-1)+1] = xn 
            cells[3*(num-1)+2] = yn
        end
    # Tension is -ve. Cell Protrusion speed is normal     
    elseif (T < 0.0)
        c = lamda*dt/(tau÷5)

        # Update 'a' and phase 
        cparams[4*(num-1)+1] += c*cparams[4*(num-1)+1]
        cparams[4*(num-1)+4] -= 1

        # Find the perpenicular line to semi-major axis 
        a1 = -1/tan(cparams[4*(num-1)+3] + 10^-8)
        b1 = -1 
        c1 = (-a1*(cparams[4*(num-1)+1]*cos(cparams[4*(num-1)+3]) + 
                 cparams[4*(num-1)+1]*sin(cparams[4*(num-1)+3])))

        # Find the rear most adhesion   
        ad = -1
        max = 0.0
        for i in LinearIndices(Adh[num, :, 1])
            if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
                dis = shortest_distance(Adh[num, i, :], a1, b1, c1)
                if dis > max
                   max = dis 
                   ad = i 
                end  
            end      
        end    
        
        # If there are adhesions, find the new center of the cell 
        if ad != -1
            theta = cparams[4*(num-1)+3]

            # Solving center 
            xsol = nlsolve((F,x) -> solve_center!(F, x, cparams[4*(num-1)+1],
            cparams[4*(num-1)+2], theta, Adh[num, ad, 1], 
            Adh[num, ad, 2], cells[3*(num-1)+1], cells[3*(num-1)+2]), 
            cells[3*(num-1)+1:3*(num-1)+2] .+ Float64[cos(theta), sin(theta)])
            x_c, y_c = xsol.zero

            dis = sqrt((cells[3*(num-1)+1] - x_c)^2.0 +
                        (cells[3*(num-1)+2] - y_c)^2.0)
            xn = cells[3*(num-1)+1] + ((vf*5*dis)/tau)*cos(theta)
            yn = cells[3*(num-1)+2] + ((vf*5*dis)/tau)*sin(theta)     
            
            # Iterate through all the adhesions 
            for i in LinearIndices(Adh[num, :, 1])     

                # If adhesion exists 
                if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
                    x = Adh[num, i, 1]
                    y = Adh[num, i, 2]

                    Adh[num, i, 1] = (x + c*(cos(theta)*cos(theta)*x+
                                    sin(theta)*cos(theta)*y) + cells[3*(num-1)+1] - xn)
                    Adh[num, i, 2] = (y + c*(sin(theta)*cos(theta)*x+
                                    sin(theta)*sin(theta)*y) + cells[3*(num-1)+2] - yn)
                end    
            end

            # Shift the center of the cells 
            cells[3*(num-1)+1] = xn 
            cells[3*(num-1)+2] = yn
        end
    # Tension is more than Stalling Tension
    # Cell is stalled
    # No translation of adhesions     
    else    
        # update phase 
        cparams[4*(num-1)+4] -= 1
    end    

    # Change the phase if the protrusion phase reaches it's
    # end (when a > a_max)
    if cparams[4*(num-1)+4] < 0 && cparams[4*(num-1)+1] > a
        cparams[4*(num-1)+4] = 0 # non -ve times -> Contraction
    end
end                    

# After the adhesions formations, check if any adheisons will 
# mature 
function mature!(cells, num, Adh, Adh0, cAdh0, cp0, k_m, fTh, dt)

    # iterate through all the adhesions 
    for i in LinearIndices(Adh[num, :, 1])
        # If adhesion exists 
        if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
            dis = sqrt((Adh[num, i, 1]+cells[3*(num-1)+1] - Adh0[num, i, 1])^2.0 + 
                        (Adh[num, i, 2]+cells[3*(num-1)+2] - Adh0[num, i, 2])^2.0)
            off_rate = k_m*exp(-dis/fTh)

            # mature adhesion site 
            if (rand() < 1-off_rate*dt)
                Adh[num, i, :] = Float64[-10^8, -10^8]
                Adh0[num, i, :] = Float64[-10^8, -10^8]
                cAdh0[num, i, :] = Float64[-10^8, -10^8]
                cp0[num, i, :] = Float64[-10^8, -10^8]
            end    
        end    
    end
end

# Detach the bonds 
function detach!(cells, num, cAdh0, cp0, Adh, Adh0, k_b, k_f, alpha, 
                dt)

    # iterate through all the adhesions 
    for i in LinearIndices(Adh[num, :, 1])
        # If adhesion exists 
        if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
            x0 = (cAdh0[num, i, 1]*cos(cp0[num, i, 2]) + 
                  cAdh0[num, i, 2]*sin(cp0[num, i, 2]))

            #detachment rate 
            k_x = k_b - (k_b - k_f)*(x0+cp0[num, i, 1])/(2.0*cp0[num, i, 1])   
            
            # off rate 
            dis = sqrt((Adh[num, i, 1]+cells[3*(num-1)+1] - Adh0[num, i, 1])^2.0 + 
                        (Adh[num, i, 2]+cells[3*(num-1)+2] - Adh0[num, i, 2])^2.0)
            # The offrate increases as cell length decreases
            off_rate = k_x*exp(alpha*dis/(2.0*cp0[num, i, 1]))

            # detach adhesion site 
            if (rand() < off_rate*dt)
                Adh[num, i, :] = Float64[-10^8, -10^8]
                Adh0[num, i, :] = Float64[-10^8, -10^8]
                cAdh0[num, i, :] = Float64[-10^8, -10^8]
                cp0[num, i, :] = Float64[-10^8, -10^8]
            end    
        end    
    end
end                
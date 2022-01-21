include("ell_ovlap_area.jl")
include("Cell_funcs.jl")

# Calculate the adhesion energy of each cell 
function one_adh_energy(cells, num, Adh, Adh0, k_s)
    
    # Temporary variables
    adh_energy = 0.0
    f_x = 0.0
    f_y = 0.0
    x = 0.0
    y = 0.0

    # Change in angle of the polarity of the cell 
    dtheta = cells[3*num]
    cos_t = cos(dtheta)
    sin_t = sin(dtheta)

    # iterate through all the adhesion sites info
    for i in LinearIndices(Adh[num, :, 1])
        if Adh[num, i, 1] != -10^8 && Adh[num, i, 2] != -10^8
            x = Adh[num, i, 1]
            y = Adh[num, i, 2]
            f_x = -(cos_t*x-sin_t*y + cells[3*(num-1)+1]
                    -Adh0[num, i, 1])
            f_y = -(sin_t*x+cos_t*y + cells[3*(num-1)+2]
                    -Adh0[num, i, 2])   
            
            # Adh energy of the bond 
            adh_energy += f_x*f_x + f_y*f_y    
        end    
    end    

    return 0.5*k_s*adh_energy
end

# Calculate the ovarlap energy of two cells 
function pairwise_ovlap_energy_i(cells, cparams, ind_i, ind_j, 
                              k_oo, k_io, k_ii)

    # temporary variables                            
    oval_energy = 0.0
    area_oo = 0.0
    area_io = 0.0 
    area_ii = 0.0

    # cell coordinates    
    x1 = cells[3*(ind_i-1)+1]
    y1 = cells[3*(ind_i-1)+2]
    a1 = cparams[4*(ind_i-1)+1]
    b1 = cparams[4*(ind_i-1)+2]
    t1 = cparams[4*(ind_i-1)+3]

    x2 = cells[3*(ind_j-1)+1]
    y2 = cells[3*(ind_j-1)+2]
    a2 = cparams[4*(ind_j-1)+1]
    b2 = cparams[4*(ind_j-1)+2]
    t2 = cparams[4*(ind_j-1)+3]


    # angles of cell polarities (all angles in range [-pi, pi])
    phi_i = t1 + cells[3*(ind_i-1)+3]
    phi_i -= 2.0*pi*floor((phi_i + pi)*(1.0/(2.0*pi)))
    phi_j = t2 + cells[3*(ind_j-1)+3]
    phi_j -= 2.0*pi*floor((phi_j + pi)*(1.0/(2.0*pi)))
    phi_c = atan(y2 - y1, x2 - x1)
    ovlap_const = 0.5*((1.0+0.5*cos(2.0*(phi_c - phi_i)))*
                   (1.0+0.5*cos(2.0*(phi_c - phi_j))))        

    # find the overlap areas of cells
    # (partial function) 
    ee_area(a1, b1, a2, b2) = ellipse_ellipse_overlap(a1, b1, x1, y1, t1, 
                                a2, b2, x2, y2, t2)       
    area_oo = ee_area(a1, b1, a2, b2)

    # outer ellipses overlap                 
    if area_oo > 0.0
        oval_energy -= k_oo*ovlap_const*area_oo

        # in-out ellipses overlap 
        area_io = ee_area(a1/2., b1/2., a2, b2)
        if area_io > 0.0
            oval_energy += k_io*ovlap_const*area_io
        end
        
        # out-in ellipses overlap 
        area_oi = ee_area(a1, b1, a2/2., b2/2.)
        if area_oi > 0.0
            oval_energy += k_io*ovlap_const*area_oi

            area_ii = ee_area(a1/2., b1/2., a2/2., b2/2.) 
            if area_ii > 0.0
                # inner ellipses overlap 
                oval_energy += k_ii*ovlap_const*area_ii 
            end    
        end
    end                

    return oval_energy
end

# Calculate the ovarlap energy of two cells 
function pairwise_ovlap_energy_e(ell_i_in, ell_i_out, ell_j_in, ell_j_out, 
            ovlap_const, k_oo, k_io, k_ii)

    # temporary variables                            
    oval_energy = 0.0
    area_oo = 0.0
    area_io = 0.0 
    area_ii = 0.0       

    # find the overlap areas of cells
    # (partial function)        
    area_oo = ellipse_ellipse_overlap_e(ell_i_out, ell_j_out)

    # outer ellipses overlap                 
    if area_oo > 0.0
        oval_energy -= k_oo*ovlap_const*area_oo

        # in-out ellipses overlap 
        area_io = ellipse_ellipse_overlap_e(ell_i_in, ell_j_out)
        if area_io > 0.0
            oval_energy += k_io*ovlap_const*area_io
        end

        # out-in ellipses overlap 
        area_oi = ellipse_ellipse_overlap_e(ell_i_out, ell_j_in)
        if area_oi > 0.0
            oval_energy += k_io*ovlap_const*area_oi

            area_ii = ellipse_ellipse_overlap_e(ell_i_in, ell_j_in) 
            if area_ii > 0.0
                # inner ellipses overlap 
                oval_energy += k_ii*ovlap_const*area_ii 
            end    
        end
    end                

    return oval_energy
end

# calculate the total overlap energy of the cell 
function one_oval_energy(cells, cparams, ind_i,
                    Ovlap_indices, k_oo, k_io, k_ii)
    
    # temporary variables 
    oval_energy = 0.0 
    
    # check if the cell overlaps with any other cells 
    pairwise_oval_energy(i, j) = pairwise_ovlap_energy_i(cells, cparams, 
                                i, j, k_oo, k_io, k_ii)
    for ind_j in LinearIndices(Ovlap_indices[ind_i, :])
        
        # check for individual overlap 
        if Ovlap_indices[ind_i, ind_j] == 1
            oval_energy += pairwise_oval_energy(ind_i, ind_j)
        end
    end    

    return oval_energy
end                    

# Virtual contraction of cells
function virtual_contraction!(cparam, Adh, h)
    
    # Temporary variables 
    phi = cparam[3]
    x = 0.0
    y = 0.0

    # Update the contractions
    c = h/cparam[1]

    # Update the "length(a)"
    cparam[1] -= h

    for i in LinearIndices(Adh[:, 1])
        
        # If the adhesion exists, translate 
        if Adh[i, 1] != -10^8 && Adh[i, 2] != -10^8
            # coordinates of adhesion
            x = Adh[i, 1]
            y = Adh[i, 2]

            Adh[i, 1] = (x - c*(cos(phi)*cos(phi)*x+
                                sin(phi)*cos(phi)*y))
            Adh[i, 2] = (y - c*(sin(phi)*cos(phi)*x+
                                sin(phi)*sin(phi)*y))                  
        end    
    end    

    return cparam, Adh
end 

# Virtual extension of cells
function virtual_extension!(cparam, Adh, h)
    
    # Temporary variables 
    phi = cparam[3]
    x = 0.0
    y = 0.0

    # Update the extension
    c = h/cparam[1]

    # Update the "length(a)"
    cparam[1] += h

    for i in LinearIndices(Adh[:, 1])
        
        # If the adhesion exists, translate 
        if Adh[i, 1] != -10^8 && Adh[i, 2] != -10^8
            # coordinates of adhesion
            x = Adh[i, 1]
            y = Adh[i, 2]

            Adh[i, 1] = (x + c*(cos(phi)*cos(phi)*x+
                                sin(phi)*cos(phi)*y))
            Adh[i, 2] = (y + c*(sin(phi)*cos(phi)*x+
                                sin(phi)*sin(phi)*y))                  
        end    
    end    

    return cparam, Adh
end

# Compute the tension on the cell 
function compute_tension(cells, num, cparams, Ovlap_indices, Adh, 
                        Adh0, k_oo, k_io, k_ii, k_s)

    # Make copies 
    cps = copy(cparams)
    adh = copy(Adh)

    # temporary variables 
    step = 10^-4
    c_ind = 4*(num-1)
    oval_energy(cparams) = one_oval_energy(cells, cparams, num, 
                            Ovlap_indices, k_oo, k_io, k_ii)
    adh_energy(Adh) = one_adh_energy(cells, num, Adh, Adh0, k_s)                        

    # Use a Mid-point formula diff formula 
    # Extension - a + h 
    cps[c_ind+1:c_ind+4], adh[num, :, :] = virtual_extension!(cparams[c_ind+1:c_ind+4],
                                        Adh[num, :, :], step)
    e2 = oval_energy(cps)
    e2 += adh_energy(adh) 
    
    # Contraction - a - h 
    cps[c_ind+1:c_ind+4], adh[num, :, :] = virtual_contraction!(cparams[c_ind+1:c_ind+4],
                                        Adh[num, :, :], step)
    e1 = oval_energy(cps)
    e1 += adh_energy(adh) 

    # Tension -> T = -dU/da    
    return -(e2 - e1)/(2*step)
end                        

# Total Energy of the system 
function total_energy(cells, cparams, Ovlap_indices, Adh, Adh0, 
                        k_s, k_oo, k_io, k_ii)

    # temporary variables
    total_energy = 0.0
    E = 0.0
    adh_energy(num) = one_adh_energy(cells, num, Adh, Adh0, k_s)
    oval_energy(ell_i_in, ell_i_out, ell_j_in, ell_j_out, 
    ovlap_const) = pairwise_ovlap_energy_e(ell_i_in, ell_i_out, ell_j_in, ell_j_out, 
                                ovlap_const, k_oo, k_io, k_ii)

    # Find the overlap indices 
    find_overlap!(cells, cparams, Ovlap_indices)
    z_size = size(Ovlap_indices)
    E_ovlaps = zeros(Real, z_size)

    # Iterate thorugh all the cells 
    Threads.@threads for ind_i in 1:(size(cells)[1]÷3)
    #for ind_i in 1:(size(cells)[1]÷3)  
    
        # cell coordinates    
        x1 = cells[3*(ind_i-1)+1]
        y1 = cells[3*(ind_i-1)+2]
        a1 = cparams[4*(ind_i-1)+1]
        b1 = cparams[4*(ind_i-1)+2]
        t1 = cparams[4*(ind_i-1)+3]

        # angles of cell polarities (all angles in range [-pi, pi])
        phi_i = t1 + cells[3*(ind_i-1)+3]
        phi_i -= 2.0*pi*floor((phi_i + pi)*(1.0/(2.0*pi))) 

        # ellipses 
        ell_i_out = ellipse(x1, y1, a1, b1, t1, 45)
        ell_i_in = ellipse(x1, y1, a1/2, b1/2, t1, 45)
    
        # adhesion energy 
        total_energy += adh_energy(ind_i)

        # overlap energy 
        for ind_j in LinearIndices(Ovlap_indices[ind_i, :])

            x2 = cells[3*(ind_j-1)+1]
            y2 = cells[3*(ind_j-1)+2]
            a2 = cparams[4*(ind_j-1)+1]
            b2 = cparams[4*(ind_j-1)+2]
            t2 = cparams[4*(ind_j-1)+3]

            # angles of cell polarities (all angles in range [-pi, pi])
            phi_j = t2 + cells[3*(ind_j-1)+3]
            phi_j -= 2.0*pi*floor((phi_j + pi)*(1.0/(2.0*pi)))
            phi_c = atan(y2 - y1, x2 - x1)
            ovlap_const = 0.5*((1.0+0.5*cos(2.0*(phi_c - phi_i)))*
                            (1.0+0.5*cos(2.0*(phi_c - phi_j))))

            # ellipses 
            ell_j_out = ellipse(x2, y2, a2, b2, t2, 45)
            ell_j_in = ellipse(x2, y2, a2/2, b2/2, t2, 45)

            if Ovlap_indices[ind_i, ind_j] == 1
                
                # If Pairwise energy isn't calculated 
                # yet, do it 
                if E_ovlaps[ind_i, ind_j] == 0.0
                    E = oval_energy(ell_i_in, ell_i_out, ell_j_in, ell_j_out, 
                    ovlap_const)
                    E_ovlaps[ind_i, ind_j] = E
                    E_ovlaps[ind_j, ind_i] = E
                end
            end    
        end    
    end    

    total_energy += sum(E_ovlaps)
    return total_energy                    
end
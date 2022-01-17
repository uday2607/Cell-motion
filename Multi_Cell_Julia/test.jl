include("Adh_funcs.jl")
using Optim
using NPZ

function main()

    a = 6.0
    b = 3.0
    L = 40.0
    Num = 5
    Nadh = 64
    k_plus = 0.5
    dt = 1
    lamda = 0.5
    tau = 30
    T_S = 100.0
    k_s = 0.5
    k_m = 1.0
    k_oo = 5.0
    k_io = k_oo*100.0
    k_ii = k_io*100.0
    fThreshold = 0.5
    k_b = 0.05
    k_f = 0.01
    alpha = 50.0
    TIME = 100
    a_min = a
    for i in 1:tau
        a_min -= a_min*lamda*dt/tau
    end    

    # Spwan cells and adhesions
    cells, cparams, Ovlap_indices = linear_preset(a, b, Num)
    Adh0, Adh, cAdh0, cp0 = random_adhesions(cells, cparams, Nadh, k_plus, dt)

    CELLS = zeros(2*TIME, 3*Num)
    CPARAMS = zeros(2*TIME, 4*Num)
    ADH = zeros(2*TIME, Num, Nadh, 2)
    ADH0 = zeros(2*TIME, Num, Nadh, 2)

    for t in 1:(TIME) 

        println("Time = ", t)

        # Find the overlap indices
        Ovlap_indices = find_overlap!(cells, cparams, Ovlap_indices)

        # Update cell coordinates 
        for num in 1:Num
            if cparams[4*num] > 0
                # contraction 
                println("Cell Contraction #", num)

                contraction!(cells, num, cparams, Ovlap_indices, Adh, Adh0, 
                lamda, tau, dt, T_S, k_oo, k_io, k_ii, k_s, a_min)                                

                # bonds mature or detach 
                if cparams[4*num] == 2
                    mature!(cells, num, Adh, Adh0, cAdh0, cp0, 
                            k_m, fThreshold/k_s, dt)
                elseif (cparams[4*num] != 1 && cparams[4*num] != 0 && cparams[4*num] != 2)
                    detach!(cells, num, cAdh0, cp0, Adh, Adh0, 
                            k_b, k_f, alpha, dt)                            
                end     
            
            elseif cparams[4*num] < 0    
                # protrusion
                println("Cell Protrusion #", num)

                protrusion!(cells, num, Adh, Adh0, cparams, Ovlap_indices, T_S, lamda, 
                            k_s, k_oo, k_io, k_ii, dt, tau, a)

                # Old bonds detach -> New bonds form  
                if cparams[4*num] != 0
                    detach!(cells, num, cAdh0, cp0, Adh, Adh0, 
                            k_b, k_f, alpha, dt)    
                end                 

                # New adhesions at a smaller rate 
                #one_cell_random_adh!(cells, num, cparams, Adh, Adh0, cAdh0, cp0,
                #                    Nadh, k_plus/10, dt, 1)
            
            else 
                # New bonds 
                println("New adhesions #", num)

                one_cell_random_adh!(cells, num, cparams, Adh, Adh0, cAdh0, cp0,
                                    Nadh, k_plus, dt, 0)
                cparams[4*num] = 1
            end     
        end

        # store the data
        CELLS[(2*t-1), :] = cells 
        CPARAMS[(2*t-1), :] = cparams 
        ADH[(2*t-1), :, :, :] = Adh[:, :, :]
        ADH0[(2*t-1), :, :, :] = Adh0[:, :, :]

        # Energy minimization
        energy(cells) = total_energy(cells, cparams, Ovlap_indices, Adh, Adh0, k_s, 
                                    k_oo, k_io, k_ii)
        #td = TwiceDifferentiable(energy, cells; autodiff = :forward)                            
        res = optimize(energy, cells, LBFGS(),
                    Optim.Options(iterations = 10^5))         
        cells = res.minimizer                   

        # Rotate and shift the Focal adhesions       
        rotation_and_shift!(cells, cparams, Adh)  
        
        # store the data
        CELLS[2*t, :] = cells 
        CPARAMS[2*t, :] = cparams 
        ADH[2*t, :, :, :] = Adh[:, :, :]
        ADH0[2*t, :, :, :] = Adh0[:, :, :]
    end

    return CELLS, CPARAMS, ADH, ADH0
end

function save_data(CELLS, CPARAMS, ADH, ADH0)

    npzwrite("DATA.npz", CELLS=CELLS, CPARAMS=CPARAMS, 
            ADH=ADH, ADH0=ADH0)
end    

CELLS, CPARAMS, ADH, ADH0 = main()
save_data(CELLS, CPARAMS, ADH, ADH0)
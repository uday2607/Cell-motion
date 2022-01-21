include("Adh_funcs.jl")
using Optim
using NPZ
using PyCall
using CMAEvolutionStrategy


function main()

    a = 6.0
    b = 3.0
    L = 40.0
    Num = 8
    Nadh = 64
    k_plus = 0.5
    dt = 1
    lamda = 0.5
    tau = 30
    T_S = 200.0
    k_s = 0.5
    k_m = 1.0
    k_oo = 5.0
    k_io = k_oo*100.0
    k_ii = k_io*100.0
    fThreshold = 0.5
    k_b = 0.005
    k_f = 0.0005
    alpha = 10.0
    TIME = 300
    a_min = a
    for i in 1:tau
        a_min -= a_min*lamda*dt/tau
    end    

    # Spwan cells and adhesions
    cells, cparams, Ovlap_indices = linear_preset(a, b, Num)
    Adh0, Adh, cAdh0, cp0 = random_adhesions(cells, cparams, Nadh, k_plus, dt)

    CELLS = zeros(TIME, 3*Num)
    CPARAMS = zeros(TIME, 4*Num)
    ADH = zeros(TIME, Num, Nadh, 2)
    ADH0 = zeros(TIME, Num, Nadh, 2)

    for t in 1:(TIME) 

        println("Time = ", t)

        # store the data
        CELLS[t, :] = cells 
        CPARAMS[t, :] = cparams 
        ADH[t, :, :, :] = Adh[:, :, :]
        ADH0[t, :, :, :] = Adh0[:, :, :]

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

                p_flag = protrusion!(cells, num, Adh, Adh0, cparams, Ovlap_indices, T_S, lamda, 
                            k_s, k_oo, k_io, k_ii, dt, tau, a)

                # Old bonds detach -> New bonds form  
                if cparams[4*num] != 0
                    detach!(cells, num, cAdh0, cp0, Adh, Adh0, 
                            k_b, k_f, 2.0*alpha, dt)    
                end                 

                # only when the cell protruded, form new adhesions 
                if p_flag == 1
                    # New adhesions at a smaller rate 
                    one_cell_random_adh!(cells, num, cparams, Adh, Adh0, cAdh0, cp0,
                    Nadh, k_plus/10, dt, 1)
                end
            
            else 
                # New bonds 
                println("New adhesions #", num)

                one_cell_random_adh!(cells, num, cparams, Adh, Adh0, cAdh0, cp0,
                                    Nadh, k_plus, dt, 0)
                cparams[4*num] = 1
            end     
        end

        # Energy minimization
        println("Energy: ", total_energy(cells, cparams, Ovlap_indices, Adh, Adh0, k_s, 
        k_oo, k_io, k_ii))
        min_bounds = Float64[]
        max_bounds = Float64[]
        for num in 1:Num
            push!(min_bounds, cells[3*(num-1)+1]-4.0*abs(cparams[3*(num-1)+1]*cos(cparams[3*(num-1)+3]))-eps(1e4))
            push!(max_bounds, cells[3*(num-1)+1]+4.0*abs(cparams[3*(num-1)+1]*cos(cparams[3*(num-1)+3]))+eps(1e4))
            push!(min_bounds, cells[3*(num-1)+2]-4.0*abs(cparams[3*(num-1)+1]*sin(cparams[3*(num-1)+3]))-eps(1e4))
            push!(max_bounds, cells[3*(num-1)+2]+4.0*abs(cparams[3*(num-1)+1]*sin(cparams[3*(num-1)+3]))+eps(1e4))
            push!(min_bounds, -pi)
            push!(max_bounds, pi)
        end
        #Julia -> SAMIN                                
        """res = optimize((cells) -> total_energy(cells, cparams, Ovlap_indices, 
                    Adh, Adh0, k_s, k_oo, k_io, k_ii), min_bounds, max_bounds,
                    cells, SAMIN(rt=0.9, f_tol = 1e-5),
                   Optim.Options(iterations = 10^6))
        cells = res.minimizer"""
        #Julia -> LBFGS                            
        """res = optimize((cells) -> total_energy(cells, cparams, Ovlap_indices, 
        Adh, Adh0, k_s, k_oo, k_io, k_ii), cells, LBFGS(),
                Optim.Options(iterations = 10^5))                 
        cells = res.minimizer"""
        #Scipy -> LBFGS
        res = so.minimize((cells) -> total_energy(cells, cparams, Ovlap_indices, 
                    Adh, Adh0, k_s, k_oo, k_io, k_ii), x0=cells, method="L-BFGS-B",
                    jac="3-point", options=Dict("maxfun"=>10^5, "maxls"=>50))
        println(res["success"])            
        cells = res["x"]
        #Scipy -> Differential Evolution
        """res = so.differential_evolution((cells) -> total_energy(cells, cparams, Ovlap_indices, 
                    Adh, Adh0, k_s, k_oo, k_io, k_ii), maxiter=10^3, bounds=cat(min_bounds, max_bounds, dims=2), 
                    disp=1, x0=cells, updating="deferred", polish=1)
        println(res["success"])            
        cells = res["x"]"""
        #Scipy -> Basin Hopping
        """minimizer_kwargs = Dict([("method", "L-BFGS-B")])
        res = so.basinhopping((cells) -> total_energy(cells, cparams, Ovlap_indices, 
                    Adh, Adh0, k_s, k_oo, k_io, k_ii), x0=cells, minimizer_kwargs=minimizer_kwargs,
                    niter=2000, niter_success=100, disp=true, stepsize=0.1, interval=10, T=0)
        println(res["message"])            
        cells = res["x"]"""
        """function f(cells)
            energy = total_energy(cells, cparams, Ovlap_indices, Adh, Adh0, k_s,
                        k_oo, k_io, k_ii)
            return energy            
        end     
        res = minimize(f=f, x0=cells, s0=1;
        lower = min_bounds, upper = max_bounds, multi_threading=true)
        cells = xbest(res)"""

        # Rotate and shift the Focal adhesions       
        rotation_and_shift!(cells, cparams, Adh)
        println("Energy: ", total_energy(cells, cparams, Ovlap_indices, Adh, Adh0, k_s, 
                            k_oo, k_io, k_ii))
    end

    return CELLS, CPARAMS, ADH, ADH0
end

function save_data(CELLS, CPARAMS, ADH, ADH0)

    npzwrite("DATA.npz", CELLS=CELLS, CPARAMS=CPARAMS, 
            ADH=ADH, ADH0=ADH0)
end    

CELLS, CPARAMS, ADH, ADH0 = main()
save_data(CELLS, CPARAMS, ADH, ADH0)
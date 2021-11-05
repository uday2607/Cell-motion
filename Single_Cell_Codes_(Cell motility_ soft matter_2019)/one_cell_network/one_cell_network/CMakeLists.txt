cmake_minimum_required(VERSION 3.12)
project(one_cell_network)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(one_cell_network
        brent.h
        contraction.h
        dbrent.h
        Deformations.h
        dF1dim.h
        dlinmin.h
        energy_grad.h
        extension.h
        f1dim.h
        frprmn.h
        ifcell.h
        in_out.h
        linmin.h
        main.cpp
        mnbrak.h
        Network.h
        nrutil.h
        rotation.h
        vector.h)

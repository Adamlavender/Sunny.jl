using Sunny, GLMakie

units = Units(:meV, :angstrom)
J1 = 1.0


function anneal!( sys , kTs , nsweeps , dt ; damping = 0.1)
    sampler = Langevin( dt ; damping , kT = kTs[ end ])
    Es = zeros(length( kTs ) ) # Buffer for saving energy as we proceed
    for (i , kT ) in enumerate( kTs )
        sampler.kT = kT
        for j in 1: nsweeps
        step!( sys , sampler ) # Perform a single step of the Langevin dynamics
            if j % 100 == 0 # Print progress every 100 sweeps
                println("Annealing step $j at kT = $kT")
            end
        end
        Es[ i ] = energy( sys ) # Query the energy
    end
    return Es # Return the energy values collected during annealing
end

J2 = 1.0
# Build FCC crystal and system
fcc = Sunny.fcc_crystal()
fig1 = view_crystal(fcc) # View the crystal structure
display(fig1) 
ratios = [0.25, 0.5, 0.75] #Get J2 ratios
sizes = 1:10

for J2_ratio in ratios
    J2 = J2_ratio * J1
    println("\nTesting J2/J1 = $J2_ratio")
    for n in sizes
        sys = System(fcc, [1 => Moment(s=1, g=2)], :dipole)
        set_exchange!(sys, J1, Bond(1,2,(0,0,0))) #set nearest neighbour coupling
        set_exchange!(sys, J2, Bond(1,1,(1,0,0))) #set next nearest neighbour coupling
        sys = repeat_periodically(sys, (n, n, n)) # change system size
        randomize_spins!(sys)
        minimize_energy!(sys) # minimize energy to find ground state
        e = energy_per_site(sys)
        println("  n = $n, energy per site = $e")
       
    end
end

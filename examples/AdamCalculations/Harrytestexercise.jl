using Sunny, GLMakie
using LinearAlgebra

units = Units(:meV, :angstrom)
J1 = 1.0

function anneal!( sys , kTs , nsweeps , dt ; damping = 0.1)
    sampler = Langevin( dt ; damping , kT = kTs[ end ])
    Es = zeros(length( kTs ) ) # Buffer for saving energy as we proceed
    for (i , kT ) in enumerate( kTs )
        sampler.kT = kT
        for j in 1: nsweeps
            step!( sys , sampler )
        end
        Es[ i ] = energy( sys )
    end
    return Es
end

fcc = Sunny.fcc_crystal()
fig1 = view_crystal(fcc)
ratios = [0.25,0.5,0.75]
sizes = 1:4

energy_vs_n = Dict{Float64, Vector{Float64}}()

for J2_ratio in ratios
    J2 = J2_ratio * J1
    energies = Float64[]
    println("\nTesting J2/J1 = $J2_ratio")
    for n in sizes
        sys = System(fcc, [1 => Moment(s=1, g=2)], :dipole)
        set_exchange!(sys, J1, Bond(1,2,[0,0,0]))
        set_exchange!(sys, J2, Bond(1,1,[1,0,0]))
        sys = repeat_periodically(sys, (n, n, n))
        randomize_spins!(sys)
        kTs = range(2.0, 0.01, length=20)
        anneal!(sys, kTs, 500, 0.05)
        minimize_energy!(sys)
        fig_2 = plot_spins(sys; color=[S[3] for S in sys.dipoles])
        display(fig_2)
        e = energy_per_site(sys)
        push!(energies, e)
        println("n = $n, Dominant wavevectors at J_2 ratio = $J2_ratio")
        print_wrapped_intensities(sys)
   
    end
    energy_vs_n[J2_ratio] = energies
end

#using results of this loop to then find smallest possible system of
#that accurately describes the ground state
println("\nSuggested magnetic supercells for J2 ratios:")
#for J_2 = 0.25, at n=2 we get 99.6 % weight to [1/2,0,0], at larger even n, this holds --> ground state, domains form for odd n
println("J2/J1 = 0.25:")
suggest_magnetic_supercell([[1/2,0,0]])
#for J_2 = 0.5, we already get two domains at n=2, choose the dominant wavevector
#The fact that the periodicity of the supercell matches the periodicty of the domain allows for a ground state to form.
println("J2/J1 = 0.5:")
suggest_magnetic_supercell([[0,0,1/2]])
#for J_2 = 0.75, at n=2 we get 100% weight to [1/2,1/2,1/2], this holds at large even n, domains form for odd n
println("J2/J1 = 0.75:")
suggest_magnetic_supercell([[1/2,1/2,1/2]])

#Check that the choice of supercell allows us to find a ground state for J_2 = 0.5
# Step 1: Get the suggested supercell for J2 = 0.5
dims = (1, 1, 2)  # For a diagonal supercell, e.g., (1, 1, 2)

# Step 2: Build and simulate the system
sys = System(fcc, [1 => Moment(s=1, g=2)], :dipole)
set_exchange!(sys, J1, Bond(1,2,[0,0,0]))
set_exchange!(sys, 0.5*J1, Bond(1,1,[1,0,0]))  # J2 = 0.5*J1
sys = repeat_periodically(sys, dims)
randomize_spins!(sys)
kTs = range(2.0, 0.01, length=20)
anneal!(sys, kTs, 500, 0.05)
minimize_energy!(sys)
println("Ground state for J_2=0.5")
print_wrapped_intensities(sys)




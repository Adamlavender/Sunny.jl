using Sunny
using GLMakie
using StaticArrays

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
 """
    Calculate the spin wave dispersion along a given path in the Brillouin zone.
    
    Parameters:
    - `lswt`: An instance of `SpinWaveTheory` representing the system.
    - `path`: A `BZPath` object defining the path in the Brillouin zone.
    
    Returns:
    - A dictionary containing the dispersion data.
"""
function spinwave_dispersion(system,path_info, J2; ordered=true)
    # To get the q-values for these high-symmetry points:
    first_path_labels = path_info.paths[1]
    first_path_qs = [path_info.points[label] for label in first_path_labels]
    path = q_space_path(system.crystal, first_path_qs, 500)
    if ordered
        lswt = SpinWaveTheory(system; measure=Sunny.ssf_perp(system))
        res = intensities_bands(lswt, path)
    else
        measure = Sunny.ssf_perp(system)
        kernel = lorentzian(fwhm=0.4)
        energies = range(0.0, 3.0, 150)
        lswt = SpinWaveTheoryKPM(system;measure,tol=0.05)
        res = intensities(lswt, path; energies, kernel)
    end
    
    plot_intensities(res; units, title = "Spin Wave dispersion for J_2=$J2")


end



function find_ground_state(sys;minimize=true)
    randomize_spins!(sys)
    kTs = range(2.0, 0.3, length=40)
    anneal!(sys, kTs, 1000, 0.02)
    minimize_energy!(sys, maxiters=10_000)

    return sys
end

# --- For a ground state with J2/J1 = 0.25, 
units = Units(:meV, :angstrom)
# --- Set up your system (example for J2/J1 = 0.25) ---
fcc = Sunny.fcc_crystal()
sys = System(fcc, [1 => Moment(s=1, g=2)], :dipole)
path_info = Sunny.irreducible_bz_paths(sys.crystal) #This returns two paths of high symmetry
J1 = 1.0   
n=12
shape1 = [2 0 0; 0 1 0; 0 0 1]
shape2 = [n 0 0; 0 n 0; 0 0 n]
shape3 = [1 0 1; 0 1 0; -1 -1 1]


set_exchange!(sys, J1, Bond(1,2,[0,0,0]))
set_exchange!(sys, 0.25, Bond(1,1,[1,0,0]))
sys1 = reshape_supercell(sys, shape1)  # Use minimal supercell for the order

set_exchange!(sys, J1, Bond(1,2,[0,0,0]))
set_exchange!(sys, 0.5, Bond(1,1,[1,0,0]))
sys2 = reshape_supercell(sys, shape2)  # Use minimal supercell for the order

set_exchange!(sys, J1, Bond(1,2,[0,0,0]))
set_exchange!(sys, 0.75, Bond(1,1,[1,0,0]))
sys3= reshape_supercell(sys, shape3)  # Use minimal supercell for the order

gs1 = find_ground_state(sys1)  # Example for J2/J1 = 0.25
gs3 = find_ground_state(sys3)  # Example for J2/J1 = 0.75
#spinwave_dispersion(sys3,path_info,0.75)
#spinwave_dispersion(sys1,path_info,0.25)

# Now for the disordered case, we can use the KPM method
# Use a large cell where n=8

gs2 = find_ground_state(sys2) # Example for J2/J1 = 0.5
#spinwave_dispersion(sys2,path_info,0.5; ordered=false)

#SCGA stuff
dq = 0.02  # grid resolution in RLU
kT = 1.0   # temperature (choose as appropriate for your study)
# Example for J2/J1 = 0.25
J2=0.5
sys_scga = System(fcc, [1 => Moment(s=1, g=2)], :dipole)
set_exchange!(sys_scga, J1, Bond(1,2,[0,0,0]))
set_exchange!(sys_scga, J2, Bond(1,1,[1,0,0]))
measure = ssf_perp(sys_scga)
scga = SCGA(sys_scga;measure,kT,dq)
grid = q_space_grid(fcc, [1, 0, 0], range(-10, 10, 200), [0, 1, 0], (-10, 10))

# Calculate and plot the instantaneous structure factor on the slice. Selecting
# `saturation=1.0` sets the color saturation point to the maximum intensity
# value. This is reasonable because we are above the ordering temperature and
# do not have sharp Bragg peaks.

res = intensities_static(scga, grid)
plot_intensities(res; saturation=1.0, title="Static Intensities at J2=$J2")


#Repeating parts 6 and 7 but with LandauLifshit
#Static intensity (static structure factor slice)
J2=0.25
expanded_gs = repeat_periodically(gs1, (8, 8, 8))
langevin = Langevin(; damping=0.2, kT=16*units.K)
suggest_timestep(expanded_gs, langevin; tol=1e-2)
langevin.dt = 0.025;
energies = [energy_per_site(expanded_gs)]
for _ in 1:1000
    step!(expanded_gs, langevin)
    push!(energies, energy_per_site(expanded_gs))
end

lines(energies, color=:blue, figure=(size=(600,300),), axis=(xlabel="Timesteps", ylabel="Energy (meV)"))
formfactors = [1 => FormFactor("Co2")]
measure = ssf_perp(expanded_gs; formfactors)
sc = SampledCorrelationsStatic(expanded_gs; measure)
add_sample!(sc, expanded_gs)

# Collect 20 additional samples. Perform 100 Langevin time-steps between
# measurements to approximately decorrelate the sample in thermal equilibrium.

for _ in 1:20
    for _ in 1:100
        step!(expanded_gs, langevin)
    end
    add_sample!(sc, expanded_gs)
end

grid = q_space_grid(fcc, [1, 0, 0], range(-10, 10, 200), [0, 1, 0], (-10, 10))
res = intensities_static(sc, grid)
plot_intensities(res; saturation=1.0, title="Static Intensities at J=$J2 with Landau Lifshitz")

#Energy momentum plot with landau lifshitz using dynamic structure factor
dt = 2*langevin.dt
energies = range(0, 6, 50)
sc = SampledCorrelations(expanded_gs; dt, energies, measure)
for _ in 1:5
    for _ in 1:100
        step!(expanded_gs, langevin)
    end
    add_sample!(sc, expanded_gs)
end
first_path_labels = path_info.paths[1]
first_path_qs = [path_info.points[label] for label in first_path_labels]
path = q_space_path(fcc, first_path_qs, 500)
res = intensities(sc, path; energies, langevin.kT)
plot_intensities(res; units, title="Intensities at J2=$J2 with Landau Lifshitz")
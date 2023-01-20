# Structure Factor Calculations

A dynamical structure factor gives a basic characterization of a spin system's
dynamical properties and is of fundamental importance when making comparisons
between theory and experimental scattering data. More specifically, it is a
function containing information about dynamical spin correlations, typically
written:

```math
𝒮^{αβ}_{jk}(𝐪, ω).
```
Given wave vector ``𝐪``, a frequency ``ω``, basis (atom) indices ``j`` and
``k``, and spin components ``α`` and ``β``, the dynamical structure factor will
yield a complex value.

Calculating the structure factor is relatively involved process. Sunny
provides a number of tools to facilitate the calculation and to extract
information from the results. These tools are briefly outlined below. "Real
life" use cases can be found in our tutorials and detailed function information
is available in the Library API.


## Basic Usage

### Calculating a dynamical stucture factor

The basic function for calculating dynamical structure factors is
[`calculate_structure_factor`](@ref). The steps for using it effectively are the
following:

1. Build a [`System`](@ref) and ensure that it is properly equilibrated at the
   temperature you wish to study. For example, if your `System` is in a ground
   state, one could use a [`LangevinHeunP`](@ref) integrator to thermalize it.
2. Set up a [`LangevinSampler`](@ref) that will generate decorrelated samples of
   spin configurations at the desired temperature.
3. Call `calculate_structure_factor(sys, sampler; kwargs...)`, which will return
   return a `StructureFactor`, containing all ``𝒮^{αβ}_{jk}(𝐪, ω)`` data.

The calculation can be configured in a number of ways, and we encourage you to
see the [`calculate_structure_factor`](@ref) documentation for a list of all
keywords. In particular, the user will likely want to specify the energy range (`ωmax`)
and resolution (`numω`) as well as the number of samples to calculate (`numsamps`).

### Extracting information

The basic function for extracting information from a `StructureFactor` at a
particular wave vector, ``𝐪``, is [`get_intensities`](@ref). It takes a
`StructureFactor` and either a single wave vector or an array of wave vectors.
For example: `get_intensities(sf, [0.0, 0.5, 0.5])`. Note that the wave vector
is specified in terms of reciprocal lattice units, though an alternative basis
may be specified by providing a transformation matrix to the keyword `newbasis`.

`get_intensities` will return a vector of intensities at different ``ω``s. The
precise ``ω`` values corresponding to each index can be retrieved by calling
`ωvals(sf)`, where `sf` is your `StructureFactor`.

Recall that the full structure contains a number of indices:
``𝒮^{αβ}_{jk}(𝐪,ω)``, but `get_intensities` only returns information
corresponding to ``ω``. By default, Sunny traces out the spin component indices
``α`` and ``β``. This behavior can be changed with the keyword argument
`contraction`. In addition to `:trace`, one may set `contraction=:depolarize` to
apply polarization corrections, or `contraction=(α,β)` to retrieve a particular
matrix element. The basis indices ``j`` and ``k`` are always reduced to out
through a phase averaging procedure. Note that information pertaining to these
indices is generally not accessible to experimental inquiry.

Since Sunny currently only calculates the structure factor on a finite lattice,
it is important to realize that exact information is only available at a
discrete set of wave vectors. Specifically, for each axis index $i$, we will get
information at ``q_i = \frac{n}{L_i}``, where $n$ runs from
``(\frac{-L_i}{2}+1)$ to $\frac{L_i}{2}`` and ``L_i`` is the linear dimension of
the lattice used for the calculation. If you request a wave vector that does not
fall in this set, Sunny will automatically round to the nearest ``𝐪`` that is
available. If `get_intensities` is given the keyword argument
`interpolation=:linear`, Sunny will use trilinear interpolation to determine the
results at the requested wave vector. 

To retrieve the intensities at all wave vectors for which there is exact data,
one can use the function [`intensity_grid`](@ref). This takes an optional
keyword argument `bzsize`, which must be given a tuple of three integers
specifying the number of Brillouin zones to calculate, e.g., `bzsize=(2,2,2)`. 

To calculate the intensities along a particular path, one may use the function
[`path`](@ref). This takes two arguments: a structure factor and a list of of
wave vectors. For example, `path(sf, [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0), (0.5,
0.5, 0.0)])`. `path` will return energy intensities along a path connecting
these points. The number of wave vectors sampled along the path is set with the
keyword `density`, which determines the number of wave vectors per inverse angstrom.

Note that all of these functions share keywords with [`get_intensities`](@ref).
In particular, they all take the keyword `kT` to set the temperature. It is
generally recommended to provided a value to `kT` corresponding to the
temperature at which measurements were taken. This allows Sunny to apply a
classical-to-quantum rescaling of the energy intensities. 

### Static structure factors

Static structure factors may be calculated simply by summing over all the
energies (i.e., the ``ω``-axis) provided by `get_intensities`. We recommend
calculating static structure factors from dynamical structure factors in this
way (rather than directly from a series of equilibrium samples). This approach
makes it possible to apply the classical-to-quantum intensity rescaling, which
is energy dependent.

For convenience, Sunny provides the function [`get_static_intensities`](@ref), which
will perform the summation for you.
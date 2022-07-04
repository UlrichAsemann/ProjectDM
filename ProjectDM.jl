#========================================================================================#
"""
	NBodies - Version 5

A system of N interacting bodies moving in dim dimensions over time T in res timesteps.

Author: Niall Palfreyman, 29/05/2022.
"""

# Additional code written by Ulrich Asemann


module NBodies

using GLMakie
using LinearAlgebra

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	NBody

An NBody system capable of containing multiple (N) bodies that gravitationally interact
with each other.
"""
mutable struct NBody
	N							# Number of bodies
	nsteps						# Duration of simulation
	dt							# Timestep resolution
	x0							# Initial positions of bodies
	p0							# Initial momenta of bodies
	m							# Masses of bodies
	x							# Current position of system
	p							# Current momentum of system
	e							# Current energy of system
	size						# Size of the body

	"Construct a new NBody"
	function NBody( T=40, resolution=20000)
		# Initialise all fields of the decoding apparatus:
		new(
			0,					# Initially no bodies in the system
			resolution,			# Duration
			T/resolution,		# Timestep resolution
			[],					# Initial positions
			[],					# Initial momenta
			[],					# Masses of bodies
			[],					# Current position empty
			[],					# Current momentum empty
			[],					# Current energy empty
			[]					# Size of body
		)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	addbody!( nbody::NBody, x0::Vector{Float64}, p0::Vector{Float64}, m::Float64=1)

Add to the system a new body with initial position and momentum x0, p0, with mass m and the diameter size.
The new structure elements nb.x and nb.p represent the current position and momentum of the
bodies, which will be constantly updated during the simulation.
"""
function addbody!( nbody::NBody, x0::Vector{Float64}, p0::Vector{Float64}, m::Float64=1.0, size::Float64=1.0)
	push!( nbody.x0, x0)
	push!( nbody.x, deepcopy(x0))
	push!( nbody.p0, p0)
	push!( nbody.p, deepcopy(p0))
	push!( nbody.m, m)
	push!( nbody.size, size)
	nbody.N  += 1
end

#-----------------------------------------------------------------------------------------
"""
	simulate( nb::NBody)

Run a simulation of the given NBody system over duration nb.nsteps * nb.dt.
"""
function simulate( nb::NBody)
	t = 0:nb.dt:nb.nsteps*nb.dt
	x = Vector{typeof(nb.x0)}(undef,nb.nsteps+1)
	p = Vector{typeof(nb.p0)}(undef,nb.nsteps+1)
	e = Vector{Float64}(undef,nb.nsteps+1)

	# Initialisation:
	x[1] = nb.x0
	p[1] = nb.p0
	e[1] = energy_calc(nb)			# using energy_calc() to calculate the energy of the system

	# Simulation using RK4 method and energy-calc():
	for n = 1:nb.nsteps
		rk4!(nb)
		x[n+1] = nb.x
		p[n+1] = nb.p
		e[n+1] = energy_calc(nb)
	end

	(t,x,e,p)
end

#-----------------------------------------------------------------------------------------
"""
	rk4!( nbody::NBody)

	Perform a single RK4-step of the given NBody system.

	Calculating the different slopes k1, k2, k3 and k4 to calcualte the slope for the current value.
	Formula for x: x[n+1] = x[n] + Δt * 1/6 * (k1 + 2*k2 + 2*k3 + k4)
"""
function rk4!( nb::NBody)

	# Half-timestep:
	dt2 = nb.dt/2

	# Evaluation of k1:
	xk1 = nb.x
	pk1 = nb.p

	# Evaluation of k2:
	xk2 = nb.x + dt2 * pk1./nb.m
	pk2 = nb.p + dt2 * forceOnMasses(xk1 ,nb.m)

	# Evaluation of k3:
	xk3 = nb.x + dt2 * pk2./nb.m
	pk3 = nb.p + dt2 * forceOnMasses(xk2, nb.m)

	# Evaluation of k4:
	xk4 = nb.x + nb.dt * pk3./nb.m
	pk4 = nb.p + nb.dt * forceOnMasses(xk3, nb.m) 
	
	# Evaluation of full-step:
	nb.x = nb.x + nb.dt * 1/6 * (pk1./nb.m + 2pk2./nb.m + 2pk3./nb.m + pk4./nb.m)
	nb.p = nb.p + nb.dt * 1/6 * (forceOnMasses(xk1, nb.m) + 2*forceOnMasses(xk2, nb.m) + 2*forceOnMasses(xk3, nb.m) + forceOnMasses(xk4, nb.m))

end

#-----------------------------------------------------------------------------------------
"""
	forceOnMasses( locations::Vector, masses::Vector)

	Internal utility function: Calculate all gravitational forces between the N bodies with
	the given mass at the given locations.
	This method is based on the following Newtonian formula:

	G =  6,67259 * 10^-11(N*m^2)/kg^2
	Using the gravitational constant for a near real system.

	F_ij = Force on i-th body due to j-th body = - (G m_i m_j / |r_i-r_j|^3) * vec(r_i-r_j)
"""
function forceOnMasses( locations::Vector, masses::Vector)
	G = 6.67259e-11										# Newton's gravitational constant
	gmm = G * masses * masses'							# Precalculate products between the masses

	locnPerBody = repeat(locations,1,length(locations))	# Matrix of body locations beside each other
	relpos = locnPerBody -								# For all bodies i and j, calculate the
					permutedims(locnPerBody)			# relative position vector r_ij between them.
	invCube = abs.(relpos'.*relpos) .^ (3/2)			# Calculate 1 / (r_ij)^3 between all bodies
	for i in 1:length(locations) invCube[i,i] = 1 end	# Prevent zero-division along diagonal

	forceFromMasses = -gmm .* relpos ./ invCube			# Calculate force contribution FROM each source

	vec(sum( forceFromMasses, dims=2))					# Calculate sum of all forces ON each source.
end

#-----------------------------------------------------------------------------------------
"""
	relpos(locations::Vector)

	Function to calculate the relative positions from each body to another

"""
function relpos(locations::Vector)

	locnPerBody = repeat(locations,1,length(locations))	
	relpos = locnPerBody - permutedims(locnPerBody)			

	relpos
end

#-----------------------------------------------------------------------------------------
"""
	animate( nb::NBody, t, x, e)

	Animate the simulation data (t,x,e) of the given NBody system.
	This implementation is identical to that of version 3.
"""
function animate( nb::NBody, t, x, e)
	# Prepare an Observable containing the initial x/y coordinate and energy for each body:
	x_current = Observable( map( bodycoords->bodycoords[1], x[1]))
	y_current = Observable( map( bodycoords->bodycoords[2], x[1]))
	timestamp = Observable( string( "t = ", round(t[1], digits=2)))
	energy  = Observable( string( "System-Energy = ", e[1]))

	# Variables for the animation
	limits = 8e8
	text_pos = 0.9*limits
	animation_size = 900

	# Prepare the animation graphics:
	fig = Figure(resolution=(animation_size, animation_size))
	ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "N-body 2D Motion")
	limits!( ax, -limits, limits, -limits, limits)
	scatter!( ax, x_current, y_current, markersize=((animation_size/limits)*nb.size), color=:blue) 	# Calculating the size of each body
																									# depending on the size of the animation
	text!( timestamp, position=(-text_pos, -text_pos), textsize=30, align=(:left,:center))
	text!( energy, position=(-text_pos, text_pos), textsize=20, align=(:left,:center))
	display(fig)


	# Run the animation:
	for i in 1:nb.nsteps+1
		x_current[] = map( bodycoords->bodycoords[1], x[i])
		y_current[] = map( bodycoords->bodycoords[2], x[i])
		timestamp[] = string( "t = ", round(t[i], digits=2))
		energy[]	= string( "System-Energy = ", e[i])

		sleep(1e-4)
	end
end

#-----------------------------------------------------------------------------------------
"""
	energy_calc(nb::NBody)

	Evaluating the potential energy in the system to check, if the energy stays the same, so the animation is realistic or not
	Values can fluctuate a bit, due to the calculation with big timesteps
	Formula: EPot = G * M * m * (1/r1 - 1/r2)
"""
function energy_calc(nb::NBody)

	# Gravitational constant
	G = 6.67259e-11		

	# Precalculating -(G * M * m)
	gmm = (G * nb.m * nb.m')

	# Calculating the relative positions of each body to another
	rel_pos = relpos(nb.x)

	# Calculating the magnitude of the distance between each body
	xvalue_ij = sqrt.(abs.(rel_pos'.*rel_pos))              

	# Preventing zero-division along diagonal
	for i in 1:length(nb.x) xvalue_ij[i,i] = 1 end        

	# Calculating the matrix with the potential energy
	epot_values = gmm ./ xvalue_ij                        

	# Adding up the potential energy of each body to another
	systemenergy = sum(triu(epot_values))

	# Returning the total systemenergy
	systemenergy
	
end

end		# of NBodies
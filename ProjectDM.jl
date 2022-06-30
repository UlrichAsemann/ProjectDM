#========================================================================================#
"""
	NBodies - Version 5

A system of N interacting bodies moving in dim dimensions over time T in res timesteps.
I noticed in version 4 that the Sun and planet gradually move further away from each other,
which of course is to be expected if I use Euler's method, which always errs outward from
any curved trajectory. So in version 5, I need to replace Euler's method by the Runge-Kutta-2
method. And this in turn means I must refactor the Euler code out into its own method, then
replace it by Runge-Kutta-2.

This solution worked so well that I then added a second, more complicated use-case (demo2()),
in which three bodies (2 Suns and a very lonely! planet) orbit around each other.

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
	function NBody( T=40, resolution=20000, G=1)
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

Add to the system a new body with initial position and momentum x0, p0, and with mass m. And size size.
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

	# Simulation using RK4 method:
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
	animate( nb::NBody, t, x, e)

	Animate the simulation data (t,x,e) of the given NBody system.
	This implementation is identical to that of version 3.
"""
function animate( nb::NBody, t, x, e)
	# Prepare an Observable containing the initial x/y coordinate and energy for each body:
	x_current = Observable( map( bodycoords->bodycoords[1], x[1]))
	y_current = Observable( map( bodycoords->bodycoords[2], x[1]))
	timestamp = Observable( string( "t = ", round(t[1], digits=2)))
	energy  = Observable( string( "energy = ", e[1]))

	# Prepare the animation graphics:
	fig = Figure(resolution=(900, 900))
	ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "N-body 2D Motion")
	limits!( ax, -1e9, 1e9, -1e9, 1e9)
	scatter!( ax, x_current, y_current, markersize=(nb.size), color=:blue)
	text!( timestamp, position=(-0.9e9, -0.9e9), textsize=30, align=(:left,:center))
	text!( energy, position=(-0.9e9, 0.9e9), textsize=20, align=(:left,:center))
	display(fig)


	# Run the animation:
	for i in 1:nb.nsteps+1
		x_current[] = map( bodycoords->bodycoords[1], x[i])
		y_current[] = map( bodycoords->bodycoords[2], x[i])
		timestamp[] = string( "t = ", round(t[i], digits=2))
		energy[]	= string( "energy = ", e[i])

		sleep(1e-4)
	end
end

#-----------------------------------------------------------------------------------------
"""
	energy(nb::NBody)

	Evaluating the energy in the system to check, if the energy stays the same, so the animation is realistic or not
	Formula: EKin = 1/2 * m * v^2
"""

function energy_calc(nb::NBody)

	# Transforming the momentum into velocity
	v = hcat(nb.p./nb.m)

	v_value = zeros(length(nb.p))  

	# Calculating the velocity for each body
	for i in 1:length(nb.p)
		
		v_value[i] = sqrt(sum(v[i].^2))
	end

	# Calculating the energy of the system
	systemenergy = sum(1/2*(nb.m.*(v_value.^2)))

	# Giving back the whole energy of the system
	systemenergy
end

#-----------------------------------------------------------------------------------------
"""
	demo()

	Demonstrate simulation of a simple 2-body problem in a simple 3-step use-case.
"""
function demo()
	# Build the 2-body system:
	nb = NBody( 20, 1000)								# 20 time units divided into 1000 steps
	addbody!( nb, [0.0, 1.0], [ 0.8,0.0], 2.0, 1.0)			# Sun (m = 2)
	addbody!( nb, [0.0,-1.0], [-0.8,0.0], 1.0, 1.0)				# Planet (m = 1)
	
	# Calculate the simulation data:
	(t,x,e) = simulate(nb)

	# Run the animation:
	animate( nb, t, x, e)
end

#-----------------------------------------------------------------------------------------
"""
	demo2()

	Demonstrate simulation of a chaotic 3-body system.
"""
function demo2()
	# Build the 3-body system:
	nb = NBody( 20, 1000)
	addbody!( nb, [0.0, 1.0],	[ 0.8, 0.0], 	2.0, 2.0)		# Sun 1 (m = 2)
	addbody!( nb, [0.0,-1.0],	[-0.8, 0.0],	1.0, 1.0)		# Sun 2 (m = 1)
	addbody!( nb, [3.0, 0.0],	[ 0.0, 0.1],	0.2, 0.2)		# Planet (m = 0.2)

	energy_calc(nb)
	
	# Run the simulation:
	t,x,e = simulate(nb)

	# Run the animation:
	animate(nb, t, x, e)	
end

function demo3()
	"""
	Demo der Erde
	Demo zur Implementierung der neuen physikalischen Größen
	Gravitatonskonstante etc.

	m_erde = 5,972 *10^24kg
	m_mond = 7,3483 * 10^22kg

	G =  6,67259 * 10^-11(N*m^2)/kg^2

	Distannz Mond zu Erde: 384 400km

	Größenverhältnis Erde zu Mond: 1 zu 3,7
	"""

	# Build the 3-body system:
	nb = NBody(10000000, 1000)

	m_sun = 1.989e30 #kg 
	v_moon = 1023 #m/s
	m_moon = 7.3483e22 #kg


	addbody!( nb, [ 0.0, 0.0],			[ 0.0, 	0.0], 			5.972e24, 		74.0)		# Earth
	addbody!( nb, [ 384400000.0, 0.0],	[ 0.0, m_moon*v_moon],	m_moon,      	20.0)		# Moon
	addbody!( nb, [ -384400000.0, 0.0],	[ 0.0, -m_moon*v_moon],	m_moon,      	20.0)
	
	energy_calc(nb)

	
	# Run the simulation:
	t,x,e = simulate(nb)

	# Run the animation:
	animate(nb, t, x, e)	
end

function demo4()
	# Build the 3-body system:
	nb = NBody( 20, 1000)
	addbody!( nb, [0.0, 1.0],	[ 0.8, 0.0], 	2.0)		# Sun 1 (m = 2)
	addbody!( nb, [0.0,-1.0],	[-0.8, 0.0],	1.0)		# Sun 2 (m = 1)
	addbody!( nb, [3.0, 0.0],	[ 0.0, 0.1],	0.2)		# Planet (m = 0.2)
	
	"""
	Testverfahren zur Implementierung von energy_calc und RK4


	nb.p./nb.m
	energy_calc(nb)

	sum(diag(diag(energy_calc(nb).*energy_calc(nb)')[3]))
	length(nb.p)
	a = (nb.p./nb.m)
	a[1]*a[1]'
	diag(a[1]*a[1]')
	a'
	b = diag(a*a')
	b[1][1]

	016*016
	64*64
	0.04*0.04
	0.25*0.25


	[[a[1] .* a[1]]; [a[2] .* a[2]]; [a[3] .* a[3]]]
	[a[1:length(a)].*a[1:length(a)]]
	"""


	# Run the simulation:
	#t,x = simulate(nb)

	# Run the animation:
	#animate(nb, t, x)	
end

end		# of NBodies
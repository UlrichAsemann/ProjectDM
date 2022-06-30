"""
    Hier wird das Projekt und die Beantwortung der Fragestellung in einer eigenen Datei beantwortet.
    Das Ganze ist work-in-progress!
"""


include("ProjectDM.jl")
using NBodies

nb = NBody(10000000, 1000)



"""
m_sun = 1.989e30 #kg 
v_moon = 1023 #m/s
m_moon = 7.3483e22 #kg


addbody!( nb, [ 0.0, 0.0],			[ 0.0, 	0.0], 			5.972e24, 		74.0)		# Earth
addbody!( nb, [ 384400000.0, 0.0],	[ 0.0, m_moon*v_moon],	m_moon,      	20.0)		# Moon
#addbody!( nb, [ -384400000.0, 0.0],	[ 0.0, -m_moon*v_moon],	m_moon,      	20.0)

energy_calc(nb)


# Run the simulation:
t,x,e = simulate(nb)

# Run the animation:
animate(nb, t, x, e)	
"""


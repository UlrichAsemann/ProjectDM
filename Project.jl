"""
    Hier wird das Projekt und die Beantwortung der Fragestellung in einer eigenen Datei beantwortet.
    Das Ganze ist work-in-progress!
"""


include("ProjectDM.jl")
using NBodies

# Build the 3-body system:
nb = Main.NBodies.NBody(100, 1000)

#m_sun = 1.989e30
#addbody!( nb, [ 0.0, 0.0],		[ 0.0, 0.0], 	5.972e24 )		# Earth
#addbody!( nb, [ 384400.0, 0.0],	[ 0.0, 2.5e27],	7.3483e22)		# Moon
#addbody!( nb, [ 0.0, 0.0],		[ 0.0, 0.0],	m_sun)			# Sun
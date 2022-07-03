"""
    Hier wird das Projekt und die Beantwortung der Fragestellung in einer eigenen Datei beantwortet.
    Das Ganze ist work-in-progress!
"""


include("ProjectDM.jl")
using .NBodies

"""
Fragestellung: 
Kann der Mond der Erde durch einen in das System eintretenden Körper ersetzt werden,
ohne das es zu einer Kollision kommt?
"""

"""
Schritt 1: Deklaration bekannter Variablen

Masse Mond, Masse Erde
Durchschnittliche Geschwindigkeit des Mondes um die Erde
Durchschnittliche Distanz zwischen Mond und Erde
Größenverhältnis von Mond und Erde: Erde 3,7 Mal so groß wie der Mond
"""

m_moon = 7.3483e22 #kg
m_earth = 5.972e24 #kg

v_moon_y = 1023 #m/s

d_moon_earth = 384400000.0 #m

size_earth = 3.7
size_moon = 1.0

"""
Schritt 2: Erstellen eines Systems

System mit Länge der Simulation und Anzahl der Simulationsschritte
"""

sim_length = 10000000
sim_steps  = 10000

nb = NBodies.NBody(sim_length, sim_steps)

"""
Schritt 3: Hinzufügen von Mond und Erde in das System
"""

NBodies.addbody!( nb, [ 0.0,            0.0],	[ 0.0,             0.0], 	m_earth, 20*size_earth)		# Earth
NBodies.addbody!( nb, [ d_moon_earth,   0.0],	[ 0.0, m_moon*v_moon_y],	m_moon,  20*size_moon)		# Moon

"""
Schritt 4: Hinzufügen eines dritten Körpers in das Systems

Deklaration der Variablen des neuen Körpers
Hinzufügen des Körpers in das System
"""

m_body = 0.0
v_body = [0.0, 0.0]
size_body = 0.0
d_body_earth = [0.0, 0.0]

#NBodies.addbody!( nb, d_body_earth, v_body, m_body, size_body)

"""
Schritt 5: Simulation und Animation
"""

# Simulation:
t,x,e = NBodies.simulate(nb)

# Animation:
NBodies.animate(nb, t, x, e)	



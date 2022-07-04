
include("ProjectDM.jl")
using .NBodies

"""
    Fragestellung: 

    Kann der Mond der Erde durch einen in das System eintretenden Körper ersetzt werden,
    ohne das es zu einer Kollision kommt?

    Fokus:
        - realistisches System
        - schnelles einstellen eines Gleichgewichtes
        - keine Kollision
"""


function demo()
"""
Schritt 1: Deklaration bekannter Variablen

    Masse Mond, Masse Erde
    Durchschnittliche Geschwindigkeit des Mondes um die Erde
    Durchschnittliche Distanz zwischen Mond und Erde
    Durchmesser von Erde und Mond
"""

    m_moon = 7.3483e22 #kg
    m_earth = 5.972e24 #kg

    v_moon = 1023 #m/s

    size_earth = 1.2742e7 #m  # diameter
    size_moon = 3.4748e6  #m  # diameter

    # Hinzufügen der Radien von Mond und Erde, da vom Zentrum der Himmelskörper berechnet wird
    d_moon_earth = 3.844e8 + 0.5*size_earth + 0.5*size_moon #m

"""
Schritt 2: Erstellen eines Systems

    System mit Länge der Simulation und Anzahl der Simulationsschritte
"""

    sim_length = 20000000
    sim_steps  = 10000

    nb = NBodies.NBody(sim_length, sim_steps)

"""
Schritt 3: Hinzufügen von Mond und Erde in das System
"""

    NBodies.addbody!( nb, [ 0.0,            0.0],	[ 0.0,           0.0], 	m_earth, size_earth)		# Erde
    NBodies.addbody!( nb, [ d_moon_earth,   0.0],	[ 0.0, m_moon*v_moon],  m_moon,  size_moon)		    # Mond

"""
Schritt 4:  Hinzufügen eines dritten Körpers in das Systems,
            der den Mond ersetzt                

    Deklaration der Variablen des neuen Körpers
    Hinzufügen des Körpers in das System
"""

m_body = 8e23
p_body = [-39.0*m_body, 600.0*m_body]
size_body = 2*size_moon
d_body_earth = [-6.8e8, -8.0e8]

NBodies.addbody!( nb, d_body_earth, p_body, m_body, size_body)

"""
Schritt 5: Simulation und Animation
"""

# Simulation:
t,x,e = NBodies.simulate(nb)

# Animation:
NBodies.animate(nb, t, x, e)	
end # of Demo
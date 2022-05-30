"""
Anfang des Projektes des 4. Semesters - NBodies Problem - Programmieren dynamischer Modelle

Frage:
Can the moon of the earth be replaced by another body without collsion and harm to the earth?

Fokus:
Effizienz z.B. Zeit oder wenig Chaos bis neues Gleichgewicht!
Kein Kollision miteinander!


Bedingungen:
Your project deliverable is a Julia module that must fulfil the following constraints:
	- the demonstration program MUST execute correctly within its own folder
	- simulates the motion of three or more gravitational bodies
	- executes an aesthetically pleasing chaotic trajectory
	- demonstrates that the motion conserves energy
	- demonstrates that the motion is indeed chaotic
	- displays explanatory text so that anyone can understand it WITHOUT studying the code!
	- the code must be neat, well documented and easily readable and comprehensible.

	- RK-4 verwenden oder höher (MUST HAVE) (eventuell RK-8)
	- konzeptuelle Gesamtgeschichte mit innovativem Aspekt

	- weitere Test

Informationen zum Projekt:

Video zu NBodies in Matlab:
	- https://www.youtube.com/watch?v=LwkvO3t1b30&t=113s
""";


"""
Ideen für Aufbau des Codes:

Benötigte Parameter für die Körper:
    - m: 		Masse (einzeln)
	- [x y]: 	Anfangskoordinaten (Vektor) (Arbeiten mit Matrixmanipulation)
	- v:		Anfangsgeschwindigkeit (Wie bewegt sich der Körper im System)
    
Benötigte Parameter für das System:
	- T:		Simulationslänge
	- Δt:		Schrittgröße zwischen den einzelnen Simulationsschritten 
	- n:		Anzahl an Körpern im System

Benötigte Methoden:
	- Berechnung der relativen Position zu den anderen Körpern
	- Berechnung der Anziehungskraft zueinander durch Newtons Gravitationsgesetz
		--> Dadurch Berechnung der neuen Positionen der Körper

Voraussetzungen:
	- Verwendung von GLMakie
	- Verwendung von RK-2 zur Berechnung der Simulationsschritte
	- Verwendung gegebener Massen von Mond und Erde
	- Finden der optimalen Parameter, um Frage zu beantworten!!!

Sonstiges (unwichtig am Anfang):
	- Überlegung der Fenstergröße der Simulation
	- Überlegung Aufbau der Simulation
	- Pfeile für Anzeige in welche Richtung das System sich bewegt	
""";
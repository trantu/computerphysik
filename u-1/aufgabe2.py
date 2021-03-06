import math
import numpy

'''
Aufgabe 1.2.1

Gleitkommaarithmetik: Operation der Gleitkommazahlen, bei der die Darstellung der Zahlen durch Angabe die Mantisse sowie das Exponent begrenzt ist. Davon hängt die Genauigkeit des Egebnisses ab. Gleitkommaarithmetik fordert in der Regel jenach Länge der Darstellung viel Rechenaufwand und Speicherbedarf, weil auch Zwischenergebnisse mit gespeichert werden. Ein Beispiel sind Berechnungen der unterschiedlichen Winkel oder Farbenänderungen bei den Grafikkarten, sodass dafür extra viel Speicher und leistungstarke Prozessoren verbaut werden.

Ausloeschung: hier spricht man von Rundungsfehlern bei der Subtraktion zweier Gleitkommazahlen, dass es entsprechende Wirkung zeigt, wenn es an einer bestimmten Stelle nach der Zahl mit kürzeren Stellen nach dem Komma  nach dem komma gerundet wird.
Bei der Reihenfolge der numerischen Sumation wird das Ergebnis beeinflusst, weil die Genauigkeit von der Zahlendarstellung abhaengt. Z.B. bei einer Reihe mit einer Wiederholung von 10 wird jedes Mal die Zahl 0.00001 summiert und wenn es ab der 4. Stelle nach dem KOmma gerundet werden soll, dann wird das Ergebnis nach der Summation gleich null sein. Wir wissen, wenn es z.B. ab 6. Stelle gerundet wird, dann erhalten wir 0.0001 als Ergebnis. Hier sieht man die Wirkung der numerischen Sumation.
'''
# mround wird benötigt um sehr kleine Zahl mit addieren zu können z.B. bei n = 6 und x = 0.000001
def mround ( x,N ) :
    if ( x==0 ) :
        return x
    return round ( x, int ( N - math.ceil ( math.log10 ( abs ( x )))))

def leibnitz(k, n):
    sum = 0.0
    rev_sum = 0.0
    single = numpy.single(0.0)
    
    # 0 bis k
    for i in range(0,k):
        sum = sum + mround(pow(-1,i) / (2.0*i+1),n)
    r = range(0,k)
    
    list.reverse(r) 
    
    # umgekerte Reihenfolge
    for i in r:
        rev_sum = rev_sum + mround(pow(-1,i) / (2.0*i+1),n)
    # single precision 
    for i in range(0,k):
        single = single + numpy.single(pow(-1,i) / (2.0*i+1))
        
    return (sum, rev_sum, single)    

def stagniert(n):
    return math.ceil(math.sqrt(1/n))

#Aufgabe 1.2.2
print "%g %g %g" % leibnitz(500,4)
#>>> 0.784842 0.784842 0.784898
print "%g %g %g" % leibnitz(500,6)
#>>> 0.784898 0.784898 0.784898

#Aufgabe 1.2.3
# n = 6
print stagniert(0.0000004)
#>>> k = 1582.0    

# n = 4
print stagniert(0.00004)
#>>> k = 159.0

__author__ = 'seb'


#!/usr/bin/python

from NumPy import *
import Gnuplot, Gnuplot.funcutils

g = Gnuplot.Gnuplot(debug=1)

a = Gnuplot.Data(([1,1], [2,2], [3,3]),title="A")
b = Gnuplot.Data(([1,1], [2,4], [3,9]),title="B")
c = Gnuplot.Data(([1,0.5], [2,1], [3,1.5]),title="C")


g('set output "/tmp/myGraph.png"')
g('set terminal png small ')
l = []
for i in (a,b,c):
    l.append(i)

g._add_to_queue(l)
g.replot()



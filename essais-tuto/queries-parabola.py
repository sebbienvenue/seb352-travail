#for queries examplefrom tutorial 
from sys import argv
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData
from aiida.orm.calculation import *

path="/home/aiida/Documents/seb352-travail/essais-tuto/res/"
StructureData = DataFactory("structure")
ParameterData = DataFactory("parameter")
#PwCalculation= DataFactory("calculation")

qb=QueryBuilder()



qb.append(ParameterData,
	project=["attributes.step0", "attributes.steps"],
	filters={"id":{"==":5615}}
	)

ok=0
a=qb.count()
file=open(path+"results-parabola-dict", 'w')



for i in qb.iterall():
	file.write("{}\n\n{}\n\n{}\n\n{}".format(i[0],i[1][0],i[1][1],i[1][2]))
	ok+=1
	print i[0]['dE']
	print len(i[0])

file.close()

new_dict={'0': i[0], '1': i[1][0], '2': i[1][1],'3' :i[1][2]}

#print new_dict
print

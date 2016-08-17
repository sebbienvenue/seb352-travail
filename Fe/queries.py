from sys import argv
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData
#from aiida.orm.calculation import Calculation, JobCalculation
import os


path=os.getcwd()
path=str(path)+"/"


StructureData=DataFactory("structure")
ParameterData=DataFactory("parameter")



qb=QueryBuilder()
f=open(path+"seb", 'w')
f.write("Element, Volume[a.u],     Energy[eV]\n")

qb.append(StructureData,
	project=["attributes.kinds"],
#	filters={"name":{"in":"attributes.kinds"}}
	tag="structure"
	)

qb.append(Calculation,
	#project=["id"],
	output_of="structure"
	)

qb.append(ParameterData,
	project=["attributes.volume", "attributes.energy",
		"attributes.energy_accuracy_units"],
	output_of="Calculation"
	)


ok=0

for i in qb.iterall():
	if( i[0][0]['symbols'][0] == 'Fe'):
		print i[0][0]['symbols'][0], i[1], i[2]
		ok+=1
		#print i
		f.write("{},\t {},\t {}\n".format(i[0][0]['symbols'][0], str(i[1]), str(i[2])))
	














f.close()
print "\nFinished\nNombre element in the query:  ", ok, "\n"

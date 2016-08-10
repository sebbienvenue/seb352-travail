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



qb.append

#essai une query sur PwClaculation
qb.append(Calculation,
        #filters={"id":{"==":4285}},
	tag="calculation",
	output_of="structure"
        )


#donne juste les nom qui sont dans les groupes
qb.append(Group,
	group_of="calculation",
        project=["name"],
	filters={"name":{"in": ["tutorial_pbesol", "tutorial_lda", "tutorial_pbe"]   }}
	)


a=qb.count()
for i in qb.iterall():
        print i


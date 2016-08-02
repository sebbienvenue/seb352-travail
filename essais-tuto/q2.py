#for queries examplefrom tutorial 

from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData

qb=QueryBuilder()
#qb.append(Node, project=["id"])
StructureData = DataFactory("structure")
ParameterData = DataFactory("parameter")

#enumerate the <pk> for each query key
#for node, in qb.iterall():
#	print node
#print
#print("Number of species "+str( qb.count()))

#qb.append(StructureData, project=["id", "uuid"], 
#	filters={"or":[
#	{"id":{"==":285}}, {"id":{"==":3512}} ] })


#	Pour etablir des liens entre etats
#qb.append(RemoteData, tag="remote", project=["*"])
#qb.append(Group, group_of="remote")

#qb.append(ParameterData, project=["attributes.energy_smearing"]) #, filters=)
qb.append(ParameterData, project=["attributes.element"])

for i in qb.iterall():
	print i



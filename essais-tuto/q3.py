#for queries examplefrom tutorial 
from sys import argv
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData
path="/home/aiida/Documents/seb352-travail/essais-tuto/res/"
StructureData = DataFactory("structure")
ParameterData = DataFactory("parameter")

qb1=QueryBuilder()
qb2=QueryBuilder()
qb3=QueryBuilder()
#qb.append(Node, project=["id"])

#enumerate the <pk> for each query key
#for node, in qb.iterall():
#	print node
#print
#print("Number of species "+str( qb.count()))

#qb.append(StructureData, project=["id", "uuid"], 
#	filters={"or":[
#	{"id":{"==":285}}, {"id":{"==":3512}} ] })


#	Pour etablir des liens entre etats
qb1.append(RemoteData, tag="remote", project=["*"])
qb1.append(Group, group_of="remote")

qb2.append(RemoteData,  project=["*"])

qb3.append(Group)


#qb.append(ParameterData, project=["attributes.energy_smearing"]) #, filters=)
#qb.append(ParameterData, project=["attributes.element"])

f1=open(path+"remoteData_Group", 'w')
f2=open(path+"remoteData", 'w')
f3=open(path+"Group", 'w')

for i in qb1.iterall():
	f1.write(str(i)+"\n")

for j in qb2.iterall():
	f2.write(str(j)+"\n")

for k in qb3.iterall():
	f3.write(str(k)+"\n")



f1.close()
f2.close()
f3.close()


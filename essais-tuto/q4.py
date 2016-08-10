#for queries examplefrom tutorial 
from sys import argv
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData
path="/home/aiida/Documents/seb352-travail/essais-tuto/res/"
StructureData = DataFactory("structure")
ParameterData = DataFactory("parameter")

qb=QueryBuilder()




#########	THis is working

qb.append(StructureData, 
	project=["extras.formula"],
#	filters={"extras.formula":"Sn2O3"},
	filters={"extras.formula":"LaCoO3"},
	tag="structure")

qb.append(
	Calculation,
	project=["id"],
	tag="calculation", 
	output_of="structure")

qb.append(ParameterData, 
	tag="results",
	filters={"attributes.energy_smearing":{"<=":-0.0001}},
	project=[ "attributes.energy_smearing", 
	"attributes.energy_smearing_units",
	"attributes.absolute_magnetization",
	"attributes.absolute_magnetization_units"],
	output_of="calculation"  )

#essai
qb.append(Group,
        group_of="results",
#        project=["name"],
        filters={"name":{"in": ["tutorial_pbesol", "tutorial_lda", "tutorial_pbe"]   }}      
        )




for i in qb.iterall():
	print i


###  	Ecriture des resultas
out=open(path+"out-query.txt", 'w')

unit_Mag="Bohrmag / cell"
for i in qb.iterall():  #res:
	if(i[3]==None):
		i[3]=0
	if(i[4]==None):
                i[4]=unit_Mag
	out.write("{0}, {1}, {2}, {3}, {4}\n".format(i[0] ,i[1], i[2], i[3], i[4]) )

out.close()

#for queries examplefrom tutorial 

from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData

StructureData=DataFactory("structure")
ParameterData=DataFactory("parameter")

qb=QueryBuilder()
qb.append(RemoteData, tag="remote", project=["*"])
qb.append(Group,group_of="remote",
	filters={"name":{"in": ["tutorial_pbesol", "tutorial_lda", "tutorial_pbe"]   }})

qb.append(ParameterData, project=["attributes.energy_smearing"]


#qb.append(ParameterData, project=["attributes.energy_smearing"],
#	 filters={"id":{"==":1}} )

#qb.append(ParameterData, project=["attributes.energy_smearing"]

qb.all()


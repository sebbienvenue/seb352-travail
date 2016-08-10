from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData
from datetime import datetime
path="/home/aiida/Documents/seb352-travail/essais-tuto/"

StructureData=DataFactory("structure")
ParameterData=DataFactory("parameter")

qb1=QueryBuilder()
qb=QueryBuilder()
#juste
#qb1.append(StructureData, project=["ctime", "id"], filters={"or":[{"ctime":{">=":datetime(2015, 1,1)}},
#	{"id":{"<=":10}} ]})

#StructureData query pour JobCalculation
qb1.append(StructureData, tag="jobs", project=["*"])
qb1.append(JobCalculation, input_of="jobs")


#pour le cas des id less than 10
#qb1.append(StructureData, )


qb.append(Stru#pas pour donner le lien
#for i in qb1.iterall():
#	print i


for i in qb.iterall():
	print i 

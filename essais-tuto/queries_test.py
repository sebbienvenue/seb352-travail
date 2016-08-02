#for queries examplefrom tutorial 

from aiida.orm.querybuilder import QueryBuilder

qb=QueryBuilder()
qb.all()
qb.append(Node)
qb.all()
qb.count()

#enumerate the <pk> for each query key
for node, in qb.iterall():
	print node

#may need this line
StructureData = DataFactory("structure")
qb=QueryBuilder()
qb.append(StructureData)	#met le pk pour chaque structure si on met qb.all()
qb.all()


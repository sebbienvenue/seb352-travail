#for queries examplefrom tutorial 

from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.remote import RemoteData

qb=QueryBuilder()
qb.append(Group, filters={"name":{"in": ["tutorial_pbesol", "tutorial_lda", "tutorial_pbe"]   }})

for i in qb.iterall():
	print i
res=qb.iterdict()

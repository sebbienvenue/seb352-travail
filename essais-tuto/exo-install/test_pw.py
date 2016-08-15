alat=5.4
alat2=alat/2
alat3=alat2/2

from aiida.orm.data.array.kpoints import KpointsData
kpoints = KpointsData()
kpoints_mesh = 2
kpoints.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])

StructureData = DataFactory("structure")
the_cell = [[alat/2,alat/2,0.],[alat/2,0.,alat/2],[0.,alat/2,alat/2]]
structure = StructureData(cell=the_cell)
structure.cell
structure.append_atom(position=(alat/4.,alat/4.,alat/4.),symbols="Si")
structure.sites

from ase.lattice.spacegroup import crystal
ase_structure = crystal('Si', [(0,0,0)], spacegroup=227,
cellpar=[alat, alat, alat, 90, 90, 90],primitive_cell=True)
structure=StructureData(ase=ase_structure)
structure.store()

nameCode='pw@local-seb'

code=Code.get_from_string(nameCode)
calc=code.new_calc()
calc.label="PW test"
calc.description="First calculus on BaTiO3"
calc.set_resources({"num_machines": 1})
calc.set_max_wallclock_seconds(30*60)

calc.use_structure(structure)
calc.use_kpoints(kpoints)

calc.use_pseudos_from_family('SSSP')

#parameters_dict = {'CONTROL': {'calculation':'scf','tstress': True,
#'tprnfor': True,},'SYSTEM': {'ecutwfc': 30.,'ecutrho': 200.,'mickeymouse': 240.,},
#'ELECTRONS': {'conv_thr': 1.e-8,},
#}

#parameters_dict = {'CONTROL': {'calculation':'scf',},'SYSTEM': {'ecutwfc': 30.,'ecutrho': 200.,},
#'ELECTRONS': {'conv_thr': 1.e-6,},
#}

parameters_dict = {'CONTROL': {'calculation':'scf','tstress': True,
'tprnfor': True,},'SYSTEM': {'ecutwfc': 30.,'ecutrho': 200.,},
'ELECTRONS': {'conv_thr': 1.e-6,'electron_maxstep': 300,},
}

ParameterData = DataFactory("parameter")
parameters = ParameterData(dict=parameters_dict)
calc.use_parameters(parameters)

calc.submit_test()

calc.store_all()
print calc.pk

#calc.set_extra("element", "Si")

#calc.submit()
#calc.res.warnings

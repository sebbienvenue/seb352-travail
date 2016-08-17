import numpy as np
from aiida.backends.utils import load_dbenv, is_dbenv_loaded

if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm import CalculationFactory, DataFactory
from aiida.orm import load_node
from aiida.orm.utils import DataFactory
from aiida.workflows2.db_types import Float, Str, NumericType, SimpleData
from aiida.orm.code import Code
from aiida.orm.data.structure import StructureData
from aiida.orm.utils import DataFactory
from aiida.workflows2.run import run
from aiida.workflows2.fragmented_wf import FragmentedWorkfunction, \
    ResultToContext, while_
from aiida.workflows2.wf import wf
from aiida.orm.calculation.job.quantumespresso.pw import PwCalculation
from ase.lattice.spacegroup import crystal
from pressure_convergence import *
#from equation_of_states import *


def fromA_to_au(AA):
    return AA*1.889725989

ParameterData = DataFactory("parameter")
KpointsData = DataFactory("array.kpoints")
StructureData = DataFactory("structure")
PwCalculation = CalculationFactory("quantumespresso.pw")


@wf
def create_diamond_fcc(element):
    """
    Workfunction to create a diamond crystal structure of a given element.
    At the moment only Si and Ge are valid elements (for C we miss the pseudo)
    :param element: The element to create the structure with.
    :return: The structure.
    """
    import numpy as np

    elem_alat = {
        'Si': 5.431,  # Angstrom
        "Ge": 5.658,
    }

    # Validate input element
    symbol = str(element)
    if symbol not in elem_alat.keys():
        raise ValueError("Valid elements are only Si and Ge")

    # Create cel starting from a protopype with alat=1
    alat = elem_alat[symbol]
    the_cell = np.array([[0., 0.5, 0.5],
                         [0.5, 0., 0.5],
                         [0.5, 0.5, 0.]]) * alat

    # Create a structure data object
    StructureData = DataFactory("structure")
    structure = StructureData(cell=the_cell)
    structure.append_atom(position=(0., 0., 0.), symbols=str(element))
    structure.append_atom(position=(0.25 * alat, 0.25 * alat, 0.25 * alat), symbols=str(element))
    return structure


@wf
def create_bcc(element):
    """
    Workfunction to create a bcc structure of a given element.
    At the moment only Fe are valid elements
    :param element: The element to create the structure with.
    :return: The structure.
    """
    import numpy as np

    elem_alat = {
        'Fe': fromA_to_au(2.86),  # Angstrom
         }

    # Validate input element
    symbol = str(element)
    if symbol not in elem_alat.keys():
        raise ValueError("Valid elements are only Si and Ge")

    # Create cel starting from a protopype with alat=1
    alat = elem_alat[symbol]
    the_cell = np.array([[1., 0., 0.],
                         [0., 1., 0.],
                         [0., 0., 1.]]) * alat

    # Create a structure data object
    StructureData = DataFactory("structure")
    structure = StructureData(cell=the_cell)
    structure.append_atom(position=(0., 0., 0.), symbols=str(element))
    structure.append_atom(position=(0.5 * alat, 0.5 * alat, 0.5 * alat), symbols=str(element))
    return structure

@wf
def rescale(structure, scale):
    """
    Workfunction to rescale a structure

    :param structure: An AiiDA structure to rescale
    :param scale: The scale factor
    :return: The rescaled structure
    """
    the_ase = structure.get_ase()
    new_ase = the_ase.copy()
    new_ase.set_cell(the_ase.get_cell() * float(scale), scale_atoms=True)
    new_structure = DataFactory('structure')(ase=new_ase)
    return new_structure


@wf
def create_rescaled_bcc(element,  scale):
	"""Workfunction and rescale the structure given
	"""
	s0= create_bcc(element)
	return rescale(s0, scale)

def get_pseudos(structure, family_name):
    """
    Set the pseudo to use for all atomic kinds, picking pseudos from the
    family with name family_name.

    :note: The structure must already be set.

    :param family_name: the name of the group containing the pseudos
    """
    from collections import defaultdict
    from aiida.orm.data.upf import get_pseudos_from_structure

    # A dict {kind_name: pseudo_object}
    kind_pseudo_dict = get_pseudos_from_structure(structure, family_name)

    # We have to group the species by pseudo, I use the pseudo PK
    # pseudo_dict will just map PK->pseudo_object
    pseudo_dict = {}
    # Will contain a list of all species of the pseudo with given PK
    pseudo_species = defaultdict(list)

    for kindname, pseudo in kind_pseudo_dict.iteritems():
        pseudo_dict[pseudo.pk] = pseudo
        pseudo_species[pseudo.pk].append(kindname)

    pseudos = {}
    for pseudo_pk in pseudo_dict:
        pseudo = pseudo_dict[pseudo_pk]
        kinds = pseudo_species[pseudo_pk]
        for kind in kinds:
            pseudos[kind] = pseudo

    return pseudos



def generate_scf_input_params(structure, codename, pseudo_family):
    # The inputs
    inputs = PwCalculation.process().get_inputs_template()

    # The structure
    inputs.structure = structure

    inputs.code = Code.get_from_string(codename)
    # calc.label = "PW test"
    # calc.description = "My first AiiDA calculation of Silicon with Quantum ESPRESSO"
    inputs._options.resources = {"num_machines": 1}
    inputs._options.max_wallclock_seconds = 30 * 60

    # Kpoints
    KpointsData = DataFactory("array.kpoints")
    kpoints = KpointsData()
    kpoints_mesh = 8 #kmesh
    kpoints.set_kpoints_mesh([kpoints_mesh, kpoints_mesh, kpoints_mesh])
    inputs.kpoints = kpoints

    # Calculation parameters
    parameters_dict = {
        "CONTROL": {"calculation": "scf",
                    "tstress": True,  #  Important that this stays to get stress
                    "tprnfor": True,
                    "restart_mode":'from_scratch',},
        "SYSTEM": {"ecutwfc": 55.,
                   "ecutrho": 240.,
                   "starting_magnetization(1)":0.7,
                   "occupations":'smearing',
                   "degauss":0.03,
                   "smearing":'m-v',
                   "nspin":2,},
        "ELECTRONS": {"mixing_beta":0.7,
                      "conv_thr": 1.e-8,}
    }
    ParameterData = DataFactory("parameter")
    inputs.parameters = ParameterData(dict=parameters_dict)

    # Pseudopotentials
    inputs.pseudo = get_pseudos(structure, pseudo_family)

    return inputs


def get_info(calc_results):
    return (calc_results['output_parameters'].dict.volume,
            calc_results['output_parameters'].dict.energy,
            calc_results['output_parameters'].dict.energy_units)


scale_facs = (0.96, 0.98, 1.0, 1.02, 1.04)
labels = ["c1", "c2", "c3", "c4", "c5"]


class EquationOfStates(FragmentedWorkfunction):
    @classmethod
    def _define(cls, spec):
        spec.input("element", valid_type=SimpleData)
        spec.input("code", valid_type=SimpleData)
        spec.input("pseudo_family", valid_type=SimpleData)
        spec.outline(
            cls.run_pw,
            cls.return_results,
        )
        spec.dynamic_output()

    def run_pw(self, ctx):

        PwProcess = PwCalculation.process()
        ctx.s0 = create_bcc(Str(self.inputs.element))

        ctx.eos_names = []
        calcs = {}

        for label, factor in zip(labels, scale_facs):
            s = rescale(ctx.s0, Float(factor))
            inputs = generate_scf_input_params(
                s, str(self.inputs.code), str(self.inputs.pseudo_family))
            print "Running a scf for {} with scale factor {}".format(
                self.inputs.element, factor)

            # Launch the code
            future = self.submit(PwProcess, inputs)
            # Store the future
            calcs[label] = future

        # Ask the workflow to continue when the results are ready and store them
        # in the context
        return ResultToContext(**calcs)

    def return_results(self, ctx):
        eos = []
        for label in labels:
            eos.append(get_info(ctx[label]))

        # Return information to plot the EOS
        ParameterData = DataFactory("parameter")
        retdict = {
            'initial_structure': ctx.s0,
            'result': ParameterData(dict={'eos_data': eos})
        }
        for k, v in retdict.iteritems():
            self.out(k, v)























if __name__ == "__main__":
    import argparse
    path="/home/bienvenue/Documents/AiiDa/Atomistic-Fe/"
    file_log=open(path+"LOG", 'w')
    #print path
    # codename = 'pw@localhost'
    # code = Code.get_from_string(codename)
    parser = argparse.ArgumentParser(description='Energy calculation example for BCC iron .')

    structure = create_bcc(element=Str('Fe'))

    print "Initial structure:", structure
    file_log.write("Initial structure:"+ str(structure))
    args = parser.parse_args()
    #wf_results = run(PressureConvergence, structure=structure, code=Str('pw@localhost'), pseudo_family=Str('atomistic'),
    #                 volume_tolerance=Float(0.1))

    run(EquationOfStates, element=Str("Fe"),
        code=Str('pw-deneb@sebbienvenue'), pseudo_family=Str('atomistic'))

    print "Workflow results:"
    print wf_results
    file_log.write("Workflow results:"+str(wf_results))


    file_log.close()










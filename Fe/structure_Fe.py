from aiida.backends.utils import load_dbenv, is_dbenv_loaded

if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm import load_node
from aiida.orm.utils import DataFactory
from aiida.workflows2.db_types import Float, Str, NumericType, SimpleData
from aiida.orm.code import Code
from aiida.orm.data.structure import StructureData
from aiida.workflows2.run import run
from aiida.workflows2.fragmented_wf import FragmentedWorkfunction, \
    ResultToContext, while_
from aiida.workflows2.wf import wf
#from common_wf import generate_scf_input_params
#from create_rescale import rescale, create_diamond_fcc
from aiida.orm.calculation.job.quantumespresso.pw import PwCalculation
from ase.lattice.spacegroup import crystal

def fromA_to_au(AA):
    return AA*1.889725989

ParameterData = DataFactory("parameter")
KpointsData = DataFactory("array.kpoints")
StructureData = DataFactory("structure")

alatBCC= fromA_to_au(2.86)

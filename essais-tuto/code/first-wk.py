from create_rescale import create_diamond_fcc, rescale
from aiida.workflows2.db_types import Float, Str
from common_wf import get_pseudos, generate_scf_input_params
#futher imports
from aiida.workflows2.run import run 
from aiida.workflows2.wf import wf
from aiida.orm.calculation.job.quantumespresso.pw import PwCalculation
from aiida.workflows2.defaults import registry

def run_eos(codename='pw-5.1@localhost', pseudo_family='GBRV_lda', element="Si"):
	return run_eos_wf(Str(codename), Str(pseudo_family), Str(element))

@wf
def run_eos_wf(codename, pseudo_family, element):
	print "Workfunction node pk: {}".format(registry.current_calc_node)
	#Instantiate a JobCalc process and create basic structure
	JobCalc = PwCalculation.process()
	s0 = create_diamond_fcc(Str(element))
	eos=[]
	scale_facs = (0.98, 0.99, 1.0, 1.02, 1.04)
	for factor in  scale_facs:
		s = rescale(s0,Float(factor))
		inputs = generate_scf_input_params(
			s, str(codename), str(pseudo_family))
		print "Running a scf for {} with scale factor {}".format(
			element, factor)
	calc_results = run(JobCalc,**inputs)
	eos.append(get_info(calc_results))
	#Return information to plot the EOS
	ParameterData = DataFactory("parameter")
	return {'initial_structure': s0,'result': ParameterData(dict={'eos_data': eos})}

def get_info(calc_results):
	return (calc_results['output_parameters'].dict.volume,
		calc_results['output_parameters'].dict.energy,
		calc_results['output_parameters'].dict.energy_units)

###execution 

run_eos()

#to keep track of the pk registration
#from aiida.workflows2.defaults import registry
#print registry.current_calc_node

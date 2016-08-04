from common_wf import get_pseudos, generate_scf_input_params
from create_rescale import create_diamond_fcc, rescale
from aiida.workflows2.run import run
from aiida.workflows2.fragmented_wf import FragmentedWorkfunction, \
    ResultToContext
from aiida.workflows2.db_types import Float, Str, NumericType, SimpleData
from aiida.orm.code import Code

# Set up the factories
ParameterData = DataFactory("parameter")
KpointsData = DataFactory("array.kpoints")
StructureData = DataFactory("structure")
PwCalculation = CalculationFactory("quantumespresso.pw")


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
        ctx.s0 = create_diamond_fcc(Str(self.inputs.element))

        ctx.eos_names = []
        calcs = {}
 
        for label, factor in zip(labels, scale_facs):
            s = rescale(ctx.s0,Float(factor))
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

        #Return information to plot the EOS
        ParameterData = DataFactory("parameter")
        retdict = {
                'initial_structure': ctx.s0,
                'result': ParameterData(dict={'eos_data': eos})
	       }
        for k, v in retdict.iteritems():
            self.out(k, v)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Equation of states example.')
    parser.add_argument('--pseudo', type=str, dest='pseudo',
                        help='The pseudopotential family', required=True)
    parser.add_argument('--code', type=str, dest='code',
                        help='The codename to use', required=True)

    args = parser.parse_args()
  
    # Get the structure from the calculation
    run(EquationOfStates, element=Str("Si"),
        code=Str(args.code), pseudo_family=Str(args.pseudo))

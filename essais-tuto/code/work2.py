from aiida.workflows2.db_types import Float, Str
from ../create_rescale import create_diamond_fcc, rescale

s0=create_diamond_fcc(Str("Si"))
rescaled_structures= [rescale(s0, Float(factor)) for factor in (0.98 , 0.99, 1.0, 1.1, 1.2)]

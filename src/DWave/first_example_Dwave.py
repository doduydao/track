# Calling library:
from dimod import ConstrainedQuadraticModel, Integer, Binary, Real

# Define variable:
# - Integer
# - Binary
# - Real
i = Integer('i', upper_bound=4)
j = Integer('j', upper_bound=4)

# Define type of model:
# - Binary Quadratic Models: are unconstrained and have binary variables
# - Constrained Quadratic Models: can be constrained and have real, integer and binary variables.
# - Discrete Quadratic Models: are unconstrained and have discrete variables.
cqm = ConstrainedQuadraticModel()

# Define objective function
cqm.set_objective(-i * j)

# Adding constraints
cqm.add_constraint(2 * i + 2 * j <= 8, "Max perimeter")

print(cqm)

# Library for sampling:
from dwave.system import LeapHybridCQMSampler

# Define sampler:
# - Hybrid Solvers such as Leap’s hybrid_binary_quadratic_model_version<x> solver or,
# for discrete quadratic models (DQM), hybrid_discrete_quadratic_model_version<x>
# - Classical Solvers such as ExactSolver for exact solutions to small problems
# - Quantum Solvers such as the Advantage system.
sampler = LeapHybridCQMSampler()

# Sampling
answer = sampler.sample_cqm(cqm)
print(answer)

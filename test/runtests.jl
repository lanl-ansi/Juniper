using Logging

# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test
using JuMP

using Ipopt
using PowerModels

using MINLPBnB
using MINLPBnB.BnBTree

# include("load_mod.jl")

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3

# model_from_mod("data/batch0812_nc.mod")

minlpbnb_strong = MINLPBnB.MINLPBnBSolver(IpoptSolver(print_level=0);
                                    log_levels=[:Table],
                                    branch_strategy=:StrongPseudoCost,
                                    strong_branching_nvars = 3
                                )
minlpbnb_mosti = MINLPBnB.MINLPBnBSolver(IpoptSolver(print_level=0);
                                log_levels=[:Table],
                                branch_strategy=:MostInfeasible,
                            )  

minlpbnb_pseudo = MINLPBnB.MINLPBnBSolver(IpoptSolver(print_level=0);
                            log_levels=[:Table],
                            branch_strategy=:PseudoCost,
                        )                               

start = time()
include("basic.jl")
include("power_models_acp.jl")
include("power_models_socwr.jl")
println("Time for all tests: ", time()-start)
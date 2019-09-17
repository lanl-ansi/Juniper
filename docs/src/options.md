## General

The most basic configuration of Juniper is:

```
optimizer = Juniper.Optimizer
params = Dict{Symbol,Any}()
params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
```

and then creating the model with:
```
m = Model(with_optimizer(optimizer, params))
```

The argument `nl_solver` defines the solver for the relaxation here `IpoptSolver`. [Ipopt](https://projects.coin-or.org/Ipopt) is used for all our test cases. The Ipopt julia package is described [here](https://github.com/JuliaOpt/Ipopt.jl). The solver itself can have parameters i.e `print_level=0`.

JuMP supports a lot of different NLP solvers (open source as well as commercial). A list of some NLP solvers is mentioned [here](http://www.juliaopt.org/JuMP.jl/0.18/installation.html#getting-solvers)

You can add options by adding another key, value pair to the `params` dictionary:

```
params[:branch_strategy] = :StrongPseudoCost
```

In this example the strategy used for branching is defined.

In the following the options are explained. The type for the option is given after `::` and the default value in `[]`.

**Attention:**
The default values might change in the future after several tests were executed to determine the best overall options. 

## Tolerance

### atol::Float64 [1e-6]

This tolerance is used to check whether a value is integer or not and is used in the feasibility pump.
More information about the feasibility pump can be found [here](#Feasibility-Pump-1).

## Branching

### branch_strategy::Symbol [:StrongPseudoCost]

Possible values:

* `:MostInfeasible`
    * Branch on variables closest to 0.5
* `:PseudoCost`
    * Use `:MostInfeasible` first and then [Pseudo Cost Branching](https://en.wikipedia.org/wiki/Branch_and_cut#Branching_Strategies).
* `:StrongPseudoCost`
    * Use [Strong Branching](https://en.wikipedia.org/wiki/Branch_and_cut#Branching_Strategies) first and then `:PseudoCost`
    * More options for strong branching are described [here](#Options-for-strong-branching-1)
* `:Reliability`
    * Use [Reliability Branching](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf) in a slightly different version.

### traverse_strategy::Symbol [:BFS]

Determines which node is used for the next branching.

Possible values:

* `:BFS`
    * Best-first-search: always take the node with the best bound
* `:DFS`
    * Depth-first-search: always take the node with the highest level
    * Might find a feasible solution faster
* `:DBFS`
    * Use of DFS first until the first feasible solution is found then switch to BFS

## Objective Cuts
### incumbent_constr::Bool [false]

Add a constraint `objval >=/<= incumbent`. 

### obj_epsilon::0 [Float64]

Add a constraint of the following form at the root node if not `0`.

If minimizing:

$\text{obj } \leq (1+\epsilon)\text{LB}$

If maximizing:

$\text{obj } \geq (1-\epsilon)\text{UB}$

## Options for strong branching

### strong\_branching\_perc::Float64 [100]

Defines the percentage of variables to consider for strong branching. 
If set to 25 it means that strong branching is performed on 25% of all discrete variables.
Variables which are discrete in the relaxation aren't considered again but count to the number of 
all discrete variables.
If the number of variables is smaller than `2` it is fixed at `2` as strong branching doesn't make sense for one variable. **Attention:** `strong_branching_time_limit` might terminate strong branching earlier.

### strong\_branching\_nsteps::Int64 [1]

Defines the number of steps in which strong branching is used. `:PseudoCost` will be used for later steps.

### strong\_branching\_time\_limit::Float64 [100]s

For big problems with either a lot of variables or a long relaxation time it turned out to be reasonable
to reduce the number of strong branching variables. Strong branching will terminate if a variable index was chosen and strong branching took longer than the time limit allowed.

If you don't want to use this time limit you can set it to `Inf`.

### strong_restart::Bool [true]

If a child, while running strong branching, is infeasible this holds for the whole node. Therefore we can tighten the bounds and rerun the strong branch part. (This might occur more then once)
This option is also used in reliablity branching.

## Options for reliablity branching

The implemented version of reliablity branching uses the gain score as in pseudo cost branching 
and if some branching variables aren't reliable in a sense that strong branching wasn't performed 
at least `reliablility_branching_threshold` times then strong branching is performed on those.
Afterwards it will be branched on the variable with the highest gain score.

### reliablility\_branching\_perc::Float64 [25]

Defines the percentage of variables to consider for the strong branching part of reliablity branching.
If the number of variables is smaller than `2` it is fixed at `2` as strong branching doesn't make sense for one variable. 

### reliablility\_branching\_threshold::Int64 [5]

Defines whether strong branching is used to improve the reliability of the gain score.
If a variable was used less than `reliablility_branching_threshold` times for strong branching then strong branching is performed on some of those candidates. The amount of candidates used is calculated by `reliablility_branching_perc`.

### strong_restart::Bool [true]

This option is explained in strong branching but is also used for reliability branching.

## Options gain computation

The objective gain is used for pseudo cost and therefore also for strong pseudo cost branching.
$\mu$ is a parameter inside the computation of a so called score.
The version we use for the score function is described [here](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf).

### gain_mu::Float64 [0.167]

The parameter is used and a bit described in this [paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf) in (3).

## Parallel

Juniper can be run in parallel to speed up the algorithm.
You have to start julia with `julia -p P` where `P` is the number of processors available or at least the number of processors you want to use.

Then you have to specify the number of processor as an option.

### processors::Int64 [1]

The number of processors used for the branch and bound part. **Attention:** Even if you start julia using
`julia -p P` you still have to define the number of processors using this option.

### two\_processors\_per\_node::Bool [true]

In the branch and bound part of Juniper each branching node has two children these will be computed in parallel using two processors when set to `true` and Juniper is running in parallel. If `false` the two processors will be used instead to compute two open nodes in parallel.

## Feasibility Pump

Juniper has the option to find a feasible solution before the branch and bound part starts. The following options describe how to use the feasibility pump.

### feasibility_pump::Bool [Auto]

Determines whether or not the feasibility pump should be used to get a feasible solution.
The default is `true` if a mip solver is specified and `false` otherwise.
**Attention**: If set to `true` you need to also set the `mip_solver` option.

### mip_solver::JuMP.OptimizerFactory [nothing]

This has to be set to a mip solver if the feasibility pump should be used.
A list of some MIP solvers is mentioned [here](http://www.juliaopt.org/JuMP.jl/0.18/installation.html#getting-solvers).

If you want to use [Cbc](https://projects.coin-or.org/Cbc)
you would need to use

```
using Cbc
```
and set the option with 
```
params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)
```

### feasibility\_pump\_time\_limit::Int64 [60]s

The time limit of the feasibility pump in seconds. After that time limit the branch and bound part starts whether a feasible solution was found or not.

### feasibility\_pump\_tolerance\_counter::Int64 [5]

In the feasibility pump the objective is to reduce the difference between the mip and the nlp solution.
If the default tolerance [`atol`](#atol::Float64-[1e-6]-1) can't be reached for `feasibility_pump_tolerance_counter` consecutive times
but the `tolerance*10` can be reached, a warning will be thrown and the tolerance is increased.

If there is no warning like `Real objective wasn't solved to optimality` afterwards there is no need to worry at all. Normally in the end of the feasibility pump the real objective is used to improve the 
objective value. In the case that the nlp was solved (maybe only locally like ipopt in non-convex cases) the warning before can be ignored as it is given that the solution has no rounding issues. 
In the other case a warning like `Real objective wasn't solved to optimality` is thrown.
This means that the objective might be not the best possible given the mip solution and if a warning for the tolerance change was thrown there might be rounding issues. 
You can set this value to a huge number (more than 100 should be enough) if you don't want to use this option.

### tabu\_list\_length::Int64 [30]

During the run of the feasibility pump it might happen that the alternating solve steps get into a cycle.
By using a tabu list cycles can be avoided. The length determines the length of the cycle which will be avoided. If a cycle is encountered which is longer the feasibility pump terminates.

### num\_resolve\_nlp\_feasibility\_pump::Int64 [1]

If the NLP is infeasible during the feasibility pump it can be restarted with a random starting point for the NL solver. This will be done as long as it is infeasible or `num_resolve_nlp_feasibility_pump` is reached.

## User Limits

You can stop the solver before the optimal solution is found using various options explained in this section.
This is reasonable if the problem is to big to solve to (local) optimality fast.
If the solver stops because of one of the following options an appropriate status like `MOI.TIME_LIMIT` is returned.

### time_limit::Float64 [Inf]

The maximum time in seconds the solver should run. 

If this limit is reached the status will be `MOI.TIME_LIMIT`.

**Note:** The solver will check after each branching step therefore it isn't a strict limit and depends on the duration of a relaxation.

### mip_gap::Float64 [0.0001]

Stops the solver if the gap is smaller than `mip_gap`. The default is `0.01%`.

If this limit is reached the status will be `MOI.OBJECTIVE_LIMIT`.

### best\_obj\_stop::Float [NaN]

If an incumbent is found which is better than `best_obj_stop` the incumbent is returned. A warning gets thrown if `best_obj_stop` can't be reached.

If this limit is reached the status will be `MOI.OBJECTIVE_LIMIT`.

### solution_limit::Int [0]

The solver stops if the requested amount of feasible solutions is found.
If `0` the option gets ignored.

If this limit is reached the status will be `MOI.SOLUTION_LIMIT`.

## Resolve

Sometimes the non linear solver doesn't find a feasible solution in the first run.

### num\_resolve\_root\_relaxation::Int [3]
This especially bad if this happens for the root relaxation. If there is no optimal/local optimal
solution in the root relaxation you can use this option to resolve a couple of time until a solution is found or the number of resolves exceeded this value.

## Extra options

### allow\_almost\_solved\_integral::Bool [true]
The non linear solver might return the status `ALMOST_LOCALLY_SOLVED` which means: 
"The algorithm converged to a stationary point, local
  optimal solution, or could not find directions for improvement within relaxed tolerances."
Inside Juniper this is mostly considered as `LOCALLY_SOLVED` (see next option) but you can use this option to restart the search once if solution is integral but only `ALMOST_LOCALLY_SOLVED` to maybe find a `LOCALLY_SOLVED` solution. 

### allow\_almost\_solved::Bool [true]
See above option for an explanation of `ALMOST` solved status codes. You can completely disable allowing such status codes with this option.
Setting it to true means that all `ALMOST` solved status codes are considered as Infeasible/Numerical error throughout the tree search.

## Logging

### log_levels::Vector{Symbol} [[:Table,:Info,:Options]]

You can change the option `log_levels` to define what kind of logs you want to see.

The output for `[:Table,:Info]` looks something like this:

![default-logging](https://user-images.githubusercontent.com/4931746/32625934-07b7db3c-c58e-11e7-922d-18a0a8776437.png)

**:Table**

- #ONodes
    - The number of open nodes
- CLevel
    - The current node is at level ... of the tree
- Incumbent
    - Best integral solution found
- Best Bound
    - The best bound of the open nodes
- Gap 
    - The gap between `Incumbent` and `Best Bound`
- Time
    - The time spend since the beginning of branch and bound
    - Doesn't count time before branch and bound starts (i.e. feasibility pump or root relaxation)
- GainGap
    - The difference in percentage between a guessed gain and the actual gain.
    - Used if `branch_strategy = PseudoCost` or after strong branching / reliability branching.

**:Options**

includes something like this before the info is printed:

```
time_limit               : 10.0
strong_branching_nsteps  : 5
```

Possible symbols which can be added to the vector are:

- :Timing
    - Provides some more timing information
- :AllOptions
    - prints all options 
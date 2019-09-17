# Juniper.jl Changelog

### Unreleased
- Solving the two children of a node in parallel
  - Corresponding options `processors` and `two_processors_per_node` 
- Changed option `strong_branching_approx_time_limit` to `strong_branching_time_limit`
  - Default is still 100s

### 0.5.0
- Upgraded to MOI 0.9.1 and JuMP v0.20.0
  - Implementation of `TimeLimitSec`, `Silent` & `RawParameter`

### v0.4.3
- Bugfix: No bounds on binary variables in Feasibility pump [#143](https://github.com/lanl-ansi/Juniper.jl/issues/143)
- Bugfix: Feasibility pump if no objective in problem [#145](https://github.com/lanl-ansi/Juniper.jl/issues/145)

### v0.4.2
- Ability to not accept `almost` MOI solver status codes in tree search
- Better handling of MOI solver status codes internally
- Bugfix: Don't call objective value if sub problem not solved [#130](https://github.com/lanl-ansi/Juniper.jl/issues/130)
- Return ALMOST_LOCALLY_SOLVED if corresponding relaxation is only almost solved

### v0.4.1
- Support for user defined functions see issue [#118](https://github.com/lanl-ansi/Juniper.jl/issues/118)

### v0.4.0
- MPB -> MOI
- JuMP 0.19
- Dropped support for Julia versions < 1.0

### v0.3.0
- removes support for Julia < v1

### v0.2.6
- bugfix in init_strong_restart
- bugfix in mip_gap if objval=0
- bugfix if gap was 0 in table printing 
- support for primal start values

### v0.2.5
- add linear constraints using @constraint in root model
- in BFS mode use the highest depth if several nodes have the same objective

### v0.2.4
- Remove Julia version upper bound

### v0.2.3
- Again support for `mu_init` in Ipopt

### v0.2.2
- Add support for Julia v1.0

### v0.2.1
- Add support for Julia v0.6 and v0.7

### v0.2.0
- Bugfix: Reset of `mu_init` in Ipopt options to have the default `mu_init` if `solve` is called again
- Bugfix: Break on time limit in relaxation, fpump and strong branching
- Bugfix: Infeasible in Reliability branching
- Strong branching: 
    - Change bounds even if no restart
    - branch on best variable with two children

### v0.1.1
- Freemodel for commerical nlp solvers with license restrictions
- More convenient parallel options 
    - `processors = 2` now uses 2 processors for solving nodes and one thread for supervision

### v0.1.0
- Traverse options (BFS,DFS,DBFS)
- Branch options (Strong, Reliable, Pseudo, MostInfeasible)
- Parallel solving of nodes
- Basic feasibility pump

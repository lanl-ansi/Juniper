export JuniperSolver

#=
A solver for MINLP problems using a NLP solver and Branch and Bound
=#

function get_default_options()
    log_levels                          = [:Options,:Table,:Info]
    atol                                = 1e-6
    num_resolve_root_relaxation         = 3
    branch_strategy                     = :StrongPseudoCost
    gain_mu                             = 0.167 # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf
    # Strong branching  
    strong_branching_perc               = 100
    strong_branching_nsteps             = 1
    strong_branching_approx_time_limit  = 100
    strong_restart                      = true
    # Reliability branching     
    reliability_branching_threshold     = 5 # reliability param
    reliability_branching_perc          = 25
    # Obj cuts  
    incumbent_constr                    = false
    obj_epsilon                         = 0
    # :UserLimit    
    time_limit                          = Inf  
    mip_gap                             = 1e-4
    best_obj_stop                       = NaN
    solution_limit                      = 0
    all_solutions                       = false
    list_of_solutions                   = false
    # Parallel  
    processors                          = 1
    # Traversing    
    traverse_strategy                   = :BFS
    # Feasibility Pump  
    feasibility_pump                    = true # changes to false if mip_solver not provided
    feasibility_pump_time_limit         = 60
    feasibility_pump_tolerance_counter  = 5
    tabu_list_length                    = 30
    num_resolve_nlp_feasibility_pump    = 1
    mip_solver                          = nothing

    # Only for testing
    force_parallel                      = false
    debug                               = false
    debug_write                         = false
    debug_file_path                     = "debug.json"

    fixed_gain_mu                       = false

    return SolverOptions(log_levels,atol,num_resolve_root_relaxation,branch_strategy,gain_mu,
        strong_branching_perc,strong_branching_nsteps,strong_branching_approx_time_limit,strong_restart,
        reliability_branching_threshold,reliability_branching_perc,
        incumbent_constr,obj_epsilon,time_limit,mip_gap,best_obj_stop,solution_limit,all_solutions,
        list_of_solutions,processors,traverse_strategy,
        feasibility_pump,feasibility_pump_time_limit,feasibility_pump_tolerance_counter,
        tabu_list_length,num_resolve_nlp_feasibility_pump,
        mip_solver, force_parallel, debug, debug_write, debug_file_path, fixed_gain_mu)
end

function combine_options(options)
    branch_strategies = Dict{Symbol,Bool}()
    for strat in [:StrongPseudoCost,:PseudoCost,:Reliability,:MostInfeasible]
        branch_strategies[strat] = true
    end

    traverse_strategies = Dict{Symbol,Bool}()
    traverse_strategies[:BFS] = true
    traverse_strategies[:DFS] = true
    traverse_strategies[:DBFS] = true

    options_dict = Dict{Symbol,Any}()
    for kv in options
        if !in(kv[1], fieldnames(SolverOptions))
            @warn "Option "*string(kv[1])*" is not available"
        end
        options_dict[kv[1]] = kv[2]
    end
    if haskey(options_dict, :log_levels)
        if length(options_dict[:log_levels]) == 0
            options_dict[:log_levels] = Symbol[]
        end
    end
    defaults = get_default_options()
    if defaults.feasibility_pump == true && (!haskey(options_dict, :mip_solver) || options_dict[:mip_solver] == nothing)
        defaults.feasibility_pump = false
    end

    # if gain mu is specified we shouldn't change it later
    if haskey(options_dict,:gain_mu)
        options_dict[:fixed_gain_mu] = true
    end

    for fname in fieldnames(SolverOptions)
        if haskey(options_dict, fname)
            # check that mip_solver is defined if feasibile pump should be used
            if fname == :feasibility_pump && options_dict[:feasibility_pump] == true
                if !haskey(options_dict, :mip_solver) || options_dict[:mip_solver] == nothing
                    @warn "The feasibility pump can only be used if a mip solver is defined."
                    options_dict[:feasibility_pump] = false
                end
            end

            # check branch strategy
            if fname == :branch_strategy 
                if !haskey(branch_strategies, options_dict[fname])
                    error("Branching strategy "*string(options_dict[fname])* " is not supported")
                end
            end

             # check traverse strategy
             if fname == :traverse_strategy 
                if !haskey(traverse_strategies, options_dict[fname])
                    error("Traverse strategy "*string(options_dict[fname])* " is not supported")
                end
            end

            if fieldtype(SolverOptions, fname) != typeof(options_dict[fname])
                options_dict[fname] = convert(fieldtype(SolverOptions,fname), options_dict[fname])
            end
            setfield!(defaults, fname, options_dict[fname])
        end
    end
    return defaults
end

function JuniperSolver(nl_solver::MathProgBase.AbstractMathProgSolver;options...)
    options_obj = combine_options(options)
    return JuniperSolverObj(nl_solver,options_obj)
end


# Create Pavito conic model
MathProgBase.ConicModel(s::JuniperSolverObj) = MathProgBase.ConicModel(ConicNonlinearBridge.ConicNLPWrapper(nlp_solver=s))

# Create Pavito LinearQuadratic model
MathProgBase.LinearQuadraticModel(s::JuniperSolverObj) = MathProgBase.NonlinearToLPQPBridge(MathProgBase.NonlinearModel(s))

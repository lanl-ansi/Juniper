mutable struct RegisteredFunction
    s         :: Symbol
    dimension :: Integer
    f         :: Function
    gradf     :: Union{Nothing, Function}
    grad2f    :: Union{Nothing, Function}
    autodiff  :: Bool
end


# Options for the solver (more details like defaults in solver.jl)
mutable struct SolverOptions
    nl_solver                           :: Union{Nothing, JuMP.OptimizerFactory} # needs to be set
    log_levels                          :: Vector{Symbol}
    silent                              :: Bool
    atol                                :: Float64
    num_resolve_root_relaxation         :: Int64
    branch_strategy                     :: Symbol
    gain_mu                             :: Float64
    strong_branching_perc               :: Float64
    strong_branching_nsteps             :: Int64
    strong_branching_approx_time_limit  :: Float64
    strong_restart                      :: Bool
    reliability_branching_threshold     :: Int64
    reliability_branching_perc          :: Float64
    incumbent_constr                    :: Bool
    obj_epsilon                         :: Float64
    time_limit                          :: Float64
    mip_gap                             :: Float64
    best_obj_stop                       :: Float64
    solution_limit                      :: Int64
    all_solutions                       :: Bool
    list_of_solutions                   :: Bool
    processors                          :: Int64
    traverse_strategy                   :: Symbol
    feasibility_pump                    :: Bool
    feasibility_pump_time_limit         :: Float64
    feasibility_pump_tolerance_counter  :: Int64
    tabu_list_length                    :: Int64
    num_resolve_nlp_feasibility_pump    :: Int64
    mip_solver                          :: Union{Nothing, JuMP.OptimizerFactory}
    allow_almost_solved                 :: Bool  
    allow_almost_solved_integral        :: Bool  
    registered_functions                :: Union{Nothing, Vector{RegisteredFunction}}
    
    # only for testing
    force_parallel                      :: Bool
    debug                               :: Bool
    debug_write                         :: Bool
    debug_file_path                     :: String

    fixed_gain_mu                       :: Bool
end

mutable struct SolutionObj
    solution    :: Vector{Float64}
    objval      :: Float64
end

# Juniper MOI struct 

mutable struct JuniperProblem 
    nl_solver           :: JuMP.OptimizerFactory
    nl_solver_options   :: Vector{Tuple}
   
    model               :: JuMP.Model

    relaxation_status   :: MOI.TerminationStatusCode
    relaxation_objval   :: Float64
    relaxation_solution :: Vector{Float64}

    status              :: MOI.TerminationStatusCode
    objval              :: Float64
    best_bound          :: Float64

    x                   :: Vector{JuMP.VariableRef}
    primal_start        :: Vector{Real}
    num_constr          :: Int64
    num_nl_constr       :: Int64
    num_q_constr        :: Int64
    num_l_constr        :: Int64
    num_var             :: Int64
    l_var               :: Vector{Float64}
    u_var               :: Vector{Float64}

    has_nl_objective    :: Bool
    nlp_evaluator       :: MOI.AbstractNLPEvaluator

    objective           :: Union{SVF, SAF, SQF, Nothing}

    disc2var_idx        :: Vector{Int64}
    var2disc_idx        :: Vector{Int64}

    var_type            :: Vector{Symbol}
    obj_sense           :: Symbol
    num_disc_var        :: Int64

    solution            :: Vector{Float64}

    soltime             :: Float64
    options             :: SolverOptions
    solutions           :: Vector{SolutionObj}
    nsolutions          :: Int64

    mip_solver          :: JuMP.OptimizerFactory
    mip_solver_options  :: Vector{Tuple}

    relaxation_time     :: Float64
    start_time          :: Float64

    # Info  
    nintvars            :: Int64
    nbinvars            :: Int64
    nnodes              :: Int64
    ncuts               :: Int64
    nbranches           :: Int64
    nlevels             :: Int64

    fpump_info          :: Dict{Symbol,Float64}

    # debug 
    debugDict           :: Dict{Symbol,Any}

    JuniperProblem() = new()
end 

###########################################################################
########################## FPump ##########################################
###########################################################################
mutable struct TabuList
    sols      :: Vector{Vector{Float64}}
    length    :: Int64
    pointer   :: Int64

    TabuList() = new()
end

############################################################################
######################## Tree structure ####################################
############################################################################
mutable struct BnBNode
    idx                 :: Int64
    level               :: Int64
    l_var               :: Vector{Float64}
    u_var               :: Vector{Float64}
    solution            :: Vector{Float64}
    var_idx             :: Int64
    state               :: Symbol
    relaxation_state    :: MOI.TerminationStatusCode
    best_bound          :: Float64
    path                :: Vector{String}   # list of parent hashes
    hash                :: String
end

function Base.copy(b::BnBNode) 
    return BnBNode(b.idx, b.level, copy(b.l_var), copy(b.u_var), copy(b.solution), b.var_idx, b.state, b.relaxation_state, b.best_bound, copy(b.path), b.hash)
end

mutable struct Incumbent
    objval      :: Float64
    solution    :: Vector{Float64}
    only_almost :: Bool
end

mutable struct GainObj
    minus           :: Vector{Float64} # gain of objective per variable on left node
    plus            :: Vector{Float64} # gain of objective per variable on right node
    minus_counter   :: Vector{Int64} # obj_gain_m / obj_gain_mc => average gain on left node
    plus_counter    :: Vector{Int64} # obj_gain_p / obj_gain_pc => average gain on right node
    # counter of how often one child is infeasible (>= 0) subtract 1 if both Optimal
    # can be negative temporately in step_obj but not in tree.obj_gain_
    inf_counter     :: Vector{Int64} 
end

mutable struct BnBTreeObj
    m               :: Juniper.JuniperProblem
    incumbent       :: Incumbent
    limit           :: Symbol
    obj_gain        :: GainObj
    disc2var_idx    :: Vector{Int64}
    var2disc_idx    :: Vector{Int64}
    options         :: Juniper.SolverOptions
    obj_fac         :: Int64 # factor for objective 1 if max -1 if min
    start_time      :: Float64 
    nsolutions      :: Int64
    branch_nodes    :: Vector{BnBNode}
    best_bound      :: Float64
    global_solver   :: Bool

    BnBTreeObj() = new()
end

mutable struct StrongBranchStep
    var_idx             :: Int64
    l_relaxation_state  :: MOI.TerminationStatusCode
    r_relaxation_state  :: MOI.TerminationStatusCode
    init_restart        :: Bool

    StrongBranchStep() = new()
end

# the object holds information for the current step
mutable struct StepObj
    node                :: BnBNode # current branch node
    var_idx             :: Int64   # variable to branch on
    state               :: Symbol  # if infeasible => break (might be set by strong branching)
    nrestarts           :: Int64 
    gain_gap            :: Float64
    obj_gain            :: GainObj
    idx_time            :: Float64
    node_idx_time       :: Float64
    upd_gains_time      :: Float64
    node_branch_time    :: Float64
    branch_time         :: Float64
    integral            :: Vector{BnBNode}
    branch              :: Vector{BnBNode}
    l_nd                :: BnBNode
    r_nd                :: BnBNode
    counter             :: Int64
    upd_gains           :: Symbol
    strong_disc_vars    :: Vector{Int64}
    branch_strategy     :: Symbol # the strategy that is used in this step
    strong_branching    :: Union{Nothing,Vector{StrongBranchStep}} # used only if debug=true and branch_strategy of this step is :Strong

    StepObj() = new()
end

function Base.copy(s::StepObj) 
    step_obj = StepObj()
    for f in fieldnames(StepObj)
        if isdefined(s,f)
            if isa(getfield(s,f), BnBNode) || isa(getfield(s,f), Vector)
                setfield!(step_obj, f, copy(getfield(s,f)))
            else
                setfield!(step_obj, f, getfield(s,f))
            end
        end
    end

    return step_obj
end

mutable struct StrongObj
    rank                    :: Int64
    step_obj                :: StepObj
    restart                 :: Bool
    max_gain                :: Float64
    max_gain_var            :: Int64
    max_gain_disc_var       :: Int64
    need_to_resolve         :: Bool
    left_node               :: Union{Nothing, BnBNode}
    right_node              :: Union{Nothing, BnBNode}
    gains                   :: GainObj
    set_to_last_var         :: Bool
    infeasible_disc_vars    :: Vector{Int64}

    StrongObj() = new()
end

mutable struct TimeObj
    solve_leaves_get_idx :: Float64
    solve_leaves_branch :: Float64
    branch :: Float64
    get_idx :: Float64
    upd_gains :: Float64
end
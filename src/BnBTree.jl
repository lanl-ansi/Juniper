include("table_log.jl")
include("bb_inits_and_defaults.jl")
include("bb_strategies.jl")
include("bb_user_limits.jl")
include("bb_type_correct.jl")
include("bb_integral_or_branch.jl")
include("bb_gains.jl")

"""
    get_branching_disc_idx!(m, step_obj, opts, disc2var_idx, gains, counter::Int64=1)

Get the index of a variable to branch on.
"""
function get_branching_disc_idx!(m, step_obj, opts, disc2var_idx, gains, counter::Int64=1)
    start = time()
    node = step_obj.node
    idx = 0
    strong_restarts = 0
    branch_strat = opts.branch_strategy
    status = :Normal
    if branch_strat == :MostInfeasible
        idx = branch_mostinfeasible(m, node, disc2var_idx)
        step_obj.branch_strategy = :MostInfeasible
    elseif branch_strat == :PseudoCost || branch_strat == :StrongPseudoCost
        if counter == 1 && branch_strat == :PseudoCost
            idx = branch_mostinfeasible(m, node, disc2var_idx)
            step_obj.branch_strategy = :MostInfeasible
        elseif counter <= opts.strong_branching_nsteps && branch_strat == :StrongPseudoCost
            status, idx, strong_restarts = branch_strong!(m, opts, disc2var_idx, step_obj, counter)
            step_obj.branch_strategy = :Strong
        else
            idx = branch_pseudo(m, node, disc2var_idx, gains, opts.gain_mu, opts.atol)
            step_obj.branch_strategy = :Pseudo
        end
    elseif branch_strat == :Reliability
        status, idx, strong_restarts = branch_reliable!(m,opts,step_obj,disc2var_idx,gains,counter)
        step_obj.branch_strategy = :Reliability
    end
    step_obj.state = status
    step_obj.var_idx = idx
    step_obj.nrestarts = strong_restarts
    step_obj.idx_time = time()-start
    return
end

"""
    optimize_node(m1, node, parent; restart=false)

Solves one node on the current processor. Bounds are set and the starting value is the solution of the parent.
If m1 is nothing use the global variable m (if called by a processor)
Return a NodeSolution object containing the result
"""
function optimize_node(m1, node, parent; restart=false)
    if m1 == nothing
        global m
    else
        m = m1
    end

    # set bounds
    for i=1:m.num_var
        JuMP.set_lower_bound(m.x[i], node.l_var[i])
        JuMP.set_upper_bound(m.x[i], node.u_var[i])
    end
    if restart > 0
        restart_values = generate_random_restart(m)
        JuMP.set_start_value.(m.x[1:m.num_var], restart_values)
    else
        JuMP.set_start_value.(m.x[1:m.num_var], parent.solution)
    end

    old_mu_init = set_subsolver_option!(m, m.model, "nl", "Ipopt", :mu_init, 0.1 => 1e-5)                                   

    status, backend = optimize_get_status_backend(m.model)

    # reset mu_init
    reset_subsolver_option!(m, "nl", "Ipopt", :mu_init, old_mu_init)
    
    objval = NaN
    solution = fill(NaN, m.num_var)
    if state_is_optimal(status; allow_almost=m.options.allow_almost_solved)
        objval = JuMP.objective_value(m.model)
        solution = JuMP.value.(m.x)
    end

    result = NodeSolution(status, objval, solution)
    return result
end


"""
    classify_node!(m, step_obj, node::BnBNode, solution::NodeSolution; temp=false)

Change the state of the node to :Integral, :Branch, :Error, ... and at to step_obj if !temp
"""
function classify_node!(m, step_obj, node::BnBNode, solution::NodeSolution; temp=false)
    status = solution.status
    node.relaxation_state = status
    objval = solution.objval
    node.solution = solution.solution

    if !state_is_optimal(status; allow_almost=m.options.allow_almost_solved) && !state_is_infeasible(status)
        node.state = :Error
    elseif state_is_optimal(status; allow_almost=m.options.allow_almost_solved)
        node.best_bound = objval
        set_cnode_state!(node, m, step_obj, m.disc2var_idx)
        if only_almost_solved(status) && (!(node.state == :Integral) || m.options.allow_almost_solved_integral) 
            @warn "Only almost solved"
        end

        # if almost_solved_integral is not allowed but it is integral change state
        if node.state == :Integral && only_almost_solved(status) && !m.options.allow_almost_solved_integral
            node.state = :Almost_Solved
            # Todo option for restart
        end
        if !temp
            push_integral_or_branch!(step_obj, node)
        end
    else
        node.state = :Infeasible
    end
end


"""
    branch!(m, opts, step_obj, counter, disc2var_idx; temp=false)

Branch a node by using x[idx] <= floor(x[idx]) and x[idx] >= ceil(x[idx])
Solve both nodes and set current node state to done.
"""
function branch!(m, opts, step_obj, counter, disc2var_idx; temp=false)
    if myid() != 1
        global procs_available
    else
        procs_available = [1]
    end

    leaf_idx = 0
    function get_next_child_node_idx(child_nodes::Vector{BnBNode})
        leaf_idx += 1
        if leaf_idx > length(child_nodes)
            return 0
        end
        return leaf_idx
    end

    ps = opts.log_levels
    node = step_obj.node
    vidx = step_obj.var_idx

    start = time()
    # it might be already branched on
    if !temp && node.state != :Branch
        for cnode in [step_obj.l_nd,step_obj.r_nd]
            if cnode.state == :Branch || cnode.state == :Integral
                push_integral_or_branch!(step_obj, cnode)
            end
        end
        return step_obj.l_nd,step_obj.r_nd
    end

    l_nd_u_var = copy(node.u_var)
    r_nd_l_var = copy(node.l_var)
    l_nd_u_var[vidx] = floor(node.solution[vidx])
    r_nd_l_var[vidx] = ceil(node.solution[vidx])

    if opts.debug
        path_l = copy(node.path)
        path_r = copy(node.path)
        push!(path_l,node.hash)
        push!(path_r,node.hash)
    else
        path_l = []
        path_r = []
    end
    l_nd = new_left_node(node, l_nd_u_var; path=path_l)
    r_nd = new_right_node(node, r_nd_l_var; path=path_r)

    # save that this node branches on this particular variable
    node.var_idx = vidx

    check_print(ps,[:All,:FuncCall]) && println("branch")

    if !temp
        step_obj.l_nd = l_nd
        step_obj.r_nd = r_nd
        node.state = :Done
    end

    start_process = time()
    child_nodes = BnBNode[l_nd,r_nd]
    @sync begin
        for p in procs_available
            @async begin
                while true
                    child_node_idx = get_next_child_node_idx(child_nodes)
                    if child_node_idx > 0
                        child_node = child_nodes[child_node_idx]
                        m_param = myid() == 1 ? m : nothing
                        child_solution = remotecall_fetch(optimize_node, p, m_param, child_node, step_obj.node)
                        classify_node!(m, step_obj, child_node, child_solution; temp=temp)
                    else
                        break
                    end
                end
            end
        end
    end
    node_time = time() - start_process

    if !temp
        step_obj.node_branch_time += node_time
    end

    if !temp
        step_obj.branch_time += time()-start
    end

    @assert child_nodes[1].relaxation_state != MOI.OPTIMIZE_NOT_CALLED
    @assert child_nodes[2].relaxation_state != MOI.OPTIMIZE_NOT_CALLED
    return child_nodes[1], child_nodes[2]
end

"""
    update_incumbent!(tree::BnBTreeObj, node::BnBNode)

Gets called if new integral solution was found.
Check whether it's a new incumbent and update if necessary
"""
function update_incumbent!(tree::BnBTreeObj, node::BnBNode)
    ps = tree.options.log_levels
    check_print(ps,[:All,:FuncCall]) && println("update_incumbent")

    factor = tree.obj_fac
    if !isdefined(tree,:incumbent) || factor*node.best_bound > factor*tree.incumbent.objval
        objval = node.best_bound
        solution = copy(node.solution)
        status = node.relaxation_state
        @assert(state_is_optimal(status; allow_almost=false) || (only_almost_solved(status) && tree.options.allow_almost_solved_integral))
        tree.incumbent = Incumbent(objval, solution, only_almost_solved(status))
        if !tree.options.all_solutions
            bound!(tree)
        end
        return true
    end

    return false
end


"""
    bound!(tree::BnBTreeObj)
"""
function bound!(tree::BnBTreeObj)
    function isbetter(n)
        return f*n.best_bound > f*incumbent_val
    end
    incumbent_val = tree.incumbent.objval
    f = tree.obj_fac
    filter!(isbetter, tree.branch_nodes)
end

"""
    add_incumbent_constr(m, incumbent)

Add a constraint >=/<= incumbent
"""
function add_incumbent_constr(m, incumbent)
    add_obj_constraint(m, incumbent.objval)
end

"""
    add_obj_epsilon_constr(tree)

Add a constraint obj ≦ (1+ϵ)*LB or obj ≧ (1-ϵ)*UB
"""
function add_obj_epsilon_constr(tree)
    # add constr for objval
    if tree.options.obj_epsilon > 0
        ϵ = tree.options.obj_epsilon
        if tree.m.obj_sense == :Min
            rhs = (1+ϵ)*tree.m.relaxation_objval
        else
            rhs = (1-ϵ)*tree.m.relaxation_objval
        end
        add_obj_constraint(tree.m, rhs)
        tree.m.ncuts += 1
    end
end

"""
    get_next_branch_node!(tree)

Get the next branch node (sorted by best bound)
Return true,step_obj if there is a branch node and
false, nothing otherwise
"""
function get_next_branch_node!(tree)
    if length(tree.branch_nodes) == 0
        return false, false, nothing
    end
    bvalue, nidx = findmax([tree.obj_fac*n.best_bound for n in tree.branch_nodes])

    trav_strat = tree.options.traverse_strategy
    if trav_strat == :BFS
        # use the one with highest depth
        _value, nidx = findmax([if isapprox(tree.obj_fac*n.best_bound, bvalue) n.level else -Inf end for n in tree.branch_nodes])
    end

    if trav_strat == :DFS || (trav_strat == :DBFS && !isdefined(tree,:incumbent))
        _value, nidx = findmax([n.level for n in tree.branch_nodes])
    end

    tree.best_bound = tree.obj_fac*bvalue
    bbreak = isbreak_mip_gap(tree)
    if bbreak
        return false, true, nothing
    end

    node = tree.branch_nodes[nidx]
    deleteat!(tree.branch_nodes, nidx)

    return true, false, new_default_step_obj(tree.m,node)
end


"""
    one_branch_step!(m1, incumbent, opts, step_obj, disc2var_idx, gains, counter)

Get a branch variable using the specified strategy and branch on the node in step_obj
using that variable. Return the new updated step_obj
"""
function one_branch_step!(m1, incumbent, opts, step_obj, disc2var_idx, gains, counter)
    if m1 === nothing
        global m
        global is_newincumbent
        if opts.incumbent_constr && incumbent !== nothing && is_newincumbent
            is_newincumbent = false
            add_incumbent_constr(m, incumbent)
        end
    else
        m = m1
    end

    node = step_obj.node
    step_obj.counter = counter

# get branch variable
    node_idx_start = time()
    get_branching_disc_idx!(m, step_obj, opts, disc2var_idx, gains, counter)
    step_obj.node_idx_time = time()-node_idx_start
    # if no variable got selected might be true that all variables are already type correct
    # node.solution can be updated if one child is infeasible and the other optimal (can result in discrete)
    if step_obj.var_idx == 0 && are_type_correct(step_obj.node.solution, m.var_type, disc2var_idx, opts.atol)
        push!(step_obj.integral, node)
    else
        if step_obj.state != :Infeasible && step_obj.state != :PartlyInfeasible
            @assert step_obj.var_idx != 0
            branch!(m, opts, step_obj, counter, disc2var_idx)
        end
    end
    return step_obj
end

"""
    upd_time_obj!(time_obj, step_obj)

Add step_obj times to time_obj
"""
function upd_time_obj!(time_obj, step_obj)
    time_obj.solve_leaves_get_idx += step_obj.node_idx_time
    time_obj.solve_leaves_branch += step_obj.node_branch_time
    time_obj.branch += step_obj.branch_time
    time_obj.get_idx += step_obj.idx_time
    time_obj.upd_gains += step_obj.upd_gains_time
end

"""
    upd_tree_obj!(tree, step_obj, time_obj)

Update the tree obj like new incumbent or new branch nodes using the step_obj
Return true if it's the end of the algorithm (checking different break rules)
"""
function upd_tree_obj!(tree, step_obj, time_obj)
    node = step_obj.node
    still_running = true

    if step_obj.node.level+1 > tree.m.nlevels
        tree.m.nlevels = step_obj.node.level+1
    end

    if step_obj.state == :Infeasible
        # if there is no incumbent yet
        if !isdefined(tree,:incumbent)
            # this can be locally infeasible or globally infeasible depending on the solver
            tree.limit = :Infeasible
        end # it will terminate and use the current solution as optimal (might want to rerun as an option)
        still_running = false
    end

    if still_running
        upd_gains_step!(tree, step_obj)
    end

    bbreak = upd_integral_branch!(tree, step_obj)
    if bbreak
        still_running = false
    end

    upd_time_obj!(time_obj, step_obj)
    return !still_running
end

"""
    solve_sequential(tree,
        last_table_arr,
        time_bnb_solve_start,
        fields,
        field_chars,
        time_obj)

Run branch and bound on a single processor
"""
function solve_sequential(tree,
    last_table_arr,
    time_bnb_solve_start,
    fields,
    field_chars,
    time_obj)

    m = tree.m
    opts = tree.options
    disc2var_idx = tree.disc2var_idx
    counter = 0
    ps = tree.options.log_levels

    dictTree = Dict{Any,Any}()
    while true
        # the _ is only needed for parallel
        exists, _, step_obj = get_next_branch_node!(tree)
        !exists && break
        isbreak_after_step!(tree) && break
        counter += 1
        if isdefined(tree,:incumbent)
            step_obj = one_branch_step!(m, tree.incumbent, opts, step_obj, disc2var_idx, tree.obj_gain, counter)
        else
            step_obj = one_branch_step!(m, nothing, opts, step_obj, disc2var_idx, tree.obj_gain, counter)
        end
        m.nnodes += 2 # two nodes explored per branch
        node = step_obj.node

        bbreak = upd_tree_obj!(tree,step_obj,time_obj)
        tree.options.debug && (dictTree = push_step2treeDict!(dictTree,step_obj))

        if check_print(ps,[:Table])
            last_table_arr = print_table([1],tree,node,step_obj,time_bnb_solve_start,fields,field_chars;last_arr=last_table_arr)
        end

        if bbreak
            break
        end
    end

    tree.options.debug && (tree.m.debugDict[:tree] = dictTree)
    return counter
end

function sendto(p::Int; args...)
    for (nm, val) in args
        remotecall_fetch(Core.eval, p, Juniper, Expr(:(=), nm, val))
    end
end

function dummysolve()
    global m
    optimize!(m.model)
end

"""
    pmap(f, tree, last_table_arr, time_bnb_solve_start,
        fields, field_chars, time_obj)

Run the solving steps on several processors
"""
function pmap(f, tree, last_table_arr, time_bnb_solve_start,
    fields, field_chars, time_obj)
    np = nprocs()  # determine the number of processes available
    if tree.options.processors+1 < np
        np = tree.options.processors+1
    end

    # function to produce the next work item from the queue.
    # in this case it's just an index.
    ps = tree.options.log_levels
    still_running = true
    run_counter = 0
    counter = 0

    for p=2:np
        seed = Random.seed!
        remotecall_fetch(seed, p, 1)
    end
    @sync begin
        for p=2:np
            @async begin
                sendto(p, m=tree.m)
                sendto(p, is_newincumbent=false)
                if tree.options.two_processors_per_node && p % 2 == 0 && p+1 <= np
                    sendto(p, procs_available=[p+1,p])
                elseif !tree.options.two_processors_per_node || p+1 > np
                    sendto(p, procs_available=[p])
                end
            end
        end
    end

    for p=3:np
        remotecall(dummysolve, p)
    end

    proc_counter = zeros(np)

    dictTree = Dict{Any,Any}()

    branch_strat = tree.options.branch_strategy
    opts = tree.options
    # starting with two as one is master anyway
    list_procs = 2:np
    if opts.two_processors_per_node 
        list_procs = 2:2:np
    end
    @sync begin
        for p=list_procs
            @async begin
                while true
                    exists, bbreak, step_obj = get_next_branch_node!(tree)
                    if bbreak
                        still_running = false
                        break
                    end

                    while !exists && still_running
                        sleep(0.1)
                        exists, bbreak, step_obj = get_next_branch_node!(tree)
                        exists && break
                        if bbreak
                            still_running = false
                            break
                        end
                    end
                    !still_running && break

                    if isbreak_after_step!(tree)
                        still_running = false
                        break
                    end

                    run_counter += 1
                    counter += 1
                    if isdefined(tree,:incumbent)
                        step_obj = remotecall_fetch(f, p, nothing, tree.incumbent, tree.options, step_obj,
                                                    tree.disc2var_idx, tree.obj_gain, counter)
                    else
                        step_obj = remotecall_fetch(f, p, nothing, nothing, tree.options, step_obj,
                        tree.disc2var_idx, tree.obj_gain, counter)
                    end
                    tree.m.nnodes += 2 # two nodes explored per branch
                    run_counter -= 1

                    !still_running && break

                    bbreak = upd_tree_obj!(tree,step_obj,time_obj)
                    tree.options.debug && (dictTree = push_step2treeDict!(dictTree,step_obj))

                    if run_counter == 0 && length(tree.branch_nodes) == 0
                        still_running = false
                    end

                    bbreak && (still_running = false)

                    if check_print(ps,[:Table])
                        if tree.options.two_processors_per_node && p % 2 == 0 && p+1 <= np
                            procs_available=[p,p+1]
                        else
                            procs_available=[p]
                        end

                        if length(tree.branch_nodes) > 0
                            bvalue, nidx = findmax([tree.obj_fac*n.best_bound for n in tree.branch_nodes])
                            tree.best_bound = tree.obj_fac*bvalue
                        end
                        last_table_arr = print_table(procs_available,tree,step_obj.node,step_obj,time_bnb_solve_start,fields,field_chars;last_arr=last_table_arr)
                    end

                    !still_running && break
                end
            end
        end
    end
    tree.options.debug && (tree.m.debugDict[:tree] = dictTree)
    return counter
end

"""
    solvemip(tree::BnBTreeObj)

Solve the MIP part of a problem given by BnBTreeObj using branch and bound.
 - Identify the node to branch on
 - Get variable to branch on
 - Solve subproblems
"""
function solvemip(tree::BnBTreeObj)
    time_obj = init_time_obj()
    time_bnb_solve_start = time()

    ps = tree.options.log_levels

    # check if already integral
    if are_type_correct(tree.m.relaxation_solution,tree.m.var_type,tree.disc2var_idx, tree.options.atol)
        tree.nsolutions = 1
        objval = JuMP.objective_value(tree.m.model)
        sol = JuMP.value.(tree.m.x)
        bbound = try 
            JuMP.objective_bound(tree.m.model)
        catch
            JuMP.objective_value(tree.m.model)
        end
        tree.incumbent = Incumbent(objval,sol,only_almost_solved(tree.m.status))
        return tree.incumbent
    end

    # check if incumbent and if mip gap already fulfilled (i.e found by fp)
    if isbreak_mip_gap(tree)
        return tree.incumbent
    end

    last_table_arr = []
    fields = []
    field_chars = []
    # Print table init
    if check_print(ps,[:Table])
        fields, field_chars = get_table_config(tree.options)
        print_table_header(fields,field_chars)
    end


    check_print(ps,[:All,:FuncCall]) && println("Solve Tree")

    branch_strat = tree.options.branch_strategy

    add_obj_epsilon_constr(tree)

    # use pmap if more then one processor
    if tree.options.processors > 1 || tree.options.force_parallel
        counter = pmap(Juniper.one_branch_step!,tree,
            last_table_arr,
            time_bnb_solve_start,
            fields,
            field_chars,
            time_obj)
    else
        counter = solve_sequential(tree,
            last_table_arr,
            time_bnb_solve_start,
            fields,
            field_chars,
            time_obj)
    end

    if !isdefined(tree,:incumbent) && tree.limit == :None
        # infeasible
        tree.limit = :Infeasible
    end

    if tree.options.obj_epsilon != 0 && tree.limit == :Infeasible
        @warn "Maybe only infeasible because of obj_epsilon."
    end

    if !isnan(tree.options.best_obj_stop) && isdefined(tree, :incumbent)
        inc_val = tree.incumbent.objval
        bos = tree.options.best_obj_stop
        sense = tree.m.obj_sense
        if (sense == :Min && inc_val > bos) || (sense == :Max && inc_val < bos)
            @warn "best_obj_gap couldn't be reached."
        end
    end

    if length(tree.branch_nodes) > 0
        bvalue, nidx = findmax([tree.obj_fac*n.best_bound for n in tree.branch_nodes])
        tree.best_bound = tree.obj_fac*bvalue
    else
        if isdefined(tree, :incumbent) && !isnan(tree.incumbent.objval)
            tree.best_bound = tree.incumbent.objval
        end
    end

    tree.m.nbranches = counter

    time_bnb_solve = time()-time_bnb_solve_start
    check_print(ps,[:Table]) && println("")
    check_print(ps,[:All,:Info]) && println("#branches: ", counter)

    tree.options.debug && debug_set_tree_obj_gain!(tree)
    check_print(ps,[:All,:Timing]) && print_final_timing(time_bnb_solve, time_obj)
end
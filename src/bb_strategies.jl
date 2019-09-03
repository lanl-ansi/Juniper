"""
branch_mostinfeasible(m, node, disc2var_idx)

Get the index of an integer variable which is currently continuous which is most unintegral.
(nearest to *.5)
"""
function branch_mostinfeasible(m, node, disc2var_idx)
    x = node.solution
    idx = 0
    max_diff = 0
    for i in disc2var_idx
        diff = abs(x[i]-round(x[i]))
        if diff > max_diff
            idx = i
            max_diff = diff
        end
    end
    return idx
end

"""
init_strong_restart!(node, opts, var_idx, disc_var_idx, l_nd, r_nd, reasonable_disc_vars, infeasible_disc_vars,
 left_node, right_node, strong_restart, max_gain_var, atol, var_type, disc2var_idx)

Tighten the bounds for the node and check if there are variables that need to be checked for a restart.
"""
function init_strong_restart!(node, opts, var_idx, disc_var_idx, l_nd, r_nd, 
                                reasonable_disc_vars, infeasible_disc_vars, 
                                left_node, right_node, strong_restart, max_gain_var,
                                atol, var_type, disc2var_idx)
    restart = false
    set_to_last_var = false
    need_to_resolve = false
    new_infeasible_disc_vars = copy(infeasible_disc_vars)

    # set the bounds directly for the node
    # also update the best bound and the solution
    if !state_is_optimal(l_nd.relaxation_state; allow_almost=opts.allow_almost_solved)
        # the solution shouldn't be updated when the current max_gain_var is then type correct
        if max_gain_var == 0 || !is_type_correct(r_nd.solution[max_gain_var],var_type[max_gain_var],atol)
            # if all variables reasonable_disc_vars are type correct => just continue with the tree
            # no restarts and no change in bounds
            if !all_reasonable_type_correct(r_nd.solution, disc2var_idx, reasonable_disc_vars, atol)
                node.l_var[var_idx] = ceil(node.solution[var_idx])
                node.best_bound = r_nd.best_bound
                node.solution = r_nd.solution
                need_to_resolve = true
                push!(new_infeasible_disc_vars, disc_var_idx)
                strong_restart && (restart = true)
            end
        end
    else
        # the solution shouldn't be updated when the current max_gain_var is then type correct
        if max_gain_var == 0 || !is_type_correct(l_nd.solution[max_gain_var],var_type[max_gain_var],atol)
            if !all_reasonable_type_correct(l_nd.solution, disc2var_idx, reasonable_disc_vars, atol)
                node.u_var[var_idx] = floor(node.solution[var_idx])
                node.best_bound = l_nd.best_bound
                node.solution = l_nd.solution
                need_to_resolve = true
                push!(new_infeasible_disc_vars, disc_var_idx)
                strong_restart && (restart = true)
            end
        end
    end


    if length(reasonable_disc_vars) == length(new_infeasible_disc_vars)
        # basically branching on the last infeasible variable 
        set_to_last_var = true
        restart = false
    end
    return need_to_resolve, restart, new_infeasible_disc_vars, set_to_last_var
end

"""
    branch_strong_on!(m,opts,step_obj,reasonable_disc_vars, disc2var_idx, strong_restart, counter)

Try to branch on a few different variables and choose the one with highest obj_gain.
Update obj_gain for the variables tried and average the other ones.
"""
function branch_strong_on!(m,opts,step_obj,
    reasonable_disc_vars, disc2var_idx, strong_restart, counter)

    function set_temp_gains!(gains, gain_l, gain_r, disc_var_idx)
        if !isinf(gain_l)
            gains.minus[disc_var_idx] = gain_l
            gains.minus_counter[disc_var_idx] = 1
        else 
            gains.inf_counter[disc_var_idx] = 1
        end
        if !isinf(gain_r)
            gains.plus[disc_var_idx] = gain_r
            gains.plus_counter[disc_var_idx] = 1
        else 
            gains.inf_counter[disc_var_idx] = 1
        end

        if !isinf(gain_l) && !isinf(gain_r)
            gains.inf_counter[disc_var_idx] -= 1
        end
    end

    function get_current_gains(opts, node, l_nd, r_nd)
        gain_l = sigma_minus(opts, node, l_nd, node.solution[node.var_idx])
        gain_r = sigma_plus( opts, node,  r_nd, node.solution[node.var_idx])

        if isinf(gain_l)
            gain = gain_r
        elseif isinf(gain_r)
            gain = gain_l
        else
            gain = (gain_l+gain_r)/2
        end
        if isnan(gain)
            gain = Inf
        end
        return gain_l, gain_r, gain
    end

    function add_strong_step(strong_step,l_nd,r_nd; restart=false) 
        strong_step.l_relaxation_state = l_nd.relaxation_state
        strong_step.r_relaxation_state = r_nd.relaxation_state
        strong_step.init_restart = restart
        push!(step_obj.strong_branching, strong_step)
    end

    current_disc_var_idx = 0
    function get_next_disc_var(step_obj, max_gain, max_gain_var, max_gain_disc_var)
        current_disc_var_idx += 1
        if current_disc_var_idx > length(reasonable_disc_vars)
            return nothing, 0, 0
        end
        strong_obj = StrongObj()
        strong_obj.rank = current_disc_var_idx
        strong_obj.step_obj = copy(step_obj)
        strong_obj.restart = false 
        strong_obj.max_gain = max_gain
        strong_obj.max_gain_var = max_gain_var
        strong_obj.max_gain_disc_var = max_gain_disc_var
        strong_obj.need_to_resolve = false
        strong_obj.left_node = nothing
        strong_obj.right_node = nothing
        strong_obj.gains = init_gains(m.num_disc_var)
        strong_obj.set_to_last_var = false
        strong_obj.infeasible_disc_vars = zeros(Int64,0)
        return strong_obj, reasonable_disc_vars[current_disc_var_idx], current_disc_var_idx
    end

    np = nprocs()  # determine the number of processes available
    if opts.processors+1 < np
        np = opts.processors+1
    end

    strong_time = time()

    strong_restarts = -1
    restart = true
    status = :Normal
    infeasible_disc_vars = zeros(Int64,0)
    gains = init_gains(m.num_disc_var)
  
    left_node, right_node = nothing, nothing
    atol = opts.atol

    need_to_resolve = false
    opts.debug && (step_obj.strong_branching = [])

    max_gain = -Inf # then one is definitely better
    max_gain_var = 0
    max_gain_disc_var = 0

    allow_almost = opts.allow_almost_solved
    while restart 
        strong_objs = Vector{StrongObj}()
        restart = false
        current_disc_var_idx = 0
        strong_restarts += 1 # is init with -1
        restart_rank = -1
        parallel_break = false
        if strong_restarts > 0
            # println("START/RESTART")
        end
        @sync begin
            for p=1:np
                if p != myid() || np == 1
                    @async begin
                        while true
                            parallel_break && break
                            # rank is used so that the first rank is used for restart purposes
                            # this makes it like sequential  
                            strong_obj, disc_var_idx, rank = get_next_disc_var(step_obj, max_gain, max_gain_var, max_gain_disc_var)
                            # returns 0 if nothing left
                            if disc_var_idx == 0
                                break
                            end
                            strong_step_obj = strong_obj.step_obj
                            node = strong_step_obj.node

                            # don't rerun if the variable has already one infeasible node
                            if disc_var_idx in infeasible_disc_vars
                                var_idx = disc2var_idx[disc_var_idx]
                                continue
                            end
                            var_idx = disc2var_idx[disc_var_idx]
                            strong_step_obj.var_idx = var_idx
                            u_b, l_b = node.u_var[var_idx], node.l_var[var_idx]
                            # don't rerun if bounds are exact or is type correct
                            if isapprox(u_b,l_b; atol=atol) || is_type_correct(node.solution[var_idx],m.var_type[var_idx], atol)
                                continue
                            end
                            push!(strong_objs, strong_obj)
                            
                            # branch on the current variable and get the corresponding children
                            l_nd,r_nd = remotecall_fetch(branch!, p, m, opts, strong_step_obj, counter, disc2var_idx; temp=true)
                            # in parallel the above call is by value not reference :/ => update node.var_idx here
                            node.var_idx = var_idx

                            opts.debug && (strong_step = new_default_strong_branch_step_obj(var_idx))
                            
                            # no current restart => we can set max_gain and variable
                            gain_l, gain_r, gain = get_current_gains(opts, node, l_nd, r_nd)

                            if !state_is_optimal(l_nd.relaxation_state; allow_almost=allow_almost) && !state_is_optimal(r_nd.relaxation_state; allow_almost=allow_almost) && counter == 1
                                # TODO: Might be Error/UserLimit instead of infeasible
                                status = :Infeasible
                                strong_obj.left_node = l_nd
                                strong_obj.right_node = r_nd
                                opts.debug && add_strong_step(strong_step,l_nd,r_nd)
                                parallel_break = true
                                break
                            end

                            # check if one part is infeasible => update bounds & restart if strong restart is true
                            if !state_is_optimal(l_nd.relaxation_state; allow_almost=allow_almost) || !state_is_optimal(r_nd.relaxation_state; allow_almost=allow_almost)
                                if !state_is_optimal(l_nd.relaxation_state; allow_almost=allow_almost) && !state_is_optimal(r_nd.relaxation_state; allow_almost=allow_almost)
                                    # TODO: Might be Error/UserLimit instead of infeasible
                                    status = :PartlyInfeasible
                                    strong_obj.left_node = l_nd
                                    strong_obj.right_node = r_nd
                                    opts.debug && add_strong_step(strong_step,l_nd,r_nd)
                                    parallel_break = true
                                    break
                                end
                                local_need_to_resolve, local_restart, new_infeasible_disc_vars, local_set_to_last_var = init_strong_restart!(node, opts, var_idx, disc_var_idx, l_nd, r_nd, reasonable_disc_vars, infeasible_disc_vars, strong_obj.left_node, strong_obj.right_node, strong_restart, strong_obj.max_gain_var, opts.atol, m.var_type, disc2var_idx)
                                strong_obj.need_to_resolve = local_need_to_resolve
                                strong_obj.restart = local_restart
                                strong_obj.set_to_last_var = local_set_to_last_var
                                # that is already a copy
                                strong_obj.infeasible_disc_vars = new_infeasible_disc_vars

                                # only variables where one branch is infeasible => no restart and break
                                if strong_obj.set_to_last_var
                                    strong_obj.max_gain_var = var_idx
                                    strong_obj.max_gain_disc_var = disc_var_idx
                                    strong_obj.left_node = l_nd
                                    strong_obj.right_node = r_nd
                                    strong_obj.restart = false
                                    strong_obj.need_to_resolve = false
                                    gain_l, gain_r, gain = get_current_gains(opts, node, l_nd, r_nd)
                                    strong_obj.max_gain = gain
                                    set_temp_gains!(strong_obj.gains, gain_l, gain_r, disc_var_idx)
                                    opts.debug && add_strong_step(strong_step,l_nd,r_nd)
                                    parallel_break = true
                                    break
                                end

                                if strong_obj.restart && time()-strong_time > opts.strong_branching_approx_time_limit
                                    strong_obj.restart = false
                                end
                            end

                            # don't update maximum if restart
                            if gain > strong_obj.max_gain && !strong_obj.restart
                                strong_obj.max_gain = gain
                                strong_obj.max_gain_var = var_idx
                                strong_obj.max_gain_disc_var = disc_var_idx
                                strong_obj.left_node = l_nd
                                strong_obj.right_node = r_nd
                                # we are changing the bounds if one branch is infeasible then the current
                                # selection might not work anymore => we have to resolve it later
                                strong_obj.need_to_resolve = false
                            end
                            set_temp_gains!(strong_obj.gains, gain_l, gain_r, disc_var_idx)
                            opts.debug && add_strong_step(strong_step,l_nd,r_nd;restart=strong_obj.restart)
                            if strong_obj.restart
                                parallel_break = true
                                break
                            end

                            if time()-m.start_time >= opts.time_limit
                                parallel_break = true
                                break
                            end
                        end
                    end
                end
            end
        end

        # end of run (maybe terminated with restart)
        # here we combine the parallel runs to something we can actually use
        # if restart we take the one which called restart the first
        sort!(strong_objs, by= so->so.rank)
        gains = init_gains(m.num_disc_var)
        # println("RECOMBINATION")
        for so in strong_objs
            #= println("rank: ", so.rank)
            println("restart: ", so.restart)
            println("need_to_resolve: ", so.need_to_resolve) =#
            gains += so.gains

            if so.restart
                restart = true
                need_to_resolve = so.need_to_resolve
                max_gain = so.max_gain
                max_gain_var = so.max_gain_var
                max_gain_disc_var = so.max_gain_disc_var
                step_obj = copy(so.step_obj)
                infeasible_disc_vars = copy(so.infeasible_disc_vars)
                break
            end
        end

        if !restart
            max_gain = -Inf
            for so in strong_objs
                if so.max_gain > max_gain
                    #= println("Setting to: ")
                    println("so.max_gain: ", so.max_gain)
                    println("so.max_gain_var: ", so.max_gain_var)
                    println("so.need_to_resolve: ", so.need_to_resolve) =#
                    need_to_resolve = so.need_to_resolve
                    max_gain = so.max_gain
                    max_gain_var = so.max_gain_var
                    max_gain_disc_var = so.max_gain_disc_var
                    left_node = so.left_node
                    right_node = so.right_node
                    step_obj = copy(so.step_obj)
                end
            end
            if max_gain == -Inf
                @assert status == :PartlyInfeasible || status == :Infeasible
                for so in strong_objs
                    if !state_is_optimal(so.left_node.relaxation_state) && !state_is_optimal(so.right_node.relaxation_state)
                        left_node = so.left_node
                        right_node = so.right_node
                        step_obj = copy(so.step_obj)
                        break
                    end
                end
            end
        end
    end

    # check if need to resolve (not if Local/Global Infeasible)
    if need_to_resolve && status == :Normal 
        status = :Resolve
    end

    return status, max_gain_var, left_node, right_node, gains, strong_restarts
end

"""
branch_strong!(m,opts,disc2var_idx,step_obj,counter)

Try to branch on a few different variables and choose the one with highest obj_gain.
Update obj_gain for the variables tried and average the other ones.
"""
function branch_strong!(m,opts,disc2var_idx,step_obj,counter)
    node = step_obj.node

    # generate an of variables to branch on
    num_strong_var = Int(round((opts.strong_branching_perc/100)*m.num_disc_var))
    # if smaller than 2 it doesn't make sense
    num_strong_var = num_strong_var < 2 ? 2 : num_strong_var
    # use strong_branching_approx_time_limit to change num_strong_var
    if !isinf(opts.strong_branching_approx_time_limit)
        approx_time_per_node = 2*m.relaxation_time
        new_num_strong_var = Int(floor(opts.strong_branching_approx_time_limit/approx_time_per_node))
        new_num_strong_var = new_num_strong_var == 0 ? 1 : new_num_strong_var
        if new_num_strong_var < num_strong_var
            @warn "Changed num_strong_var to $new_num_strong_var because of strong_branching_approx_time_limit"
            num_strong_var = new_num_strong_var
        end
    end

    # get reasonable candidates (not type correct and not already perfectly bounded)
    disc_vars = m.num_disc_var
    atol = opts.atol
    reasonable_disc_vars = get_reasonable_disc_vars(node, m.var_type, disc_vars,  disc2var_idx, atol)
    if num_strong_var < length(reasonable_disc_vars)
        shuffle!(reasonable_disc_vars)
        reasonable_disc_vars = reasonable_disc_vars[1:num_strong_var]
    end

    # compute the gain for each reasonable candidate and choose the highest
    status, max_gain_var, left_node, right_node, gains, strong_restarts = branch_strong_on!(m,opts,step_obj,
        reasonable_disc_vars, disc2var_idx, opts.strong_restart, counter)

    step_obj.obj_gain += gains

    if status != :Resolve
        step_obj.l_nd = left_node
        step_obj.r_nd = right_node
    
        # set the variable to branch (best gain)
        node.state = :Done
        node.var_idx = max_gain_var
    end

    @assert max_gain_var != 0 || status == :PartlyInfeasible || status == :Infeasible || node.state == :Infeasible
    return status, max_gain_var, strong_restarts
end

function branch_reliable!(m,opts,step_obj,disc2var_idx,gains,counter) 
    idx = 0
    node = step_obj.node
    mu = opts.gain_mu
    reliability_param = opts.reliability_branching_threshold
    reliability_perc = opts.reliability_branching_perc
    num_strong_var = Int(round((reliability_perc/100)*m.num_disc_var))
    # if smaller than 2 it doesn't make sense
    num_strong_var = num_strong_var < 2 ? 2 : num_strong_var

    # use strong_branching_approx_time_limit to change num_strong_var
    if !isinf(opts.strong_branching_approx_time_limit)
        approx_time_per_node = 2*m.relaxation_time
        new_num_strong_var = Int(floor(opts.strong_branching_approx_time_limit/approx_time_per_node))
        new_num_strong_var = new_num_strong_var == 0 ? 1 : new_num_strong_var
        if new_num_strong_var < num_strong_var
            @warn "Changed num_strong_var to $new_num_strong_var because of strong_branching_approx_time_limit"
            num_strong_var = new_num_strong_var
        end
    end


    gmc_r = gains.minus_counter .< reliability_param
    gpc_r = gains.plus_counter  .< reliability_param

    strong_restarts = 0
    reasonable_disc_vars = []
    atol = opts.atol
    for i=1:length(gmc_r)
        if gmc_r[i] || gpc_r[i]
            idx = disc2var_idx[i]
            u_b = node.u_var[idx]
            l_b = node.l_var[idx]
            if isapprox(u_b,l_b; atol=atol) || is_type_correct(node.solution[idx],m.var_type[idx],atol)
                continue
            end
            push!(reasonable_disc_vars,i)
        end
    end
    if length(reasonable_disc_vars) > 0
        unrealiable_idx = sortperm(gains.minus_counter[reasonable_disc_vars])
        reasonable_disc_vars = reasonable_disc_vars[unrealiable_idx]
        num_reasonable = num_strong_var < length(reasonable_disc_vars) ? num_strong_var : length(reasonable_disc_vars)
        reasonable_disc_vars = reasonable_disc_vars[1:num_reasonable]
        
        status, max_gain_var, left_node, right_node, gains, strong_restarts = branch_strong_on!(m,opts,step_obj,
            reasonable_disc_vars, disc2var_idx, opts.strong_restart, counter)
        
        step_obj.obj_gain += gains
        
        step_obj.upd_gains = :GainsToTree
        new_gains = copy(step_obj.obj_gain)
        
        if status == :Infeasible
            return :Infeasible, 0, strong_restarts
        end
    else 
        step_obj.upd_gains = :GuessAndUpdate
        new_gains = copy(gains)
    end
    idx = branch_pseudo(m, node, disc2var_idx, new_gains, mu, atol)
    return :Normal, idx, strong_restarts
end

function branch_pseudo(m, node, disc2var_idx, obj_gain, mu, atol)
    # use the one with highest obj_gain which is currently continuous
    idx = 0
    scores, sort_idx = sorted_score_idx(node.solution, obj_gain, disc2var_idx, mu)
    for l_idx in sort_idx
        var_idx = disc2var_idx[l_idx]
        if !is_type_correct(node.solution[var_idx], m.var_type[var_idx], atol)
            u_b = node.u_var[var_idx]
            l_b = node.l_var[var_idx]
            # if the upper bound is the lower bound => should be type correct
            #= println("disc_var_idx: ", l_idx)
            println("node.solution[var_idx]: ", node.solution[var_idx])
            println("u_b: ", u_b)
            println("l_b: ", l_b) =#
            @assert !isapprox(u_b, l_b; atol=atol)
            idx = var_idx
            break
        end
    end
    return idx
end

function sorted_score_idx(x, gains, i2v, mu)
    g_minus, g_minus_c = gains.minus, gains.minus_counter
    g_plus, g_plus_c = gains.plus, gains.plus_counter
    g_minus_c += map(i -> (i == 0) && (i = 1), g_minus_c)
    g_plus_c += map(i -> (i == 0) && (i = 1), g_plus_c)
    scores = [score(f_minus(x[i2v[i]])*g_minus[i]/g_minus_c[i],f_plus(x[i2v[i]])*g_plus[i]/g_plus_c[i],mu) for i=1:length(g_minus)]
    sortedidx = sortperm(scores; rev=true)
    infsortedidx = sortperm(gains.inf_counter[sortedidx]; rev=true)
    return scores,sortedidx[infsortedidx]
end

"""
    Score function from 
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf
"""
function score(q_m, q_p, mu) 
    minq = q_m < q_p ? q_m : q_p
    maxq = q_m < q_p ? q_p : q_m
    return (1-mu)*minq+mu*maxq
end

function diff_obj(opts, node, cnode)
    if state_is_optimal(cnode.relaxation_state; allow_almost=opts.allow_almost_solved)
        return abs(node.best_bound - cnode.best_bound)
    else
        return Inf
    end
end

f_plus(x) = ceil(x)-x
f_minus(x) = x-floor(x)
sigma_plus(opts, node,r_nd,x) = diff_obj(opts, node,r_nd)/f_plus(x)
sigma_minus(opts, node,l_nd,x) = diff_obj(opts, node,l_nd)/f_minus(x)
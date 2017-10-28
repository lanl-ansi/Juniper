module BnBTree
import MINLPBnB
using JuMP
using Ipopt

rtol = 1e-6
atol = 1e-6
srand(1)
time_solve_leafs_get_idx = 0.0
time_solve_leafs_branch = 0.0

type BnBNode
    parent      :: Union{Void,BnBNode}
    idx         :: Int64
    level       :: Int64
    m           :: MINLPBnB.MINLPBnBModel
    var_idx     :: Int64
    left        :: Union{Void,BnBNode}
    right       :: Union{Void,BnBNode}
    state       :: Symbol
    hasbranchild :: Bool # has child where to branch or is :Branch
    best_bound  :: Union{Void,Float64}
end

type BnBTreeObj
    root        :: BnBNode
    incumbent   :: Union{Void,MINLPBnB.MINLPBnBModel}
    obj_gain    :: Vector{Float64} # gain of objective per variable
    obj_gain_c  :: Vector{Float64} # obj_gain / obj_gain_c => average gain
    int2var_idx :: Vector{Int64}
    var2int_idx :: Vector{Int64}
    options     :: MINLPBnB.SolverOptions
end

function init(m)
    node = BnBNode(nothing,1,1,m,0,nothing,nothing,:Branch,true,m.objval)
    obj_gain = zeros(m.num_int_bin_var)
    obj_gain_c = zeros(m.num_int_bin_var)
    int2var_idx = zeros(m.num_int_bin_var)
    var2int_idx = zeros(m.num_var)
    int_i = 1
    for i=1:m.num_var
        if m.var_type[i] != :Cont
            int2var_idx[int_i] = i
            var2int_idx[i] = int_i
            int_i += 1
        end
    end
    return BnBTreeObj(node,nothing,obj_gain,obj_gain_c,int2var_idx,var2int_idx,m.options)
end

function new_default_node(parent,idx,level,m;
                            var_idx=0,left=nothing,right=nothing,
                            state=:Solve,hasbranchild=true,best_bound=nothing)

    return BnBNode(parent,idx,level,m,var_idx,left,right,state,hasbranchild,best_bound)     
end

function check_print(vec::Vector{Symbol}, ps::Vector{Symbol})
    for v in vec
        if v in ps
            return true
        end
    end
    return false
end

"""
    branch_mostinfeasible(tree,node,num_var,var_type,x)

Get the index of an integer variable which is currently continuous which is most unintegral.
(nearest to *.5)
"""
function branch_mostinfeasible(tree,node,num_var,var_type,x)
    idx = 0
    max_diff = 0
    for i=1:num_var
        if var_type[i] != :Cont
            diff = abs(x[i]-round(x[i]))
            if diff > max_diff
                idx = i
                max_diff = diff
            end
        end
    end
    return idx
end

"""
    branch_strong((tree,node,num_var,var_type,x)

Try to branch on a few different variables and choose the one with highest obj_gain.
Update obj_gain for the variables tried and average the other ones.
"""
function branch_strong(tree,node,num_var,var_type,counter)
    # generate an of variables to branch on
    num_strong_var = tree.options.strong_branching_nvars

    # get reasonable candidates (not type correct and not already perfectly bounded)
    int_vars = tree.root.m.num_int_bin_var
    reasonable_int_vars = zeros(Int64,0)
    for i=1:int_vars
        idx = tree.int2var_idx[i]
        u_b = node.m.u_var[idx]
        l_b = node.m.l_var[idx]
        if isapprox(u_b,l_b,atol=atol) || BnBTree.is_type_correct(node.m.solution[idx],var_type[idx])
            continue
        end
        push!(reasonable_int_vars,i)
    end
    shuffle!(reasonable_int_vars)
    reasonable_int_vars = reasonable_int_vars[1:minimum([num_strong_var,length(reasonable_int_vars)])]

    # compute the gain for each reasonable candidate and choose the highest
    max_gain = 0.0
    max_gain_var = 0
    av_gain = 0.0
    strong_int_vars = zeros(Int64,0)
    left_node = nothing
    right_node = nothing
    
    for int_var_idx in reasonable_int_vars
        push!(strong_int_vars, int_var_idx)
        var_idx = tree.int2var_idx[int_var_idx]
        l_nd,r_nd = BnBTree.branch!(tree,node,var_idx;map_to_node=false)
        gain = BnBTree.compute_gain(node;l_nd=l_nd,r_nd=r_nd)
        if gain > max_gain
            max_gain = gain
            max_gain_var = var_idx
            left_node = l_nd
            right_node = r_nd
            if gain == Inf
                break
            end
        end
        av_gain += gain
        tree.obj_gain[int_var_idx] = gain
    end
    node.left = left_node
    node.right = right_node
    node.var_idx = max_gain_var
    av_gain /= int_vars
    rest = filter(i->!(i in strong_int_vars),1:int_vars)
    if counter == 1
        tree.obj_gain[rest] = av_gain
        tree.obj_gain_c += 1
    else
        tree.obj_gain_c[strong_int_vars] += 1
    end

    @assert max_gain_var != 0
    return max_gain_var
end

"""
    get_int_variable_idx(tree,node,num_var,var_type,x,counter=1)

Get the index of a variable to branch on.
"""
function get_int_variable_idx(tree,node,num_var,var_type,x,counter::Int64=1)    
    idx = 0
    branch_strat = tree.options.branch_strategy
    if branch_strat == :MostInfeasible
        return BnBTree.branch_mostinfeasible(tree,node,num_var,var_type,x)
    elseif branch_strat == :PseudoCost || branch_strat == :StrongPseudoCost
        if counter == 1 && branch_strat == :PseudoCost
            idx = BnBTree.branch_mostinfeasible(tree,node,num_var,var_type,x)
        elseif counter <= tree.options.strong_branching_nlevels && branch_strat == :StrongPseudoCost
            idx = BnBTree.branch_strong(tree,node,num_var,var_type,counter)
        else
            # use the one with highest obj_gain which is currently continous
            obj_gain_average = tree.obj_gain./tree.obj_gain_c
            sort_idx = tree.int2var_idx[sortperm(obj_gain_average, rev=true)]
            for l_idx in sort_idx
                if !is_type_correct(x[l_idx],var_type[l_idx])
                    u_b = node.m.u_var[l_idx]
                    l_b = node.m.l_var[l_idx]
                    # if the upper bound is the lower bound => no reason to branch
                    if isapprox(u_b,l_b,atol=atol)
                        continue
                    end
                    return l_idx
                end
            end
        end
    end
    @assert idx != 0
    return idx
end

"""
    is_type_correct(x,var_type)

Check whether a variable x has the correct type
"""
function is_type_correct(x,var_type)
    if var_type != :Cont
        if !isapprox(abs(round(x)-x),0, atol=atol, rtol=rtol)
           return false
        end
    end
    return true
end

"""
    are_type_correct(sol,types)

Check whether all variables have the correct type
"""
function are_type_correct(sol,types)
    for i=1:length(sol)
        if types[i] != :Cont
            if !isapprox(abs(round(sol[i])-sol[i]),0, atol=atol, rtol=rtol)
                return false
            end
        end
    end
    return true
end

"""
    solve_leaf(leaf)

Solve a leaf by relaxation leaf is just a node.
Set the state,hasbranchild and best_bound property
Return state
"""
function solve_leaf(leaf)
    status = JuMP.solve(leaf.m.model)
    leaf.m.objval   = getobjectivevalue(leaf.m.model)
    leaf.m.solution = getvalue(leaf.m.x)
    leaf.m.status = status
    if status == :Error
        # println(leaf.m.model)
        println(Ipopt.ApplicationReturnStatus[internalmodel(leaf.m.model).inner.status])
        # error("...")
        leaf.state = :Error
        leaf.hasbranchild = false
    elseif status == :Optimal
        # check if all int vars are int
        if BnBTree.are_type_correct(leaf.m.solution,leaf.m.var_type)
            leaf.state = :Integral
            leaf.hasbranchild = false
            leaf.best_bound = leaf.m.objval
        else
            leaf.state = :Branch
            leaf.best_bound = leaf.m.objval
        end
    else
        leaf.state = :Infeasible
        leaf.hasbranchild = false
    end
    return leaf.state
end

"""
    branch!(node::BnBNode,idx,ps)

Branch a node by using x[idx] <= floor(x[idx]) and x[idx] >= ceil(x[idx])
Solve both nodes and set current node state to done.
"""
function branch!(tree::BnBTreeObj,node::BnBNode,idx;map_to_node=true)
    global time_solve_leafs_get_idx, time_solve_leafs_branch
    ps = tree.options.log_levels
    l_m = Base.deepcopy(node.m)
    r_m = Base.deepcopy(node.m)

    # save that this node branches on this particular variable
    node.var_idx = idx

    l_x = l_m.x
    l_cx = l_m.solution[idx]
    r_x = r_m.x
    r_cx = r_m.solution[idx]
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("branch")
    
    if isapprox(l_m.u_var[idx],r_m.l_var[idx],atol=atol)
        error("Shouldn't solve again")
    end
    JuMP.setupperbound(l_x[idx], floor(l_cx))
    JuMP.setlowerbound(r_x[idx], ceil(r_cx))

    l_nd = BnBTree.new_default_node(node,node.idx*2,node.level+1,l_m)
    r_nd = BnBTree.new_default_node(node,node.idx*2+1,node.level+1,r_m)

    if map_to_node
        node.left = l_nd
        node.right = r_nd
    end
    node.state = :Done

    leaf_start = time()
    l_state = solve_leaf(l_nd)
    r_state = solve_leaf(r_nd)
    leaf_time = time()-leaf_start
    if map_to_node
        time_solve_leafs_branch += leaf_time
    else
        time_solve_leafs_get_idx += leaf_time
    end

    if BnBTree.check_print(ps,[:All])
        println("State of left leaf: ", l_state)
        println("State of right leaf: ", r_state)
        println("l sol: ", l_nd.m.solution)
        println("r sol: ", r_nd.m.solution)
    end
    return l_nd, r_nd
end

function compute_gain(node;l_nd::BnBNode=node.left,r_nd::BnBNode=node.right)
    gain = 0.0
    gc = 0
    frac_val = node.m.solution[node.var_idx]
    if l_nd.state == :Integral || r_nd.state == :Integral || l_nd.state == :Infeasible || r_nd.state == :Infeasible
        return Inf
    end
    if l_nd.state == :Error && r_nd.state == :Error
        return 0.0
    end
    if l_nd.state == :Branch || l_nd.state == :Integral
        int_val = floor(frac_val)
        gain += abs(node.best_bound-l_nd.best_bound)/abs(frac_val-int_val)
        gc += 1
    end
    if r_nd.state == :Branch || r_nd.state == :Integral
        int_val = ceil(frac_val)
        gain += abs(node.best_bound-r_nd.best_bound)/abs(frac_val-int_val)
        gc += 1
    end
    gc == 0 && return Inf
    gain /= gc
    return gain
end
"""
    update_gains(tree::BnBTreeObj,node::BnBNode,counter)

Update the objective gains for the branch variable used for node
"""
function update_gains(tree::BnBTreeObj,node::BnBNode,counter)
    gain = BnBTree.compute_gain(node)

    # update all (just average of the one branch we have)
    if counter == 1
        tree.obj_gain += gain
    else
        idx = tree.var2int_idx[node.var_idx]
        tree.obj_gain[idx] += gain
        tree.obj_gain_c[idx] += 1
    end
end

"""
    update_incumbent!(tree::BnBTreeObj,node::BnBNode)

Update the incumbent if there is a new Integral solution which is better.
Return true if updated false otherwise
"""
function update_incumbent!(tree::BnBTreeObj,node::BnBNode)
    ps = tree.options.log_levels
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("update_incumbent")

    l_nd = node.left
    r_nd = node.right
    l_state, r_state = l_nd.state, r_nd.state
    factor = 1
    if tree.root.m.obj_sense == :Min
        factor = -1
    end

    if l_state == :Integral || r_state == :Integral
        # both integral => get better
        if l_state == :Integral && r_state == :Integral
            if factor*l_nd.m.objval > factor*r_nd.m.objval
                possible_incumbent = l_nd.m
            else
                possible_incumbent = r_nd.m
            end
        elseif l_state == :Integral
            possible_incumbent = l_nd.m
        else
            possible_incumbent = r_nd.m
        end
        if tree.incumbent == nothing || factor*possible_incumbent.objval > factor*tree.incumbent.objval
            tree.incumbent = possible_incumbent
            return true
        end
    end
 
    return false
end

"""
    update_branch!(tree::BnBTreeObj,node::BnBNode)

Update the branch tree. If on both children can't be branched on
=> set hasbranchild = false and check the parents as well (bubble up)
If one of both children can be branched on => bubble up the best_bound
"""
function update_branch!(tree::BnBTreeObj,node::BnBNode)
    ps = tree.options.log_levels
    l_nd = node.left
    r_nd = node.right
    l_state, r_state = l_nd.state, r_nd.state
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("update branch")
    BnBTree.check_print(ps,[:All]) && println(l_state, " ", r_state)
    if l_state != :Branch && r_state != :Branch
        local_node = node
        local_node.hasbranchild = false

        # both children aren't branch nodes
        # bubble up to check where to set node.hasbranchild = false
        while local_node.parent != nothing
            local_node = local_node.parent
            BnBTree.check_print(ps,[:All]) && println("local_node.level: ", local_node.level)
            if local_node.left.hasbranchild || local_node.right.hasbranchild
                BnBTree.check_print(ps,[:All]) && println("break")
                break
            else
                local_node.hasbranchild = false
            end
        end
    end
    # Bubble up the best bound of the children
    # => The root has always the best bound of all of it's children
    factor = 1
    if tree.root.m.obj_sense == :Min
        factor = -1
    end

    while node != nothing
        BnBTree.check_print(ps,[:All]) && println("Node idx: ", node.idx)
        l_nd = node.left
        r_nd = node.right
        if l_nd.best_bound == nothing
            node.best_bound = r_nd.best_bound
        elseif r_nd.best_bound == nothing
            node.best_bound = l_nd.best_bound
        elseif factor*l_nd.best_bound > factor*r_nd.best_bound
            node.best_bound = l_nd.best_bound
        else
            node.best_bound = r_nd.best_bound
        end
        node = node.parent
    end
end

"""
    get_best_branch_node(tree::BnBTreeObj)

Get the index of the breach node which should be used for the next branch.
Currently get's the branch with the best best bound
"""
function get_best_branch_node(tree::BnBTreeObj)
    node = tree.root
    obj_sense = tree.root.m.obj_sense
    factor = 1
    if obj_sense == :Min
        factor = -1
    end

    if node.state == :Branch
        return node
    end

    last_best_bound = node.best_bound
    while true
        l_nd = node.left
        r_nd = node.right
        if node.best_bound != last_best_bound
            error("Best bound should be the same as the root bound")
        end
        if node.hasbranchild == true
            # get into best branch
            if l_nd.hasbranchild && r_nd.hasbranchild
                if factor*l_nd.best_bound > factor*r_nd.best_bound
                    node = l_nd
                else
                    node = r_nd
                end
            elseif !l_nd.hasbranchild && !r_nd.hasbranchild
                println("node idx: ", node.idx)
                print(tree)
                error("Infeasible")
            elseif l_nd.hasbranchild
                node = l_nd
            else
                node = r_nd
            end
        end
        if node.state == :Branch
            return node
        end
    end
end

"""
    prune!(node::BnBNode, value)

Get rid of nodes which have a worse best bound then specified by value. 
Is recursive
"""
function prune!(node::BnBNode, value)
    obj_sense = node.m.obj_sense
    factor = 1
    if obj_sense == :Min
        factor = -1
    end
    if node.hasbranchild && factor*value >= factor*node.best_bound
        node.hasbranchild = false 
        node.left = nothing
        node.right = nothing
    else
        if node.left != nothing
            prune!(node.left, value)
        end
        if node.right != nothing
            prune!(node.left, value)
        end
    end
end

"""
    prune!(tree::BnBTreeObj)

Call prune! for the root node using the incumbent value
"""
function prune!(tree::BnBTreeObj)
    incumbent_val = tree.incumbent.objval
    ps = tree.options.log_levels
    BnBTree.check_print(ps,[:All,:Incumbent]) && println("incumbent_val: ", incumbent_val)

    prune!(tree.root, incumbent_val)
end

function print(node::BnBNode,int2var_idx)
    indent = (node.level-1)*2
    indent_str = repeat(" ",indent)
    println(indent_str*"idx"*": "*string(node.idx))
    println(indent_str*"var_idx"*": "*string(node.var_idx))
    println(indent_str*"state"*": "*string(node.state))
    println(indent_str*"hasbranchild"*": "*string(node.hasbranchild))
    println(indent_str*"best_bound"*": "*string(node.best_bound))
    int_idx = zeros(Int,0)
    for i=1:node.m.num_int_bin_var
        push!(int_idx,int2var_idx[i])
    end
    
    println(indent_str*"u_var"*": "*string(node.m.u_var[int_idx]))
    println(indent_str*"l_var"*": "*string(node.m.l_var[int_idx]))
    return hcat(node.m.u_var[int_idx],node.m.l_var[int_idx])
end

function print_rec(node::BnBNode,int2var_idx;remove=false,bounds=[])
    if remove != :hasnobranchild || node.hasbranchild
        a = print(node,int2var_idx)
        push!(bounds,a)
        if node.left != nothing
            print_rec(node.left,int2var_idx;remove=remove,bounds=bounds)
        end
        if node.right != nothing
            print_rec(node.right,int2var_idx;remove=remove,bounds=bounds)
        end
    end
end

function print(tree::BnBTreeObj;remove=false)
    node = tree.root
    print_rec(node,tree.int2var_idx;remove=remove)
end

function print_table_header(fields, field_chars)
    ln = ""
    i = 1
    for f in fields
        padding = field_chars[i]-length(f)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= f
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        i += 1
    end
    println(ln)
    println(repeat("=", sum(field_chars)))
end

function is_table_diff(last_arr,new_arr)
    if length(last_arr) != length(new_arr)
        return true
    end    
    for i=1:length(last_arr)
        last_arr[i] != new_arr[i] && return true
    end
    return false 
end

function print_table(tree,start_time,fields,field_chars;last_arr=[])
    arr = []
    
    i = 1
    ln = ""
    for f in fields
        val = ""
        if f == "Incumbent"
            val = tree.incumbent != nothing ? string(round(tree.incumbent.objval,2)) : "-"
        elseif f == "Best Bound"
            val = string(round(tree.root.best_bound,2))
        elseif f == "Gap"
            if tree.incumbent != nothing
                b = tree.root.best_bound
                f = tree.incumbent.objval
                val = string(round(abs(b-f)/abs(f)*100,1))*"%"
            else
                val = "-"
            end
        elseif f == "Time"
            val = string(round(time()-start_time,1))
        end
        padding = field_chars[i]-length(val)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= val
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        push!(arr,val)
        i += 1
    end
    BnBTree.is_table_diff(last_arr[1:end-1],arr[1:end-1]) && println(ln)
    return arr
end

"""
    solve(tree::BnBTreeObj)

Solve the MIP part of a problem given by BnBTreeObj using branch and bound.
 - Identify the node to branch on
 - Get variable to branch on
 - Solve subproblems
"""
function solve(tree::BnBTreeObj)
    global time_solve_leafs_get_idx, time_solve_leafs_branch
    time_solve_leafs_get_idx = 0.0
    time_solve_leafs_branch = 0.0

    fields = ["Incumbent","Best Bound","Gap","Time"]
    field_chars = [28,28,7,8]
    
    if BnBTree.are_type_correct(tree.root.m.solution,tree.root.m.var_type)
        return tree.root.m
    end

    ps = tree.options.log_levels
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("Solve Tree")
    # get variable where to split
    node = tree.root
    counter = 1    

    branch_strat = tree.options.branch_strategy
    time_upd_gains = 0.0
    time_get_idx = 0.0
    time_branch = 0.0
    time_solve_leafs = 0.0
    
    print_table_header(fields,field_chars)

    time_bnb_solve_start = time()
    last_table_arr = []
    while true
        m = node.m
        get_idx_start = time()
        v_idx = BnBTree.get_int_variable_idx(tree,node,m.num_var,m.var_type,m.solution,counter)
        time_get_idx += time()-get_idx_start
    
        BnBTree.check_print(ps,[:All]) && println("v_idx: ", v_idx)

        branch_start = time()
        
        if node.left == nothing
            l_nd,r_nd = BnBTree.branch!(tree,node,v_idx)
        end
        time_branch += time()-branch_start

        if branch_strat == :PseudoCost || (branch_strat == :StrongPseudoCost && counter > tree.options.strong_branching_nlevels)
            upd_start = time()
            BnBTree.update_gains(tree,node,counter)    
            time_upd_gains += time()-upd_start
        end

        BnBTree.update_branch!(tree,node)

        # update incumbent
        if BnBTree.update_incumbent!(tree,node)
            BnBTree.check_print(ps,[:All]) && println("Prune")
            BnBTree.prune!(tree)
            BnBTree.check_print(ps,[:All]) && println("pruned")            
            BnBTree.check_print(ps,[:All]) && print(tree)
        end
        # check if best
        if tree.incumbent != nothing && tree.incumbent.objval == tree.root.best_bound
            break
        end
        if !tree.root.hasbranchild
            error("no child to branch on")
            break
        end
    
        # println("Best bound: ", tree.root.best_bound)
        # println("Node level: ", node.level)

        # if node.level == 3
            # print(tree)
            # error("t")
        # end
        # get best branch node
        node = BnBTree.get_best_branch_node(tree)
        if BnBTree.check_print(ps,[:Table]) 
            last_table_arr = print_table(tree,time_bnb_solve_start,fields,field_chars;last_arr=last_table_arr)
        end
        counter += 1
    end

    # print(tree)
    println("Incumbent status: ", tree.incumbent.status)

    time_bnb_solve = time()-time_bnb_solve_start
    println("#branches: ", counter)
    println("BnB time: ", round(time_bnb_solve,2))
    println("Solve leaf time get idx: ", round(time_solve_leafs_get_idx,2))
    println("Solve leaf time branch: ", round(time_solve_leafs_branch,2))
    println("Branch time: ", round(time_branch,2))
    println("Get idx time: ", round(time_get_idx,2))
    println("Upd gains time: ", round(time_upd_gains,2))
    return tree.incumbent
end

end
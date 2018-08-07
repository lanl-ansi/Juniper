type Aff
    sense     :: Symbol
    var_idx   :: Vector{Int64}
    coeff     :: Vector{Float64}
    rhs       :: Float64

    Aff() = new()
end

"""
    check_print(vec::Vector{Symbol}, ps::Vector{Symbol})

Check whether something should get printed. Only if the log_levels option
include the necessary level 
"""
function check_print(log_levels::Vector{Symbol}, necessary::Vector{Symbol})
    for v in log_levels
        if v in necessary
            return true
        end
    end
    return false
end

function generate_random_restart(m; cont=true)
    values = []
    for i=1:m.num_var
        lbi_def = true
        ubi_def = true
        if m.l_var[i] > typemin(Int64) 
            lbi = m.l_var[i]
        else
            lbi = typemin(Int64)
            lbi_def = false
        end

        if m.u_var[i] < typemax(Int64) 
            ubi = m.u_var[i]
        else
            ubi = typemin(Int64)
            ubi_def = false
        end

        if !ubi_def && !lbi_def 
            ubi = 10
            lbi = -10
        elseif !ubi_def
            ubi = lbi+20
        elseif !lbi_def
            lbi = ubi-20
        end             

        if m.var_type[i] == :Cont || cont
            push!(values,(ubi-lbi)*rand()+lbi)
        else
            lbi = Int(round(lbi))
            ubi = Int(round(ubi))
            push!(values, rand(lbi:ubi))
        end
    end
    return values
end

"""
    construct_affine_vector(m)

Construct a vector of affine expressions for all linear functions using the derivative
"""
function construct_affine_vector(m)
    js = MathProgBase.jac_structure(m.d)

    jg = zeros(length(js[1]))
    MathProgBase.eval_jac_g(m.d, jg, ones(m.num_var))

    # Construct the data structure for our affine constraints
    aff = Vector{Aff}(m.num_l_constr)
    for i=1:m.num_l_constr
        aff[i] = Aff()
        aff[i].var_idx = []
        aff[i].coeff = []
        constr_expr = MathProgBase.constr_expr(m.d,i)
        aff[i].rhs = constr_expr.args[3]
        aff[i].sense = constr_expr.args[1]
    end

    # if linear constraint the derivative are the coeffs
    idx = 1
    lconstr2constr = Vector{Int64}()
    constr2lconstr = Vector{Int64}()
    c = 1
    for i=1:m.num_constr
        if m.isconstrlinear[i]
            push!(lconstr2constr,i)
            push!(constr2lconstr,c)
            c+=1
        else
            push!(constr2lconstr,0)
        end
    end

    for row in js[1]
        if m.isconstrlinear[row]
            col = js[2][idx]
            aidx = constr2lconstr[row]
            push!(aff[row].var_idx, col)
            push!(aff[row].coeff, jg[idx])
        end
        idx += 1
    end
    return aff
end


""" 
    construct_complete_affine_matrix(m)

Construct full affine matrix by using m.affs and all variables 
use construct_disc_affine_matrix if only interested in discrete variables
"""
function construct_complete_affine_matrix(m)
    mat = zeros(length(m.affs),m.num_var)
    i = 1
    for aff in m.affs
        for part_idx in 1:length(aff.var_idx)
            var = aff.var_idx[part_idx]
            coeff = aff.coeff[part_idx]
            mat[i,var] = coeff
        end
        i += 1
    end
    return mat
end

""" 
    construct_disc_affine_matrix(m; only_non_zero=true)

Construct full affine matrix by using m.affs and only discrete variables 
use construct_complete_affine_matrix if interested in all variables
if only_non_zero is set to true all rows with all zeros are removed
"""
function construct_disc_affine_matrix(m; only_non_zero=true)
    mat = zeros(length(m.affs),m.num_disc_var)
    i = 1
    non_zero_idx = []
    for aff in m.affs
        non_zero = false
        for part_idx in 1:length(aff.var_idx)
            var = aff.var_idx[part_idx]
            if m.var_type[var] != :Cont
                coeff = aff.coeff[part_idx]
                int_var = m.var2disc_idx[var]
                mat[i,int_var] = coeff
                non_zero = true
            end
        end
        if non_zero
            push!(non_zero_idx,i)
        end
        i += 1
    end

    if only_non_zero
        return mat[non_zero_idx,:]
    end

    return mat
end

"""
    get_reasonable_int_vars(node, var_type, int_vars, disc2var_idx, atol)

Get all discrete variables which aren't close to discrete yet based on atol 
"""
function get_reasonable_int_vars(node, var_type, int_vars, disc2var_idx, atol)
    reasonable_int_vars = zeros(Int64,0)
    for i=1:int_vars
        idx = disc2var_idx[i]
        u_b = node.u_var[idx]
        l_b = node.l_var[idx]
        if isapprox(u_b,l_b; atol=atol) || is_type_correct(node.solution[idx],var_type[idx],atol)
            continue
        end
        push!(reasonable_int_vars,i)
    end
    return reasonable_int_vars
end
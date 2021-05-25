# MathOptInterface results
	
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.inner === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    else
        return model.inner.status
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    return string(model.inner.status)
end

function MOI.get(model::Optimizer, ::MOI.PrimalStatus)
    if state_is_optimal(model.inner.status; allow_almost=model.inner.options.allow_almost_solved) 
        return MOI.FEASIBLE_POINT
    else
        return MOI.INFEASIBLE_POINT
    end
end

function MOI.get(model::Optimizer, ::MOI.DualStatus)
    if state_is_optimal(model.inner.status; allow_almost=model.inner.options.allow_almost_solved) 
        return MOI.FEASIBLE_POINT
    else
        return MOI.INFEASIBLE_POINT
    end
end

	
function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
        @error "optimize! not called"
    end
    return model.inner.objval
end
	
function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
        @error "optimize! not called"
    end
    return model.inner.best_bound
end

function MOI.get(model::Optimizer, ::MOI.RelativeGap)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
        @error "optimize! not called"
    end
    if isnan(model.inner.objval) || isnan(model.inner.best_bound)
        return NaN
    end
    return abs(model.inner.best_bound-model.inner.objval)/abs(model.inner.objval)
end

function MOI.get(model::Optimizer, ::MOI.SolveTime)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
        @error "optimize! not called"
    end
    return model.inner.soltime
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
        @error "optimize! not called"
    end
    MOI.throw_if_not_valid(model, vi)
    return model.inner.solution[vi.value]
end

include("POD_experiment/blend029.jl")
include("POD_experiment/nous1.jl")

@testset "POD instances" begin

@testset "blend029 full strong branching" begin
    println("==================================")
    println("blend029 full strong branching")
    println("==================================")

    m,objval = get_blend029()

    JuMP.setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_perc = 100,
            strong_branching_nsteps = 100,
            strong_restart = true
    ))
    status = JuMP.solve(m)

    @test status == :Optimal

    juniper_val = JuMP.getobjectivevalue(m)
    best_bound_val = JuMP.getobjbound(m)
    gap_val = JuMP.getobjgap(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test isapprox(juniper_val, objval, atol=1e0)
    @test isapprox(best_bound_val, objval, atol=1e0)
    @test isapprox(gap_val, 0, atol=1e-2)
end


@testset "blend029 break strong branching time limit" begin
    println("==================================")
    println("blend029 full strong branching")
    println("==================================")

    m,objval = get_blend029()

    JuMP.setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_approx_time_limit=0.01,
            time_limit = 4,
            strong_restart = false
    ))
    status = JuMP.solve(m)

    @test status == :Optimal || status == :UserLimit

    juniper_val = JuMP.getobjectivevalue(m)
    best_bound_val = JuMP.getobjbound(m)
    gap_val = JuMP.getobjgap(m)

    # maximization problem
    @test best_bound_val >= juniper_val || isnan(juniper_val)
end

@testset "nous1 restart" begin
    println("==================================")
    println("nous1 restart")
    println("==================================")

    m = get_nous1()

    JuMP.setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_restart = true,
            mip_solver=Cbc.CbcSolver(logLevel=0)
    ))
    status = JuMP.solve(m)

    @test status == :Optimal
end

@testset "nous1 no restart" begin
    println("==================================")
    println("nous1 no restart")
    println("==================================")

    m = get_nous1()

    JuMP.setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_restart = false,
            mip_solver=Cbc.CbcSolver(logLevel=0)
    ))
    status = JuMP.solve(m)

    @test status == :Optimal
end



@testset "blend029 reliability" begin
    println("==================================")
    println("blend029 reliability")
    println("==================================")

    m,objval = get_blend029()

    JuMP.setsolver(m, DefaultTestSolver(
            branch_strategy=:Reliability,
            reliability_branching_perc = 50,
            reliability_branching_threshold = 5,
            strong_restart = true
    ))
    status = JuMP.solve(m)

    @test status == :Optimal

    juniper_val = JuMP.getobjectivevalue(m)
    best_bound_val = JuMP.getobjbound(m)
    gap_val = JuMP.getobjgap(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test isapprox(juniper_val, objval, atol=1e0)
    @test isapprox(best_bound_val, objval, atol=1e0)
    @test isapprox(gap_val, 0, atol=1e-2)
    end

end
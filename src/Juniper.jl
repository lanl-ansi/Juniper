module Juniper

import MathProgBase
import JuMP
import JSON

import LinearAlgebra: dot
import Random: shuffle!, seed!
import Distributed: @spawnat, myid, nprocs, nworkers, remotecall, remotecall_fetch
import Statistics: median


include("types.jl")

function Base.show(io::IO, opts::SolverOptions)
    longest_field_name = maximum([length(string(fname)) for fname in fieldnames(SolverOptions)])+2
    for name in fieldnames(SolverOptions)
        sname = string(name)
        pname = sname*repeat(" ", longest_field_name-length(sname))
        if getfield(opts,name) == nothing
            println(io, pname, ": NA")
        else
            println(io, pname, ": ", getfield(opts,name))
        end
    end
end

include("util.jl")
include("printing.jl")
include("solver.jl")
include("model.jl")
include("BnBTree.jl")


end

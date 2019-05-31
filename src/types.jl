# Types for Boundary Conditions
# Abstracts
abstract type BoundaryCondition end
abstract type HomogeneousBoundaryCondition <: BoundaryCondition end
abstract type NonhomogeneousBoundaryCondition <: BoundaryCondition end

# Types for Differences
abstract type DiscretizationMethod end
struct DefaultFirstDifference <: DiscretizationMethod end
struct ForwardFirstDifference <: DiscretizationMethod end
struct BackwardFirstDifference <: DiscretizationMethod end
struct CentralSecondDifference <: DiscretizationMethod end
struct JumpProcess{T} <: DiscretizationMethod where T <: AbstractArray
    jumps::T
end

"""
JumpProcess(x̄, jumps::AbstractArray, truncate = (:interior, :interior))

Returns a DiscretizationMethod object that can be used to construct 
a discretized operator jump process. 

`jumps` is a `(length(x̄)-2)`-vector whose ith element is an integer that represents 
a jump size in index and direction (by sign) from `i`th element of `interiornodes(x̄)` and
the first and second elements of `truncate` represent truncation location for 
the lower bound and upper bound when the jump is out of the truncated boundary 
(`:interior` for the first/last element of `interiornodes(̄x)` and `:exterior` for 
the first/last element of `x̄`). The default parameter is `(:interior, :interior)`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> JumpProcess(x̄, [-1; -1; -1; -1], (:interior, :interior))
JumpProcess{Array{Int64,1}}([0, -1, -1, -1])
```
"""
function JumpProcess(x̄::AbstractArray, jumps::AbstractArray, 
                    truncate = (:interior, :interior))
    # each ith element of jumps defines the jump direction (negative/positive) and size
    # where i ranges from 1 to M (interior points)
    @assert length(x̄) == length(jumps) + 2 
    M = length(x̄) - 2

    # define default destinations and truncation indices
    destinations = Array(1:M) + jumps
    cutoff_lb = (truncate[1] == :interior) ? 1 : 0 
    cutoff_ub = (truncate[2] == :interior) ? M : M+1

    # perform truncations
    destinations[destinations .< cutoff_lb] .= cutoff_lb
    destinations[destinations .> cutoff_ub] .= cutoff_ub
    jumps_after_truncation = Array{Int64}(destinations - Array(1:M))

    return JumpProcess(jumps_after_truncation)
end

"""
JumpProcess(x̄, uniform_jump::AbstractArray, truncate = (:interior, :interior))

Returns a DiscretizationMethod object that can be used to construct 
a discretized operator jump process. 

`uniform_jump` is a scalar Int64 that represents a jump size in index 
and direction (by sign) from all elements of `interiornodes(x̄)` and
the first and second elements of `truncate` represent truncation location for 
the lower bound and upper bound when the jump is out of the truncated boundary 
(`:interior` for the first/last element of `interiornodes(̄x)` and `:exterior` for 
the first/last element of `x̄`). The default parameter is `(:interior, :interior)`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> JumpProcess(x̄, -1)
JumpProcess{Array{Int64,1}}([0, -1, -1, -1])
```
"""
JumpProcess(x̄::AbstractArray, uniform_jump::Int64, truncate = (:interior, :interior)) = JumpProcess(x̄, uniform_jump*ones(length(x̄) - 2), truncate)

"""
JumpProcess(x̄, jumpf::Function, truncate = (:interior, :interior))

Returns a DiscretizationMethod object that can be used to construct 
a discretized operator jump process. 

`jumpf` is a function that takes an element of `interiornodes(x̄)` and returns 
a scalar that represents the corresponding nominal jump size and direction (by sign) and
the first and second elements of `truncate` represent truncation location for 
the lower bound and upper bound when the jump is out of the truncated boundary 
(`:interior` for the first/last element of `interiornodes(̄x)` and `:exterior` for 
the first/last element of `x̄`). The default parameter is `(:interior, :interior)`.
The code uses nearest-neighbour rule to determine the indices of destinations according
to `jumpf.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> jumpf(x) = -1.4
jumpf (generic function with 1 method)

julia> JumpProcess(x̄, jumpf)
JumpProcess{Array{Int64,1}}([0, -1, -1, -1])
```
"""
function JumpProcess(x̄::AbstractArray, jumpf::Function, 
    truncate = (:interior, :interior))
    # each ith element of jumps defines the jump direction (negative/positive) and size
    # where i ranges from 1 to M (interior points)
    x = interiornodes(x̄)
    M = length(x)

    # define default destinations and truncation indices
    destinations = x + jumpf.(x)

    cutoff_lb = (truncate[1] == :interior) ? x[1] : x̄[1] 
    cutoff_ub = (truncate[2] == :interior) ? x[end] : x̄[end]

    # perform truncations
    destinations[destinations .< cutoff_lb] .= cutoff_lb
    destinations[destinations .> cutoff_ub] .= cutoff_ub

    # find the corresponding destination indices and jump indices 
    destinations_indices_on_x̄ = (destination -> findnearestindex(x̄, destination)).(destinations)
    current_indices_on_x̄ = Array(2:(M+1))
    jumps_after_truncation = destinations_indices_on_x̄ - current_indices_on_x̄
    
    return JumpProcess(jumps_after_truncation)
end


"""
JumpProcess(x̄, uniform_jump_size::Real, truncate = (:interior, :interior))

Returns a DiscretizationMethod object that can be used to construct 
a discretized operator jump process. 

`uniform_jump_size` is a scalar in Real that represents a nominal jump size 
and direction (by sign) from all elements of `interiornodes(x̄)` and
the first and second elements of `truncate` represent truncation location for 
the lower bound and upper bound when the jump is out of the truncated boundary 
(`:interior` for the first/last element of `interiornodes(̄x)` and `:exterior` for 
the first/last element of `x̄`). The default parameter is `(:interior, :interior)`.
The code uses nearest-neighbour rule to determine the indices of destinations according
to `jumpf.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> JumpProcess(x̄, -1.4)
JumpProcess{Array{Int64,1}}([0, -1, -1, -1])
```
"""
JumpProcess(x̄::AbstractArray, uniform_jump_size::Real, truncate = (:interior, :interior)) = JumpProcess(x̄, (x -> uniform_jump_size), truncate)

# Concretes
struct Reflecting <: HomogeneousBoundaryCondition end
@with_kw struct Mixed{T} <: HomogeneousBoundaryCondition where T <: Real
    ξ::T = 0.0
    direction = :auto # :forward/:backward/:auto
end
@with_kw struct Absorbing{T} <: HomogeneousBoundaryCondition where T <: Int64
    loc::T = 0
end
@with_kw struct NonhomogeneousAbsorbing{T} <: NonhomogeneousBoundaryCondition where T <: Real 
    S::T  = 0.0
end
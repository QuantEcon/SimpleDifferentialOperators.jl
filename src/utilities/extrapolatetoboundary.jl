"""
    extrapolatetoboundary(v, x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})

Returns a `length(x̄)`-vector whose `2:(length(x̄)-1)` elements are `v`,
the first and last element are extrapolated `v` on the boundaries of `x̄` according to
boundary conditions `bc` given.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = -2:2
0:5

julia> x = interiornodes(x̄)
  1
  2
  3

julia> v = (x -> x^2).(x)
3-element Array{Int64,1}:
  1
  0
  1

julia> extrapolatetoboundary(v, x̄, (Absorbing(), Absorbing()))
5-element Array{Int64,1}:
  0
  1
  0
  1
  0

julia> extrapolatetoboundary(v, x̄, (Absorbing(), Reflecting()))
 5-element Array{Int64,1}:
  0
  1
  0
  1
  1

julia> extrapolatetoboundary(v, x̄, (Mixed(ξ = 3.0), Reflecting()))
5-element Array{Float64,1}:
 -0.5
  1.0
  0.0
  1.0
  1.0
  
```
"""
function extrapolatetoboundary(v, x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})
    # setup for extension
    T = eltype(v)

    # default extrapolated boundaries (absorbing for both sides)
    lb_extended = zero(T)
    ub_extended = zero(T)

    # extrapolate on the lower bound
    if (typeof(bc[1]) <: Reflecting)
        lb_extended = v[1]
    elseif (typeof(bc[1]) <: Mixed)
        ξ_lb = bc[1].ξ
        Δ_1_minus = x̄[2] - x̄[1]

        # apply extrapolation based on the directions
        if (bc[1].direction != :backward)
          lb_extended = v[1]/(1-ξ_lb*Δ_1_minus)
        else
          lb_extended = (1+ξ_lb*Δ_1_minus)*v[1]
        end
    end

    # extrapolate on the upper bound
    if (typeof(bc[2]) <: Reflecting)
        ub_extended = v[end]
    elseif (typeof(bc[2]) <: Mixed)
        ξ_ub = bc[2].ξ
        Δ_M_plus = x̄[end] - x̄[end-1]

        # apply extrapolation based on the directions
        if (bc[2].direction != :forward)
          ub_extended = v[end]/(1+ξ_ub*Δ_M_plus)
        else
          ub_extended = (1-ξ_ub*Δ_M_plus)*v[end]
        end
    end

    return [lb_extended; v; ub_extended]
end
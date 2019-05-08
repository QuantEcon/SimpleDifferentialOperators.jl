"""
    extrapolatetoboundary(v, x̄, bc::Tuple{Mixed, Mixed})

Returns a `length(x̄)`-vector whose `2:(length(x̄)-1)` elements are `v`,
the first and last element are extrapolated `v` on the boundaries of `x̄` according to
boundary conditions `bc` given.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> x = interiornodes(x̄)
1
2
3

julia> v = (x -> x^2).(x)
5-element Array{Int64,1}:
  1
  4
  9
 16
 25

julia> extrapolatetoboundary(v, x̄, (Mixed(ξ = 1), Mixed(ξ = 1)))
7-element Array{Int64,1}:
  1
  1
  4
  9
 16
 25
 25
```
"""
function extrapolatetoboundary(v, x̄, bc::Tuple{Mixed, Mixed})
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ

    Δ_1_minus = x̄[2] - x̄[1]
    Δ_M_plus = x̄[end] - x̄[end-1]

    return [v[1]/(1-ξ_lb*Δ_1_minus); v; v[end]/(1+ξ_ub*Δ_M_plus)]
end

"""
    extrapolatetoboundary(v, x̄, bc::Tuple{Reflecting, Reflecting})

Returns a `length(x̄)`-vector whose `2:(length(x̄)-1)` elements are `v`,
the first and last element are extrapolated `v` on the boundaries of `x̄` according to
boundary conditions `bc` given.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> x = interiornodes(x̄)
1
2
3

julia> v = (x -> x^2).(x)
5-element Array{Int64,1}:
  1
  4
  9
 16
 25

julia> extrapolatetoboundary(v, x̄, (Reflecting(), Reflecting()))
7-element Array{Int64,1}:
  1
  1
  4
  9
 16
 25
 25
```
"""
function extrapolatetoboundary(v, x̄, bc::Tuple{Reflecting, Reflecting})
    return [v[1]; v; v[end]]
end

"""
    extrapolatetoboundary(v, x̄, bc::Tuple{Absorbing, Absorbing})

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
```
"""
function extrapolatetoboundary(v, x̄, bc::Tuple{Absorbing, Absorbing})
    T = eltype(v)
    return [zero(T); v; zero(T)]
end


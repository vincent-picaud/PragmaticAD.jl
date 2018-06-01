#+Tape 
struct ∂_j{T<:Real}
    value::T
    j::Int 
end
#+Tape 
struct Tape{T<:AbstractFloat}
    i_offset::Vector{Int}
    dϕ::Vector{∂_j{T}}
end
#+Tape 
function Tape{T<:AbstractFloat}(::Type{T},tape_initial_size::Int = 100)
    t=Tape{T}(Vector{Int}(tape_initial_size),Vector{∂_j{T}}(tape_initial_size))
    resize!(t.i_offset,1)
    t.i_offset[1]=1
    resize!(t.dϕ,0)
    return t
end

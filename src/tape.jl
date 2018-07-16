# +API
export afloat_count

# +Tape
# Used to store partial derivatives $\partial_j$
struct ∂_j{T<:Real}
    value::T
    j::Int 
end

# +Tape
#
Base.show(io::IO, x::∂_j) = print(io,"\partial_{v$(x.j)}ϕ=$(x.value)")

# +Tape, Internal L:Tape 
#
# Used to store $d\Phi^{(m)}$ using a storage scheme close to the
# Compressed Row Storage (CRS)
#
# *Note:* we use =mutable= to have a "by reference" semantic to avoid
#  having a tape copy (which would be catastrophic).
#
mutable struct Tape{T<:AbstractFloat}
    i_offset::Vector{Int}
    dϕ::Vector{∂_j{T}}
end

# +Tape, Internal
#
# Tape constructor
# !PragmaticAD.Tape(Float64)
#
function Tape{T<:AbstractFloat}(::Type{T},tape_initial_size::Int = 100)
    t=Tape{T}(Vector{Int}(tape_initial_size),Vector{∂_j{T}}(tape_initial_size))
    resize!(t.i_offset,1)
    t.i_offset[1]=1
    resize!(t.dϕ,0)
    return t
end




# +Tape
#
# Returns a new [[AFloat][]] index. Generally used at creation time,
# see [[create_tape_record_value][]]
#
afloat_next_index{T}(tape::Tape{T})::Int = length(tape.i_offset)

# +Tape,API
#
# Returns how many [[AFloat][]] are stored in the current [[Tape][]]
#
afloat_count{T}(tape::Tape{T})::Int = afloat_next_index(tape)-1

# +Tape
#
# A helper function that push a NTuple into an array
#
function Base.push!{T,N}(v::Array{T,1},x::NTuple{N,T})::Array{T,1}
    const n = length(v)
    resize!(v,n+N)
    for i in 1:N
        v[n+i]=x[i]
    end
    return v
end

# +Tape
#
# Creates a new differential $d\phi = \sum \partial_j \phi d x_j$ and
# returns the associated [[AFloat][]].
#
function create_tape_record{T,N}(tape::Tape{T},value::T,dϕ::NTuple{N,∂_j{T}})::AFloat{T} 
    const afloat_index = afloat_next_index(tape)
    const afloat=AFloat{T}(value,afloat_index)
    push!(tape.i_offset,tape.i_offset[end]+N)
    push!(tape.dϕ,dϕ)
    return afloat
end

# +Tape L:create_tape_record_value
#
# Creates a new [[AFloat][]]
#
function create_tape_record{T}(tape::Tape{T},value::T)::AFloat{T} 
    const afloat_index = afloat_next_index(tape)
    return create_tape_record(tape,value,(∂_j{T}(T(1),afloat_index),))
end

# +Tape L:tape_position
#
# Returns tape position.
#
# Usage example:
#
# #+BEGIN_SRC julia :eval never
# tape=PragmaticAD.get_tape(Float64)
# tpos=PragmaticAD.tape_position(tape)
# # some computations
# rewind_tape!(tape,tpos)
# #+END_SRC
#
# See: [[rewind_tape][]]
tape_position{T}(tape::Tape{T}) = afloat_next_index(tape)

# +Tape L:rewind_tape
#
# Rewinds tape
#
# See: [[tape_position][]]
#
function rewind_tape!{T}(tape::Tape{T},tape_position::Int) 
    resize!(tape.i_offset,tape_position)
    resize!(tape.dϕ,tape.i_offset[end]-1)
    return nothing
end

# +Tape                                   L:tape_f_gradient
#
# Computes differential adjoint vector action (reverse mode)
#
function f_gradient{T}(tape::Tape{T},
                       index::Int,
                       stop_at_tape_position::Int)
    @inline global2local_idx(i_glob) = i_glob - stop_at_tape_position + 1 
    const n_global = afloat_count(tape)
    const n_local = global2local_idx(n_global)
    grad = zeros(T,n_local)
    grad[global2local_idx(index)]=1
    for i in n_global:-1:stop_at_tape_position
        #    println(" $i $(global2local_idx(i))")
        const grad_i = grad[global2local_idx(i)]
        if grad_i != T(0)
            grad[global2local_idx(i)] = T(0)
            for j in tape.i_offset[i]:tape.i_offset[i+1]-1
                grad[global2local_idx(tape.dϕ[j].j)] += grad_i*tape.dϕ[j].value
            end
        end
    end

    return grad
end 



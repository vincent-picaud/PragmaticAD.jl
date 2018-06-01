# +Tape, Internal
# Used to store partial derivatives $\partial_j$
struct ∂_j{T<:Real}
    value::T
    j::Int 
end

# +Tape, Internal
#
# Used to store $d\Phi^{(m)}$ using a storage scheme close to the
# Compressed Row Storage (CRS)
#
struct Tape{T<:AbstractFloat}
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
struct AFloat{T<:AbstractFloat} <: Real
    value::T
    j::Int 
end

const AFloat32 = AFloat{Float32}
tape_Float32 = Tape(Float32)
get_tape(::Type{Float32}) = tape_Float32
const AFloat64 = AFloat{Float64}
tape_Float64 = Tape(Float64)
get_tape(::Type{Float64}) = tape_Float64

#+Tape
# 
afloat_next_index{T}(tape::Tape{T})::Int = length(tape.i_offset)
afloat_count{T}(tape::Tape{T})::Int = afloat_next_index(tape)-1

function Base.push!{T,N}(v::Array{T,1},x::NTuple{N,T})::Array{T,1}
    const n = length(v)
    resize!(v,n+N)
    for i in 1:N
        v[n+i]=x[i]
    end
    return v
end

function create_tape_record{T,N}(tape::Tape{T},value::T,dϕ::NTuple{N,∂_j{T}})::AFloat{T} 
    const afloat_index = afloat_next_index(tape)
    const afloat=AFloat{T}(value,afloat_index)
    push!(tape.i_offset,tape.i_offset[end]+N)
    push!(tape.dϕ,dϕ)
    return afloat
end
function create_tape_record{T}(tape::Tape{T},value::T)::AFloat{T} 
    const afloat_index = afloat_next_index(tape)
    return create_tape_record(tape,value,(∂_j{T}(T(1),afloat_index),))
end
AFloat{T}(value::S) where {T<:AbstractFloat,S<:Real} = create_tape_record(get_tape(T),convert(T,value))
tape_position{T}(tape::Tape{T}) = afloat_next_index(tape)
function rewind_tape!{T}(tape::Tape{T},tape_position::Int) 
    resize!(tape.i_offset,tape_position)
    resize!(tape.dϕ,tape.i_offset[end]-1)
    return nothing
end 
function Base.sin{T}(x::AFloat{T})::AFloat{T}
    const value = sin(x.value)
    const dϕ  = (∂_j{T}(cos(x.value),x.j),)
    return create_tape_record(get_tape(T),value,dϕ)
end 
function Base.:+{T}(x::AFloat{T},y::AFloat{T})::AFloat{T}
    const value = x.value+y.value
    const dϕ  = (∂_j{T}(T(1),x.j),∂_j{T}(T(1),y.j))
    return create_tape_record(get_tape(T),value,dϕ)
end
function Base.:*{T}(x::AFloat{T},y::AFloat{T})::AFloat{T}
    const value = x.value*y.value
    const dϕ  = (∂_j{T}(y.value,x.j),∂_j{T}(x.value,y.j))
    return create_tape_record(get_tape(T),value,dϕ)
end
function f_gradient{T}(tape::Tape{T},index::Int)
    const n = afloat_count(tape)
    grad = zeros(T,n)
    grad[index]=1
    for i in n:-1:1
        const grad_i = grad[i]
        if grad_i != T(0)
            grad[i] = T(0)
            for j in tape.i_offset[i]:tape.i_offset[i+1]-1
                grad[tape.dϕ[j].j] += grad_i*tape.dϕ[j].value
            end
        end
    end

    return grad
end 
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
f_gradient{T}(y::AFloat{T})=f_gradient(get_tape(T),y.j)
f_gradient{T}(y::AFloat{T},stop_at_tape_position::Int)=f_gradient(get_tape(T),y.j,stop_at_tape_position)

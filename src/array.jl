export AArray

# +Internal
import Base: setindex!,getindex,size,IndexStyle
    
# +AArray, API L:AArray
#
# This is the array type to use in place of built-in Julia =Array{T,N}= types.
#
# *Parameters*:
#
# - =AIDX= is an extra parameter defining how [[AFloat][]] indices are
#   stored. This type can be different from =P=. For instance in case
#   of sparse matrix, only an uni-dimensional array to store indices
#   associated to non-zero components is required, no need to the
#   store sparsity pattern twice. 
#
# *Design*:
# - [[id:c39435d9-532c-4390-8089-bde4c5e53f3e][Multiple parameters in the struct definition]]
# 
struct AArray{AT<:AFloat,N,P<:AbstractArray,AIDX}  <: AbstractArray{AT,N}
    parent::P
    aidx::AIDX

    # +AArray, API
    #
    # Construction from an "usual" Julia array of real numbers.
    #
    # Usage example:
    #
    # !AArray(rand(2,1))
    #
    # *Design*:
    # - [[id:b384f347-8c27-42ae-9759-2914d67cad4d][Restricted set of constructors]]
    #
    function AArray(p::Array{T,N}) where {T<:AbstractFloat,N} 
        aidx=similar(p,Int)
        const n_chunk = length(p)
        i_start=create_tape_chunk(get_tape(T),n_chunk)
        # broacast trick
        aidx[:].=i_start:i_start+n_chunk-1
        return new{AFloat{T},N,typeof(p),typeof(aidx)}(p,aidx)
    end 
end

# +AArray, Internal
#
# This function increase array size by =positive_integer= and returns
# the initial size
#
function increase_size!(v::Array{T,1},positive_integer::Int)::Int where {T}
    n=length(v)
    resize!(v,n+positive_integer)
    return n
end 

# +AArray, Internal
#
# This function allocate a chunk of =chunk_n= new [[AFloat][]].
# It returns the index of the first one.
#
function create_tape_chunk(tape::Tape{T},chunk_n::Int) where {T}
    const init_offset_n = increase_size!(tape.i_offset,chunk_n)
    const init_dϕ_n = increase_size!(tape.dϕ,chunk_n)
    
    for i in 1:chunk_n
        j=tape.i_offset[i+init_offset_n-1]
        tape.dϕ[i+init_dϕ_n]=∂_j{T}(T(1),j)
        tape.i_offset[i+init_offset_n]=1+j
    end 
    
    return init_offset_n
end



# +AArray, Internal
#
# The recommended way to access aidx member *type*
#
# See [[id:c39435d9-532c-4390-8089-bde4c5e53f3e][Multiple parameters in the struct definition]]
#
aidx_type(::Type{AArray{AT,N,P,AIDX}}) where {AT,N,P,AIDX} = AIDX

# +AArray, Internal
#
# The recommended way to access aidx member
#
# See [[id:c39435d9-532c-4390-8089-bde4c5e53f3e][Multiple parameters in the struct definition]]
#
aidx(aa::AArray) = aa.aidx

# +AArray, Internal
#
# The recommended way to access parent member *type*
#
# See [[id:c39435d9-532c-4390-8089-bde4c5e53f3e][Multiple parameters in the struct definition]]
#
parent_type(::Type{AArray{AT,N,P,AIDX}}) where {AT,N,P,AIDX} = P

# +AArray, Internal
#
# The recommended way to access parent membe
#
# See [[id:c39435d9-532c-4390-8089-bde4c5e53f3e][Multiple parameters in the struct definition]]
#
parent(aa::AArray) = aa.parent




# +AArray
#
# *Design*:
# - [[id:1b16ffef-47fa-473c-b033-de4a864dcaf3][Interface, main methods to redefine]]
#
size(aa::AArray) = size(parent(aa))

# +AArray
#
# *Design*:
# - [[id:1b16ffef-47fa-473c-b033-de4a864dcaf3][Interface, main methods to redefine]]
#
IndexStyle(::Type{<:AArray{AT,N,P,AIDX}}) where {AT,N,P,AIDX} = IndexStyle(P)

# +AArray
#
# *Design*:
# - [[id:1b16ffef-47fa-473c-b033-de4a864dcaf3][Interface, main methods to redefine]]
#
getindex(aa::AArray{AT}, i::Int) where {AT} = AT(getindex(parent(aa),i),getindex(aidx(aa),i))


# +AArray
#
# *Design*:
# - [[id:1b16ffef-47fa-473c-b033-de4a864dcaf3][Interface, main methods to redefine]]
# - [[id:ffb7408e-cb7b-4c7e-a1f8-39c3af3d37d5][Remark]]
#
@inline getindex(aa::AArray{AT,N,P}, I::Vararg{Int, N}) where {AT,N,P} = AT(getindex(parent(aa),I...),getindex(aidx(aa),I...))

# +AArray
#
# *Design*:
# - [[id:1b16ffef-47fa-473c-b033-de4a864dcaf3][Interface, main methods to redefine]]
#
@inline function setindex!(aa::AArray{AT,N}, v::AT, I::Vararg{Int, N}) where {AT,N} 
    checkbounds(parent(aa), i)
    # @inbound
    setindex!(aidx(aa),v.j,I...)
    setindex!(parent(aa),v.value,I...)
    return v
end

# +AArray
#
# *Design*:
# - [[id:1b16ffef-47fa-473c-b033-de4a864dcaf3][Interface, main methods to redefine]]
#
@inline function setindex!(aa::AArray{AT,N}, v, I::Vararg{Int, N}) where {AT,N} 
    av=AT(v)
    setindex!(aa,av,I)
    return av
end 



# function similar(aa::AArray{T,N,P}, ::Type{TS}, dims::Dims) where {T,TS,N,P}
#     return similar(parent(A), T, dims)
# end


Base.convert(::Type{AFloat{T}}, x::Integer) where {T<:AbstractFloat} = AFloat{T}(x)
Base.promote_rule(::Type{AFloat{T}}, ::Type{S}) where {T<:AbstractFloat,S<:AbstractFloat} = AFloat{T}
Base.promote_rule(::Type{AFloat{T}}, ::Type{S}) where {T<:AbstractFloat,S<:Integer} =  AFloat{T}


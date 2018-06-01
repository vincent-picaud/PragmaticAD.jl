struct AArray{AT<:AFloat,N,P<:AbstractArray,AIDX}  <: AbstractArray{AT,N}
    parent::P
    aidx::AIDX

    function AArray(p::Array{T,N}) where {T<:AbstractFloat,N} 
        aidx=similar(p,Int)
        const n_chunk = length(p)
        i_start=create_tape_chunk(get_tape(T),n_chunk)
        # broacast trick
        aidx[:].=i_start:i_start+n_chunk-1
        return new{AFloat{T},N,typeof(p),typeof(aidx)}(p,aidx)
    end 
end
# return initial size 
function increase_size!(v::Array{T,1},positive_integer::Int)::Int where {T}
    n=length(v)
    resize!(v,n+positive_integer)
    return n
end 

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
aidx_type(::Type{AArray{AT,N,P,AIDX}}) where {AT,N,P,AIDX} = AIDX
parent_type(::Type{AArray{AT,N,P,AIDX}}) where {AT,N,P,AIDX} = P
parent(aa::AArray) = aa.parent
aidx(aa::AArray) = aa.aidx
Base.size(aa::AArray) = size(parent(aa))
Base.IndexStyle(::Type{<:AArray{AT,N,P,AIDX}}) where {AT,N,P,AIDX} = IndexStyle(P)
Base.getindex(aa::AArray{AT}, i::Int) where {AT} = AT(getindex(parent(aa),i),getindex(aidx(aa),i))
# CAVEAT: my first attempt was to write getindex(parent(aa),I)
# ERROR: ArgumentError: invalid index: (1, 1)
# Stacktrace:
#  [1] getindex(::Array{Float64,2}, ::Tuple{Int64,Int64}) at ./abstractarray.jl:883
@inline Base.getindex(aa::AArray{AT,N,P}, I::Vararg{Int, N}) where {AT,N,P} = AT(getindex(parent(aa),I...),getindex(aidx(aa),I...))

@inline function Base.setindex!(aa::AArray{AT,N}, v::AT, I::Vararg{Int, N}) where {AT,N} 
    checkbounds(parent(aa), i)
    # @inbound
    setindex!(aidx(aa),v.j,I...)
    setindex!(parent(aa),v.value,I...)
    return v
end

@inline function Base.setindex!(aa::AArray{AT,N}, v, I::Vararg{Int, N}) where {AT,N} 
    av=AT(v)
    setindex!(aa,av,I)
    return av
end 
# function Base.similar(aa::AArray{T,N,P}, ::Type{TS}, dims::Dims) where {T,TS,N,P}
#     return similar(parent(A), T, dims)
# end
Base.convert(::Type{AFloat{T}}, x::Integer) where {T<:AbstractFloat} = AFloat{T}(x)
Base.promote_rule(::Type{AFloat{T}}, ::Type{S}) where {T<:AbstractFloat,S<:AbstractFloat} = AFloat{T}
Base.promote_rule(::Type{AFloat{T}}, ::Type{S}) where {T<:AbstractFloat,S<:Integer} =  AFloat{T}


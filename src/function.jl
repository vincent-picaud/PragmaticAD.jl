#
# COMPARISON OPERATORS 
#
# See: [[id:66703b8c-0f31-49db-bf41-268758ac27a9][Specialization of common functions]]
#
# Also: to be compared with ugly C++ macros.
#
for op = (:(==), :(!=), :(<), :(>), :(<=), :(>=))
    @eval begin
        import Base: ($op)
        ($op)(x1::AFloat{T},x2::AFloat{T}) where {T} = ($op)(x1.value,x2.value)
        ($op)(x1::AFloat{T},x2::Real) where {T} = ($op)(x1.value,x2)
        ($op)(x1::Real,x2::AFloat{T}) where {T} = ($op)(x1,x2.value)
    end
end


###################
# Unary functions #
###################


function Base.sin{T}(x::AFloat{T})::AFloat{T}
    const value = sin(x.value)
    const dϕ  = (∂_j{T}(cos(x.value),x.j),)
    return create_tape_record(get_tape(T),value,dϕ)
end


####################
# Binary operators #
####################

function Base.:+{T}(x::AFloat{T},y::AFloat{T})::AFloat{T}
    const value = x.value+y.value
    const dϕ  = (∂_j{T}(T(1),x.j),∂_j{T}(T(1),y.j))
    return create_tape_record(get_tape(T),value,dϕ)
end

# +API
# Overload op: (* AF AF)
#
function Base.:*{T}(x::AFloat{T},y::AFloat{T})::AFloat{T}
    const value = x.value*y.value
    const dϕ  = (∂_j{T}(y.value,x.j),∂_j{T}(x.value,y.j))
    return create_tape_record(get_tape(T),value,dϕ)
end
# +API
# Overload op: (* F AF)
#
function Base.:*{T}(x::Real,y::AFloat{T})::AFloat{T}
    const value = x*y.value
    const dϕ  = (∂_j{T}(x,y.j),)
    return create_tape_record(get_tape(T),value,dϕ)
end
# +API
# Overload op: (* AF F)
#
function Base.:*{T}(x::AFloat{T},y::Real)::AFloat{T}
    const value = x.value*y
    const dϕ  = (∂_j{T}(y,x.j),)
    return create_tape_record(get_tape(T),value,dϕ)
end

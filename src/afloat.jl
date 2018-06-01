# +API
export AFloat

# +AFloat,API L:AFloat
#
# Specialization of number alllowing to track and record operations in the [[Tape][]]
#
# - TODO [ ] support complex numbers (use [[https://en.wikipedia.org/wiki/Wirtinger_derivatives][Wirtinger_derivatives]]).
struct AFloat{T<:AbstractFloat} <: Number
    value::T
    j::Int 
end

# +AFloat,API
#
# Creates a new [[AFloat][]] from its value.
#
# !AFloat(-2.0)
#
AFloat{T}(value::S) where {T<:AbstractFloat,S<:Real} = create_tape_record(get_tape(T),convert(T,value))
# +AFloat,API
AFloat(value::T) where {T<:AbstractFloat} = AFloat{T}(value)

# +AFloat 
#
# Defines pretty print as explained [[https://docs.julialang.org/en/latest/manual/types/#man-custom-pretty-printing-1][in the official doc]].
Base.show(io::IO, x::AFloat) = print(io,x.value,"_",x.j)

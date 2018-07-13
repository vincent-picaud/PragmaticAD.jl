# This file defines function wrapper
#
export f_∇f



# +Tape
#
# Convenience function that call [[tape_f_gradient][]]
f_gradient{T}(y::AFloat{T},stop_at_tape_position::Int)=f_gradient(get_tape(T),y.j,stop_at_tape_position)

# +Internal
# Transforms regular input into autodiff 
# - Array -> AArray
# - Float -> AFloat
_f_∇f_input_float_to_afloat(x::Array) = AArray(x)
# +Internal
_f_∇f_input_float_to_afloat(x::Real) = AFloat(x)

# +Internal
# Transform local gradient into gradient (output)
function _f_∇f_output_grad(ax::AArray,grad_ay::Array,tpos::Int)
    @assert size(ax)==size(grad_ay)
    const grad_y=map(ax_i->grad_ay[ax_i.j-tpos+1],ax)
    return grad_y
end 
# +Internal
function _f_∇f_output_grad(ax::AFloat,grad_ay::Array,tpos::Int)
    return grad_ay[ax.j-tpos+1]
end

# +API
#
# Wraps a function to compute its value and gradient
#
# *Note*: inspired from [[https://github.com/denizyuret/AutoGrad.jl/blob/master/src/core.jl][AutoGrad.jl/core.jl]]
#
# Scalar example:
#
# !f(x::Number,y::Number) = y*sin(x)
# !af=PragmaticAD.f_∇f(f,2)
# !af(2,5.)
#
# 
function f_∇f(f::Function, argnum::Int=1) 

    function local_f_∇f(args...; kwargs...)

        arg_wrt = args[argnum]
        T = eltype(arg_wrt)
        
        #  check_tape_invariant(get_tape(T))

        tpos = tape_position(get_tape(T))
        args = Any[args...] # to make args writable

        args[argnum] = _f_∇f_input_float_to_afloat(arg_wrt)

        ay=f(args...,kwargs...)

        grad_ay=f_gradient(ay,tpos)
        y=ay.value
        grad_y=_f_∇f_output_grad(args[argnum],grad_ay,tpos)
        
        rewind_tape!(get_tape(T),tpos)
        
        return y,grad_y
    end
    return local_f_∇f
end



# +Tape
#
# Convenience function that call [[tape_f_gradient][]]

f_gradient{T}(y::AFloat{T},stop_at_tape_position::Int)=f_gradient(get_tape(T),y.j,stop_at_tape_position)

#
# Define global tapes
#

# +AFloat, API
const AFloat32 = AFloat{Float32}
tape_Float32 = Tape(Float32)
get_tape(::Type{Float32}) = tape_Float32
# +AFloat, API
const AFloat64 = AFloat{Float64}
tape_Float64 = Tape(Float64)
get_tape(::Type{Float64}) = tape_Float64

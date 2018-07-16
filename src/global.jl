#
# Define global tapes
#

# +AFloat, API
const AFloat32 = AFloat{Float32}
tape_AFloat32 = Tape(Float32)
get_tape(::Type{Float32}) = tape_AFloat32

# +AFloat, API
const AFloat64 = AFloat{Float64}
tape_AFloat64 = Tape(Float64)
get_tape(::Type{Float64}) = tape_AFloat64

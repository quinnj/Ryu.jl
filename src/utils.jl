uint(x::Float16) = Core.bitcast(UInt16, x)
uint(x::Float32) = Core.bitcast(UInt32, x)
uint(x::Float64) = Core.bitcast(UInt64, x)

mantissabits(::Type{Float16}) = 10
mantissabits(::Type{Float32}) = 23
mantissabits(::Type{Float64}) = 52

exponentbits(::Type{Float16}) = 5
exponentbits(::Type{Float32}) = 8
exponentbits(::Type{Float64}) = 11

bias(::Type{Float16}) = 15
bias(::Type{Float32}) = 127
bias(::Type{Float64}) = 1023

pow5_bitcount(::Type{Float16}) = 0 # TODO!!
pow5_bitcount(::Type{Float32}) = 61
pow5_bitcount(::Type{Float64}) = 121

pow5_inv_bitcount(::Type{Float16}) = 0 # TODO!!
pow5_inv_bitcount(::Type{Float32}) = 59
pow5_inv_bitcount(::Type{Float64}) = 122

qinvbound(::Type{Float16}) = 0 # TODO!!
qinvbound(::Type{Float32}) = 9
qinvbound(::Type{Float64}) = 21

qbound(::Type{Float16}) = 0 # TODO!!
qbound(::Type{Float32}) = 31
qbound(::Type{Float64}) = 63

pow5invsplit(::Type{Float16}, q) = 0 # TODO!!
pow5invsplit(::Type{Float32}, q) = FLOAT_POW5_INV_SPLIT[q + 1]
pow5invsplit(::Type{Float64}, q) = DOUBLE_POW5_INV_SPLIT[q + 1]

pow5split(::Type{Float16}, i) = 0 # TODO!!
pow5split(::Type{Float32}, i) = FLOAT_POW5_SPLIT[i + 1]
pow5split(::Type{Float64}, i) = DOUBLE_POW5_SPLIT[i + 1]

# Returns floor(log_10(2^e)).
@inline function log10pow2(e)
    # The first value this approximation fails for is 2^1651 which is just greater than 10^297.
    return (Core.bitcast(UInt32, e) * UInt32(78913)) >> 18
end

@inline function log10pow5(e)
    return ((e % UInt32) * UInt32(732923)) >> 20
end

# Returns e == 0 ? 1 : ceil(log_2(5^e)).
@inline function pow5bits(e::Int32)::Int32
    # This approximation works up to the point that the multiplication overflows at e = 3529.
    # If the multiplication were done in 64 bits, it would fail at 5^4004 which is just greater
    # than 2^9297.
    return Core.bitcast(Int32, ((Core.bitcast(UInt32, e) * UInt32(1217359)) >> 19) + UInt32(1))
end

@inline function pow5factor(value::T) where {T}
    count = 0
    while true
        q = div(value, T(5))
        r = value - T(5) * q
        r != 0 && break
        value = q
        count += 1
    end
    return count
end

@inline function multipleOfPowerOf5(value, p)
    return pow5factor(value) >= p
end

@inline function multipleOfPowerOf2(value::T, p) where {T}
    return (value & ((T(1) << p) - 1)) == 0
end

@inline function decimallength(v::UInt64)
    v >= 10000000000000000 && return 17
    v >= 1000000000000000 && return 16
    v >= 100000000000000 && return 15
    v >= 10000000000000 && return 14
    v >= 1000000000000 && return 13
    v >= 100000000000 && return 12
    v >= 10000000000 && return 11
    v >= 1000000000 && return 10
    v >= 100000000 && return 9
    v >= 10000000 && return 8
    v >= 1000000 && return 7
    v >= 100000 && return 6
    v >= 10000 && return 5
    v >= 1000 && return 4
    v >= 100 && return 3
    v >= 10 && return 2
    return 1
end

@inline function decimallength(v::UInt32)
    v >= 100000000 && return 9
    v >= 10000000 && return 8
    v >= 1000000 && return 7
    v >= 100000 && return 6
    v >= 10000 && return 5
    v >= 1000 && return 4
    v >= 100 && return 3
    v >= 10 && return 2
    return 1
end

@inline function decimallength(v::UInt16)
    v >= 10000 && return 5
    v >= 1000 && return 4
    v >= 100 && return 3
    v >= 10 && return 2
    return 1
end

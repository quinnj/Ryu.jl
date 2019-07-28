@inline function umul128(a::UInt64, b::UInt64, productHi::UInt64)::UInt64
    return _umul128(a, b, productHi)
end

@inline function shiftright128(lo::UInt64, hi::UInt64, dist::UInt32)::UInt64
    return __shiftright128(lo, hi, dist % UInt)
end

@inline function div5(x::UInt64)::UInt64
    return x / 5
end

@inline function div10(x::UInt64)::UInt64
    return x / 10
end

@inline function div100(x::UInt64)::UInt64
    return x / 100
end

@inline function div1e8(x::UInt64)::UInt64
    return x / 100000000
end

@inline function div1e9(x::UInt64)::UInt64
    return x / 1000000000
end

@inline function mod1e9(x::UInt64)::UInt32
    return (x - 1000000000 * div1e9(x)) % UInt32
end

@inline function pow5Factor(value::UInt64)::UInt32
    count = UInt32(0)
    while true
        q = div5(value)
        r = (value % UInt32) - 5 * (q % UInt32)
        r != 0 && break
        value = q
        count += UInt32(1)
    end
    return count
end

@inline function multipleOfPowerOf5(value::UInt64, p::UInt32)::Bool
    return pow5Factor(value) >= p
end

@inline function multipleOfPowerOf2(value::UInt64, p::UInt32)::Bool
    return (value & ((UInt64(1) << p) - 1)) == 0
end

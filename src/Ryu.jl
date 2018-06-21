module Ryu

import Base: IEEEFloat

include("full_table.jl")

const DOUBLE_LOG10_2_DENOMINATOR = UInt64(10000000)
const DOUBLE_LOG10_2_NUMERATOR = UInt64(3010299) # DOUBLE_LOG10_2_DENOMINATOR * log_10(2)
const DOUBLE_LOG10_5_DENOMINATOR = UInt64(10000000)
const DOUBLE_LOG10_5_NUMERATOR = UInt64(6989700) # DOUBLE_LOG10_5_DENOMINATOR * log_10(5)
const DOUBLE_LOG2_5_DENOMINATOR = UInt64(10000000)
const DOUBLE_LOG2_5_NUMERATOR = UInt64(23219280) # DOUBLE_LOG2_5_DENOMINATOR * log_2(5)

const DOUBLE_POW5_INV_BITCOUNT = 122
const DOUBLE_POW5_BITCOUNT = 121

mantissabits(::Type{Float16}) = UInt32(0)
mantissabits(::Type{Float32}) = UInt32(0)
mantissabits(::Type{Float64}) = UInt32(52)

exponentbits(::Type{Float16}) = UInt32(0)
exponentbits(::Type{Float32}) = UInt32(0)
exponentbits(::Type{Float64}) = UInt32(11)

uint_t(::Type{Float16}) = UInt16
uint_t(::Type{Float32}) = UInt32
uint_t(::Type{Float64}) = UInt64

function mulShift(m::UInt64, mul::Tuple{UInt64, UInt64}, j::Int32)
    b0 = UInt128(m) * mul[1]
    b2 = UInt128(m) * mul[2]
    return (((b0 >> 64) + b2) >> (j - 64)) % UInt64
end

function mulShiftAll(m, mul, j, mmShift)
    vp = mulShift(4 * m + 2, mul, j)
    vm = mulShift(4 * m - 1 - mmShift, mul, j)
    return vp, vm, mulShift(4 * m, mul, j)
end

function pow5Factor(value)
    for count = 0:value
        if value - 5 * div(value, 5) != 0
            return count
        end
        value = div(value, 5)
    end
    return 0
end

multipleOfPowerOf5(value, p) = pow5Factor(value) >= p

function double_pow5bits(e)
    return e == 0 ? 1 :
        div(e * DOUBLE_LOG2_5_NUMERATOR + DOUBLE_LOG2_5_DENOMINATOR - 1, DOUBLE_LOG2_5_DENOMINATOR)
end

function decimalLength(v::UInt64)::UInt32
    # This is slightly faster than a loop. For a random set of numbers, the
    # average length is 17.4 digits, so we check high-to-low.
    v >= 1000000000000000000 && return 19
    v >= 100000000000000000 && return 18
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

function printfloat(f::T, debug=false) where {T <: IEEEFloat}
    result = zeros(UInt8, 25)
    mantissaBits::UInt32 = mantissabits(T)
    exponentBits::UInt32 = exponentbits(T)
    offset::UInt32 = UInt32((1 << (exponentBits - 1)) - 1)

    bits = Core.bitcast(uint_t(T), f)
    sign = ((bits >> (mantissaBits + exponentBits)) & 1) != 0
    ieeeMantissa::UInt64 = bits & ((UInt64(1) << mantissaBits) - 1)
    ieeeExponent::UInt32 = UInt32((bits >> mantissaBits) & ((1 << exponentBits) - 1))
    
    debug && println("IN=$(bitstring(bits))")
    
    e2::Int32 = 0
    m2::UInt64 = 0
    if ieeeExponent == ((1 << exponentBits) - 1)
        return ieeeMantissa != 0 ? "NaN" : sign ? "-Inf" : "Inf"
    elseif ieeeExponent == 0
        if ieeeMantissa == 0
            return sign ? "-0.0" : "0.0"
        end
        e2 = (1 - offset - mantissaBits - 2) % Int32
        m2 = UInt64(ieeeMantissa)
    else
        e2 = (ieeeExponent - offset - mantissaBits - Int32(2)) % Int32
        m2 = UInt64((1 << mantissaBits) | ieeeMantissa)
    end
    even = (m2 & 1) == 0
    acceptBounds = even

    debug && println("S=$(sign ? '-' : '+') E=$(e2 + 2) M=$(m2)")

    mv::UInt64 = UInt64(4 * m2)
    mmShift::UInt32 = UInt32((m2 != (1 << mantissaBits)) || (ieeeExponent <= 1))
    vmIsTrailingZeros = false
    vrIsTrailingZeros = false
    q::Int32 = 0
    e10::Int32 = 0
    k::Int32 = 0
    i::Int32 = 0
    j::Int32 = 0
    vr::UInt64 = 0
    vp::UInt64 = 0
    vm::UInt64 = 0
    if e2 >= 0
        q = Int32(max(0, div(e2 * DOUBLE_LOG10_2_NUMERATOR, DOUBLE_LOG10_2_DENOMINATOR) - 1))
        e10 = q
        k = Int32(DOUBLE_POW5_INV_BITCOUNT) * double_pow5bits(q) - Int32(1)
        i = Int32(-e2 + q + k)
        vp, vm, vr = mulShiftAll(m2, DOUBLE_POW5_INV_SPLIT[q+1], i, mmShift)

        if debug
            println("$mv * 2^$e2 / 10^$q")
            println("V+=$vp\nV =$vr\nV-=$vm")
        end

        if q <= 21
            if mv % 5 == 0
                vrIsTrailingZeros = multipleOfPowerOf5(mv, q)
            else
                if acceptBounds
                    vmIsTrailingZeros = multipleOfPowerOf5(mv - 1 - mmShift, q)
                else
                    vp -= multipleOfPowerOf5(mv + 2, q)
                end
            end
        end
    else
        q = Int32(max(0, div(-e2 * DOUBLE_LOG10_5_NUMERATOR, DOUBLE_LOG10_5_DENOMINATOR) - 1))
        e10 = Int32(q + e2)
        i = Int32(-e2 - q)
        k = (double_pow5bits(i) - DOUBLE_POW5_BITCOUNT) % Int32
        j = Int32(q - k)
        vp, vm, vr = mulShiftAll(m2, DOUBLE_POW5_SPLIT[i+1], j, mmShift)

        if debug
            println("$mv * 5^$(-e2) / 10^$q")
            println("$q $i $k $j")
            println("V+=$vp\nV =$vr\nV-=$vm")
        end

        if q <= 1
            vrIsTrailingZeros = (~(UInt32(mv) & 1)) >= UInt32(q)
            if acceptBounds
                vmIsTrailingZeros = (~(UInt32(mv - 1 - mmShift)) & 1) >= UInt32(q)
            else
                vp -= 1
            end
        elseif q < 63
            vrIsTrailingZeros = (mv & ((UInt64(1) << (q - 1)) - 1)) == 0
            debug && println("vr is trailing zeros=$vrIsTrailingZeros")
        end
    end

    if debug
        println("e10=$e10")
        println("V+=$vp\nV =$vr\nV-=$vm")
        println("vr is trailing zeros=$vrIsTrailingZeros")
        println("vm is trailing zeros=$vmIsTrailingZeros")
    end

    vplength::UInt32 = decimalLength(vp)
    exp::Int32 = Int32(e10 + vplength - 1)
    removed::UInt32 = UInt32(0)
    lastRemovedDigit::UInt8 = UInt8(0)
    if vmIsTrailingZeros || vrIsTrailingZeros
        # rare
        while div(vp, 10) > div(vm, 10)
            vmIsTrailingZeros &= vm - div(vm, 10) * 10 == 0
            vrIsTrailingZeros &= lastRemovedDigit == 0
            nvr = UInt64(div(vr, 10))
            lastRemovedDigit = UInt8(vr - 10 * nvr)
            vr = nvr
            vp = div(vp, 10)
            vm = div(vm, 10)
            removed += 1
        end

        debug && println("V+=$vp\nV =$vr\nV-=$vm")

        if vmIsTrailingZeros
            while vm - div(vm, 10) * 10 == 0
                vrIsTrailingZeros &= lastRemovedDigit == 0
                nvr = div(vr, 10)
                lastRemovedDigit = UInt8(vr - 10 * nvr)
                vr = nvr
                vp = div(vp, 10)
                vm = div(vm, 10)
                removed += 1
            end
        end

        debug && println("$vr $lastRemovedDigit")

        if vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0
            lastRemovedDigit = 4
        end
        output = vr + ((vr == vm && (!acceptBounds || !vmIsTrailingZeros)) || (lastRemovedDigit >= 5))
    else
        # common
        while div(vp, 10) > div(vm, 10)
            nvr = div(vr, 10)
            lastRemovedDigit = UInt8(vr - 10 * nvr)
            vr = nvr
            vp = div(vp, 10)
            vm = div(vm, 10)
            removed += 1
        end

        debug && println("$vr $lastRemovedDigit")

        output = vr + ((vr == vm) || (lastRemovedDigit >= 5))
    end
    olength::UInt32 = UInt32(vplength - removed)
    
    if debug
        println("V+=$vp\nV =$vr\nV-=$vm")
        println("O=$output")
        println("OLEN=$olength")
        println("EXP=$exp")
    end

    index = 0
    if sign
        result[index+1] = '-'
        index += 1
    end
    i = UInt32(0)
    output2::UInt64 = output
    while output2 >= 10000
        c = UInt32(output2 - 10000 * div(output2, 10000))
        output2 = div(output2, 10000)
        c0 = UInt32((c % 100) << 1)
        c1 = UInt32(div(c, 100) << 1)
        result[index + olength - i] = DIGIT_TABLE[c0+1]
        result[index + olength - i + 1] = DIGIT_TABLE[c0+2]
        result[index + olength - i - 2] = DIGIT_TABLE[c1+1]
        result[index + olength - i - 1] = DIGIT_TABLE[c1+2]
        i += 4
    end
    if output2 >= 100
        c = UInt32((output2 - 100 * div(output2, 100)) << 1)
        output2 = div(output2, 100)
        result[index + olength - i] = DIGIT_TABLE[c+1]
        result[index + olength - i + 1] = DIGIT_TABLE[c+2]
        i += 2
    end
    if output2 >= 10
        c = UInt32(output2 << 1)
        result[index + olength - i + 1] = DIGIT_TABLE[c+2]
        result[index+1] = DIGIT_TABLE[c+1]
    else
        result[index+1] = '0' + output2
    end
    if olength > 1
        result[index + 2] = '.'
        index += olength + 1
    else
        index += 1
    end
    result[index+1] = 'e'
    index += 1
    if exp < 0
        result[index+1] = '-'
        index += 1
        exp = -exp
    end
    if exp >= 100
        result[index+1] = '0' + div(exp, 100)
        index += 1
        exp = exp - 100 * div(exp, 100)
        result[index+1] = DIGIT_TABLE[2exp+1]
        result[index+2] = DIGIT_TABLE[2exp+2]
        index += 2
    elseif exp >= 10
        result[index+1] = DIGIT_TABLE[2exp+1]
        result[index+2] = DIGIT_TABLE[2exp+2]
        index += 2
    else
        result[index+1] = '0' + exp
        index += 1
    end
    return unsafe_string(pointer(result))
end

end # module

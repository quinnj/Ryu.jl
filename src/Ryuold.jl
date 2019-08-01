module RyuOld

include("utils.jl")
include("tables.jl")

mutable struct Float{T}
    mantissa::T
    exponent::Int32
end

@inline function mulShift(m::UInt64, mul, j)
    @inbounds b0 = UInt128(m) * mul[1]
    @inbounds b2 = UInt128(m) * mul[2]
    return (((b0 >> 64) + b2) >> (j - 64)) % UInt64
end

@inline function mulShift(m::UInt32, factor, shift)
    factorLo = factor % UInt32
    factorHi = (factor >> 32) % UInt32
    bits0 = UInt64(m) * factorLo
    bits1 = UInt64(m) * factorHi
    sum = (bits0 >> 32) + bits1
    shiftedSum = sum >> (shift - 32)
    return shiftedSum % UInt32
end

@inline function tofloat(::Type{T}, ieeeMantissa, ieeeExponent) where {T}
    if ieeeExponent == 0
        e2 = (1 - bias(T) - mantissabits(T) - 2) % Int32
        m2 = ieeeMantissa
    else
        e2 = (ieeeExponent - bias(T) - mantissabits(T) - 2) % Int32
        m2 = (oftype(ieeeMantissa, 1) << mantissabits(T)) | ieeeMantissa
    end
    even = (m2 & 1) == 0
    acceptBounds = even

    mv = oftype(m2, 4) * m2
    mmShift = ieeeMantissa != 0 || ieeeExponent <= 1

    vmIsTrailingZeros = false
    vrIsTrailingZeros = false
    lastRemovedDigit = 0x00
    if e2 >= 0
        q = log10pow2(e2) - (T == Float64 ? (e2 > 3) : false)
        e10 = Core.bitcast(Int32, q)
        k = pow5_inv_bitcount(T) + pow5bits(e10) - 1
        i = -e2 + e10 + k
        @show q
        mul = pow5invsplit(T, q)
        vr = mulShift(mv, mul, i)
        vp = mulShift(mv + oftype(mv, 2), mul, i)
        vm = mulShift(mv - oftype(mv, 1) - mmShift, mul, i)
        if T == Float32
            if q != 0 && div(vp - oftype(vp, 1), 10) <= div(vm, 10)
                l = pow5_inv_bitcount(T) + pow5bits(e10 - Int32(1)) - 1
                lastRemovedDigit = (mulShift(mv, pow5invsplit(T, q - 1), -e2 + Int32(q) - 1 + l) % 10) % UInt8
            end
        end
        if q <= qinvbound(T)
            if T == Float64
                mvMod5 = (mv % UInt32) - UInt32(5) * (div(mv, 5) % UInt32)
            else
                mvMod5 = mv % 5
            end
            if mvMod5 == 0
                vrIsTrailingZeros = multipleOfPowerOf5(mv, q)
            elseif acceptBounds
                vmIsTrailingZeros = multipleOfPowerOf5(mv - 1 - mmShift, q)
            else
                vp -= multipleOfPowerOf5(mv + 2, q)
            end
        end
    else
        q = log10pow5(-e2) - (T == Float64 ? (-e2 > 1) : false)
        e10 = Core.bitcast(Int32, q) + e2
        i = -e2 - Core.bitcast(Int32, q)
        k = pow5bits(i) - pow5_bitcount(T)
        j = q - k
        @show i
        mul = pow5split(T, i)
        @show mv, mv + oftype(mv, 2), mv - oftype(mv, 1) - mmShift
        vr = mulShift(mv, mul, j)
        vp = mulShift(mv + oftype(mv, 2), mul, j)
        vm = mulShift(mv - oftype(mv, 1) - mmShift, mul, j)
        if T == Float32
            if q != 0 && div(vp - oftype(vp, 1), 10) <= div(vm, 10)
                j = q - 1 - (pow5bits(i + Int32(1)) - pow5_bitcount(T))
                lastRemovedDigit = (mulShift(mv, pow5split(T, i + 1), j) % 10) % UInt8
            end
        end
        if q <= 1
            vrIsTrailingZeros = true
            if acceptBounds
                vmIsTrailingZeros = mmShift == 1
            else
                vp -= 1
            end
        elseif q < qbound(T)
            vrIsTrailingZeros = multipleOfPowerOf2(mv, q - (T == Float32))
        end
    end
    @show vp, vm, vr

    removed = 0
    if vmIsTrailingZeros || vrIsTrailingZeros
        if T == Float32
            while div(vp, 10) > div(vm, 10)
                vmIsTrailingZeros &= vm % 10 == 0
                vrIsTrailingZeros &= lastRemovedDigit == 0
                lastRemovedDigit = (vr % 10) % UInt8
                vr = div(vr, 10)
                vp = div(vp, 10)
                vm = div(vm, 10)
                removed += 1
            end
        else
            while true
                vpDiv10 = div(vp, 10)
                vmDiv10 = div(vm, 10)
                vpDiv10 <= vmDiv10 && break
                vmMod10 = (vm % UInt32) - UInt32(10) * (vmDiv10 % UInt32)
                vrDiv10 = div(vr, 10)
                vrMod10 = (vr % UInt32) - UInt32(10) * (vrDiv10 % UInt32)
                vmIsTrailingZeros &= vmMod10 == 0
                vrIsTrailingZeros &= lastRemovedDigit == 0
                lastRemovedDigit = vrMod10 % UInt8
                vr = vrDiv10
                vp = vpDiv10
                vm = vmDiv10
                removed += 1
            end
        end
        if vmIsTrailingZeros
            if T == Float32
                while vm % 10 == 0
                    vrIsTrailingZeros &= lastRemovedDigit == 0
                    lastRemovedDigit = (vr % 10) % UInt8
                    vr = div(vr, 10)
                    vp = div(vp, 10)
                    vm = div(vm, 10)
                    removed += 1
                end
            else
                while true
                    vmDiv10 = div(vm, 10)
                    vmMod10 = (vm % UInt32) - UInt32(10) * (vmDiv10 % UInt32)
                    vmMod10 != 0 && break
                    vpDiv10 = div(vp, 10)
                    vrDiv10 = div(vr, 10)
                    vrMod10 = (vr % UInt32) - UInt32(10) * (vrDiv10 % UInt32)
                    vrIsTrailingZeros &= lastRemovedDigit == 0
                    lastRemovedDigit = vrMod10 % UInt8
                    vr = vrDiv10
                    vp = vpDiv10
                    vm = vmDiv10
                    removed += 1
                end
            end
        end
        if vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0
            lastRemovedDigit = UInt8(4)
        end
        output = vr + ((vr == vm && (!acceptBounds || !vmIsTrailingZeros)) || lastRemovedDigit >= 5)
    else
        if T == Float32
            while div(vp, 10) > div(vm, 10)
                lastRemovedDigit = (vr % 10) % UInt8
                vr = div(vr, 10)
                vp = div(vp, 10)
                vm = div(vm, 10)
                removed += 1
            end
            output = vr + (vr == vm || lastRemovedDigit >= 5)
        else
            roundUp = false
            vpDiv100 = div(vp, 100)
            vmDiv100 = div(vm, 100)
            if vpDiv100 > vmDiv100
                vrDiv100 = div(vr, 100)
                vrMod100 = (vr % UInt32) - UInt32(100) * (vrDiv100 % UInt32)
                roundUp = vrMod100 >= 50
                vr = vrDiv100
                vp = vpDiv100
                vm = vmDiv100
                removed += 2
            end
            while true
                vpDiv10 = div(vp, 10)
                vmDiv10 = div(vm, 10)
                vpDiv10 <= vmDiv10 && break
                vrDiv10 = div(vr, 10)
                vrMod10 = (vr % UInt32) - UInt32(10) * (vrDiv10 % UInt32)
                roundUp = vrMod10 >= 5
                vr = vrDiv10
                vp = vpDiv10
                vm = vmDiv10
                removed += 1
            end
            output = vr + (vr == vm || roundUp)
        end
    end
    exp = (e10 + removed) % Int32
    return Float(output, exp)
end

@inline function tochars(buf, pos, v::Float{T}, sign) where {T}
    if sign > 0
        @inbounds buf[pos] = UInt8('-')
        pos += 1
    end

    output = v.mantissa
    olength = decimallength(output)

    i = 0
    if T == UInt64
        if (output >> 32) != 0
            q = output ÷ 100000000
            output2 = (output % UInt32) - UInt32(100000000) * (q % UInt32)
            output = q

            c = output2 % 10000
            output2 = div(output2, 10000)
            d = output2 % 10000
            c0 = (c % 100) << 1
            c1 = (c ÷ 100) << 1
            d0 = (d % 100) << 1
            d1 = (d ÷ 100) << 1
            unsafe_copyto!(buf, pos + olength - 1, DIGIT_TABLE, c0 + 1, 2)
            unsafe_copyto!(buf, pos + olength - 3, DIGIT_TABLE, c1 + 1, 2)
            unsafe_copyto!(buf, pos + olength - 5, DIGIT_TABLE, d0 + 1, 2)
            unsafe_copyto!(buf, pos + olength - 7, DIGIT_TABLE, d1 + 1, 2)
            i += 8
        end
        output2 = output % UInt32
    else
        output2 = output % UInt32
    end
    while output2 >= 10000
        c = output2 % 10000
        output2 = div(output2, 10000)
        c0 = (c % 100) << 1
        c1 = (c ÷ 100) << 1
        unsafe_copyto!(buf, pos + olength - i - 1, DIGIT_TABLE, c0 + 1, 2)
        unsafe_copyto!(buf, pos + olength - i - 3, DIGIT_TABLE, c1 + 1, 2)
        i += 4
    end
    if output2 >= 100
        c = (output2 % 100) << 1
        output2 = div(output2, 100)
        unsafe_copyto!(buf, pos + olength - i - 1, DIGIT_TABLE, c + 1, 2)
        i += 2
    end
    if output2 >= 10
        c = output2 << 1
        #=@inbounds=# buf[pos + olength - i] = DIGIT_TABLE[c + 2]
        #=@inbounds=# buf[pos] = DIGIT_TABLE[c + 1]
    else
        #=@inbounds=# buf[pos] = UInt8('0') + (output2 % UInt8)
    end

    if olength > 1
        #=@inbounds=# buf[pos + 1] = UInt8('.')
        pos += olength + 1
    else
        pos += 1
    end

    #=@inbounds=# buf[pos] = UInt8('e')
    pos += 1
    exp = v.exponent + olength - 1
    if exp < 0
        #=@inbounds=# buf[pos] = UInt8('-')
        pos += 1
        exp = -exp
    end

    if exp >= 100
        c = exp % 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * div(exp, 10) + 1, 2)
        #=@inbounds=# buf[pos + 2] = UInt8('0') + (c % UInt8)
        pos += 3
    elseif exp >= 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * exp + 1, 2)
        pos += 2
    else
        #=@inbounds=# buf[pos] = UInt8('0') + (exp % UInt8)
        pos += 1
    end

    return buf, pos
end

function write(x::T) where {T <: Base.IEEEFloat}
    buf, pos = write(x, zeros(UInt8, 25), 1)
    return String(buf[1:pos-1])
end

function write(x::T, buf::Vector{UInt8}, pos) where {T <: Base.IEEEFloat}
    bits = uint(x)
    ieeeSign = signbit(x)
    ieeeMantissa = bits & ((oftype(bits, 1) << mantissabits(T)) - oftype(bits, 1))
    ieeeExponent = (bits >> mantissabits(T)) & ((oftype(bits, 1) << exponentbits(T)) - oftype(bits, 1))

    # special cases
    if ieeeExponent == ((oftype(bits, 1) << exponentbits(T)) - oftype(bits, 1)) || (ieeeExponent == 0 && ieeeMantissa == 0)
        if ieeeMantissa > 0
            unsafe_copyto!(buf, pos, UInt8['N', 'a', 'N'], 1, 3)
            return buf, pos + 3
        end
        if ieeeSign > 0
            @inbounds buf[pos] = UInt8('-')
            pos += 1
        end
        if ieeeExponent > 0
            unsafe_copyto!(buf, pos, UInt8['I', 'n', 'f'], 1, 3)
            return buf, pos + 3
        end
        unsafe_copyto!(buf, pos, UInt8['0', '.', '0'], 1, 3)
        return buf, pos + 3
    end

    m2 = (oftype(bits, 1) << mantissabits(T)) | ieeeMantissa
    e2 = ieeeExponent - bias(T) - mantissabits(T)
    mask = (oftype(bits, 1) << -e2) - 1
    fraction = m2 & mask
    if e2 > 0 || e2 < -52 || fraction != 0
        v = tofloat(T, ieeeMantissa, ieeeExponent)
    else
        # small int case
        v = Float(m2 >> -e2, Int32(0))
        while true
            q = div(v.mantissa, 10)
            r = (v.mantissa % UInt32) - UInt32(10) * (q % UInt32)
            r != 0 && break
            v.mantissa = q
            v.exponent += Int32(1)
        end
    end

    return tochars(buf, pos, v, ieeeSign)
end

end # module
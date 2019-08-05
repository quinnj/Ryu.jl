module Ryu

include("utils.jl")

const MANTISSA_MASK = 0x000fffffffffffff
const EXP_MASK = 0x00000000000007ff

function write(x::T) where {T <: Base.IEEEFloat}
    buf, pos = write(x, zeros(UInt8, 25), 1)
    return String(buf[1:pos-1])
end

function writefixed(x::T, precision) where {T <: Base.IEEEFloat}
    buf, pos = writefixed(x, precision, zeros(UInt8, 2000), 1)
    return String(buf[1:pos-1])
end

function writeexp(x::T, precision) where {T <: Base.IEEEFloat}
    buf, pos = writeexp(x, precision, zeros(UInt8, 2000), 1)
    return String(buf[1:pos-1])
end

function writex(x, buf, pos)
    write(x, buf, pos)
    return
end

function writex(x, prec, buf, pos)
    writefixed(x, prec, buf, pos)
    return
end

@inline function write(x::T, buf::Vector{UInt8}, pos) where {T}
    neg = signbit(x)
    # special cases
    @inbounds if x == 0
        if neg
            buf[pos] = UInt8('-')
        end
        buf[pos + neg] = UInt8('0')
        buf[pos + neg + 1] = UInt8('e')
        buf[pos + neg + 2] = UInt8('0')
        return buf, pos + neg + 3
    elseif isnan(x)
        buf[pos] = UInt8('N')
        buf[pos + 1] = UInt8('a')
        buf[pos + 2] = UInt8('N')
        return buf, pos + 3
    elseif !isfinite(x)
        if neg
            buf[pos] = UInt8('-')
        end
        buf[pos + neg] = UInt8('I')
        buf[pos + neg + 1] = UInt8('n')
        buf[pos + neg + 2] = UInt8('f')
        return buf, pos + neg + 3
    end

    bits = uint(x)
    mant = bits & (oftype(bits, 1) << mantissabits(T) - oftype(bits, 1))
    exp = Int((bits >> mantissabits(T)) & ((1 << exponentbits(T)) - 1))

    m2 = oftype(bits, 1 << mantissabits(T)) | mant
    e2 = exp - bias(T) - mantissabits(T)
    fraction = m2 & ((oftype(bits, 1) << -e2) - 1)
    if e2 > 0 || e2 < -52 || fraction != 0
        if exp == 0
            e2 = 1 - bias(T) - mantissabits(T) - 2
            m2 = mant
        else
            e2 -= 2
        end
        even = (m2 & 1) == 0
        mv = oftype(m2, 4 * m2)
        mp = oftype(m2, mv + 2)
        mmShift = mant != 0 || exp <= 1
        mm = oftype(m2, mv - 1 - mmShift)
        vmIsTrailingZeros = false
        vrIsTrailingZeros = false
        lastRemovedDigit = 0x00
        if e2 >= 0
            q = log10pow2(e2) - (T == Float64 ? (e2 > 3) : 0)
            e10 = q
            k = pow5_inv_bitcount(T) + pow5bits(q) - 1
            i = -e2 + q + k
            vr, vp, vm = mulshiftinvsplit(T, mv, mp, mm, q, i)
            if T == Float32 || T == Float16
                if q != 0 && div(vp - 1, 10) <= div(vm, 10)
                    l = pow5_inv_bitcount(T) + pow5bits(q - 1) - 1
                    @inbounds mul = T == Float32 ? FLOAT_POW5_INV_SPLIT[q] : HALF_POW5_INV_SPLIT[q]
                    lastRemovedDigit = (mulshift(mv, mul, -e2 + q - 1 + l) % 10) % UInt8
                end
            end
            if q <= qinvbound(T)
                if ((mv % UInt32) - 5 * div(mv, 5)) == 0
                    vrIsTrailingZeros = pow5(mv, q)
                elseif even
                    vmIsTrailingZeros = pow5(mm, q)
                else
                    vp -= pow5(mp, q)
                end
            end
        else
            q = log10pow5(-e2) - (T == Float64 ? (-e2 > 1) : 0)
            e10 = q + e2
            i = -e2 - q
            k = pow5bits(i) - pow5_bitcount(T)
            j = q - k
            vr, vp, vm = mulshiftsplit(T, mv, mp, mm, i, j)
            if T == Float32 || T == Float16
                if q != 0 && div(vp - 1, 10) <= div(vm, 10)
                    j = q - 1 - (pow5bits(i + 1) - pow5_bitcount(T))
                    @inbounds mul = T == Float32 ? FLOAT_POW5_SPLIT[i + 2] : HALF_POW5_SPLIT[i + 2]
                    lastRemovedDigit = (mulshift(mv, mul, j) % 10) % UInt8
                end
            end
            if q <= 1
                vrIsTrailingZeros = true
                if even
                    vmIsTrailingZeros = mmShift
                else
                    vp -= 1
                end
            elseif q < qbound(T)
                vrIsTrailingZeros = pow2(mv, q - (T != Float64))
            end
        end
        removed = 0
        if vmIsTrailingZeros || vrIsTrailingZeros
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
            if vmIsTrailingZeros
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
            if vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0
                lastRemovedDigit = UInt8(4)
            end
            output = vr + ((vr == vm && (!even || !vmIsTrailingZeros)) || lastRemovedDigit >= 5)
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
            output = vr + (vr == vm || roundUp || lastRemovedDigit >= 5)
        end
        nexp = e10 + removed
    else
        output = m2 >> -e2
        nexp = 0
        while true
            q = div(output, 10)
            r = (output % UInt32) - UInt32(10) * (q % UInt32)
            r != 0 && break
            output = q
            nexp += 1
        end
    end

    if neg
        @inbounds buf[pos] = UInt8('-')
        pos += 1
    end
    olength = decimallength(output)

    i = 0
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
    exp2 = nexp + olength - 1
    if exp2 < 0
        #=@inbounds=# buf[pos] = UInt8('-')
        pos += 1
        exp2 = -exp2
    end

    if exp2 >= 100
        c = exp2 % 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * div(exp2, 10) + 1, 2)
        #=@inbounds=# buf[pos + 2] = UInt8('0') + (c % UInt8)
        pos += 3
    elseif exp2 >= 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * exp2 + 1, 2)
        pos += 2
    else
        #=@inbounds=# buf[pos] = UInt8('0') + (exp2 % UInt8)
        pos += 1
    end

    return buf, pos
end

@inline function writefixed(v::T, precision, buf, pos) where {T <: Base.IEEEFloat}
    x = Float64(v)
    neg = signbit(x)
    # special cases
    @inbounds if x == 0
        if neg
            buf[pos] = UInt8('-')
            pos += 1
        end
        buf[pos] = UInt8('0')
        pos += 1
        if precision > 0
            buf[pos] = UInt8('.')
            pos += 1
            for _ = 1:precision
                buf[pos] = UInt8('0')
                pos += 1
            end
        end
        return buf, pos
    elseif isnan(x)
        buf[pos] = UInt8('N')
        buf[pos + 1] = UInt8('a')
        buf[pos + 2] = UInt8('N')
        return buf, pos + 3
    elseif !isfinite(x)
        if neg
            buf[pos] = UInt8('-')
        end
        buf[pos + neg] = UInt8('I')
        buf[pos + neg + 1] = UInt8('n')
        buf[pos + neg + 2] = UInt8('f')
        return buf, pos + neg + 3
    end

    bits = Core.bitcast(UInt64, x)
    mant = bits & MANTISSA_MASK
    exp = Int((bits >> 52) & EXP_MASK)

    if exp == 0
        e2 = 1 - 1023 - 52
        m2 = mant
    else
        e2 = exp - 1023 - 52
        m2 = (1 << 52) | mant
    end
    nonzero = false
    if neg
        buf[pos] = UInt8('-')
        pos += 1
    end
    if e2 >= -52
        idx = e2 < 0 ? 0 : indexforexp(e2)
        p10bits = pow10bitsforindex(idx)
        len = lengthforindex(idx)
        i = len - 1
        while i >= 0
            j = p10bits - e2
            #=@inbounds=# mula, mulb, mulc = POW10_SPLIT[POW10_OFFSET[idx + 1] + i + 1]
            digits = mulshiftmod1e9(m2 << 8, mula, mulb, mulc, j + 8)
            if nonzero
                pos = append_nine_digits(digits, buf, pos)
            elseif digits != 0
                olength = decimallength(digits)
                pos = append_n_digits(olength, digits, buf, pos)
                nonzero = true
            end
            i -= 1
        end
    end
    if !nonzero
        buf[pos] = UInt8('0')
        pos += 1
    end
    if precision > 0
        buf[pos] = UInt8('.')
        pos += 1
    end
    if e2 < 0
        idx = div(-e2, 16)
        blocks = div(precision, 9) + 1
        roundUp = 0
        i = 0
        if blocks <= MIN_BLOCK_2[idx + 1]
            i = blocks
            for _ = 1:precision
                buf[pos] = UInt8('0')
                pos += 1
            end
        elseif i < MIN_BLOCK_2[idx + 1]
            i = MIN_BLOCK_2[idx + 1]
            for _ = 1:(9 * i)
                buf[pos] = UInt8('0')
                pos += 1
            end
        end
        while i < blocks
            j = 120 + (-e2 - 16 * idx)
            p = POW10_OFFSET_2[idx + 1] + UInt32(i) - MIN_BLOCK_2[idx + 1]
            if p >= POW10_OFFSET_2[idx + 2]
                for _ = 1:(precision - 9 * i)
                    buf[pos] = UInt8('0')
                    pos += 1
                end
                break
            end
            #=@inbounds=# mula, mulb, mulc = POW10_SPLIT_2[p + 1]
            digits = mulshiftmod1e9(m2 << 8, mula, mulb, mulc, j + 8)
            if i < blocks - 1
                pos = append_nine_digits(digits, buf, pos)
            else
                maximum = precision - 9 * i
                lastDigit = 0
                k = 0
                while k < 9 - maximum
                    # global digits, lastDigit, k
                    lastDigit = digits % 10
                    digits = div(digits, 10)
                    k += 1
                end
                if lastDigit != 5
                    roundUp = lastDigit > 5
                else
                    requiredTwos = -e2 - precision - 1
                    trailingZeros = requiredTwos <= 0 || (requiredTwos < 60 && pow2(m2, requiredTwos))
                    roundUp = trailingZeros ? 2 : 1
                end
                if maximum > 0
                    pos = append_c_digits(maximum, digits, buf, pos)
                end
                break
            end
            i += 1
        end
        if roundUp != 0
            roundPos = pos
            dotPos = 1
            while true
                roundPos -= 1
                if roundPos == 0 || (buf[roundPos] == UInt8('-'))
                    buf[roundPos + 1] = UInt8('1')
                    if dotPos > 1
                        buf[dotPos] = UInt8('0')
                        buf[dotPos + 1] = UInt8('.')
                    end
                    buf[pos] = UInt8('0')
                    pos += 1
                    break
                end
                c = roundPos > 0 ? buf[roundPos] : 0x00
                if c == UInt8('.')
                    dotPos = roundPos
                    continue
                elseif c == UInt8('9')
                    buf[roundPos] = UInt8('0')
                    roundUp = 1
                    continue
                else
                    if roundUp == 2 && UInt8(c) % 2 == 0
                        break
                    end
                    buf[roundPos] = c + 1
                    break
                end
            end
        end
    else
        for _ = 1:precision
            buf[pos] = UInt8('0')
            pos += 1
        end
    end
    return buf, pos
end

@inline function writeexp(v::T, precision, buf, pos) where {T <: Base.IEEEFloat}
    x = Float64(v)
    neg = signbit(x)
    # special cases
    @inbounds if x == 0
        if neg
            buf[pos] = UInt8('-')
            pos += 1
        end
        buf[pos] = UInt8('0')
        pos += 1
        if precision > 0
            buf[pos] = UInt8('.')
            pos += 1
            for _ = 1:precision
                buf[pos] = UInt8('0')
                pos += 1
            end
        end
        buf[pos] = UInt8('e')
        buf[pos + 1] = UInt8('+')
        buf[pos + 2] = UInt8('0')
        buf[pos + 3] = UInt8('0')
        return buf, pos + 4
    elseif isnan(x)
        buf[pos] = UInt8('N')
        buf[pos + 1] = UInt8('a')
        buf[pos + 2] = UInt8('N')
        return buf, pos + 3
    elseif !isfinite(x)
        if neg
            buf[pos] = UInt8('-')
        end
        buf[pos + neg] = UInt8('I')
        buf[pos + neg + 1] = UInt8('n')
        buf[pos + neg + 2] = UInt8('f')
        return buf, pos + neg + 3
    end

    bits = Core.bitcast(UInt64, x)
    mant = bits & MANTISSA_MASK
    exp = Int((bits >> 52) & EXP_MASK)

    if exp == 0
        e2 = 1 - 1023 - 52
        m2 = mant
    else
        e2 = exp - 1023 - 52
        m2 = (1 << 52) | mant
    end
    nonzero = false
    precision += 1
    if neg
        buf[pos] = UInt8('-')
        pos += 1
    end
    digits = 0
    printedDigits = 0
    availableDigits = 0
    e = 0
    if e2 >= -52
        idx = e2 < 0 ? 0 : indexforexp(e2)
        p10bits = pow10bitsforindex(idx)
        len = lengthforindex(idx)
        i = len - 1
        while i >= 0
            j = p10bits - e2
            #=@inbounds=# mula, mulb, mulc = POW10_SPLIT[POW10_OFFSET[idx + 1] + i + 1]
            digits = mulshiftmod1e9(m2 << 8, mula, mulb, mulc, j + 8)
            if printedDigits != 0
                if printedDigits + 9 > precision
                    availableDigits = 9
                    break
                end
                pos = append_nine_digits(digits, buf, pos)
                printedDigits += 9
            elseif digits != 0
                availableDigits = decimallength(digits)
                e = i * 9 + availableDigits - 1
                if availableDigits > precision
                    break
                end
                if precision > 1
                    pos = append_d_digits(availableDigits, digits, buf, pos)
                else
                    buf[pos] = UInt8('0') + digits
                    pos += 1
                end
                printedDigits = availableDigits
                availableDigits = 0
            end
            i -= 1
        end
    end
    if e2 < 0 && availableDigits == 0
        idx = div(-e2, 16)
        i = MIN_BLOCK_2[idx + 1]
        while i < 200
            j = 120 + (-e2 - 16 * idx)
            p = POW10_OFFSET_2[idx + 1] + i - MIN_BLOCK_2[idx + 1]
            # @show i, p, pos
            #=@inbounds=# mula, mulb, mulc = POW10_SPLIT_2[p + 1]
            digits = p >= POW10_OFFSET_2[idx + 2] ? 0 : mulshiftmod1e9(m2 << 8, mula, mulb, mulc, j + 8)
            if printedDigits != 0
                if printedDigits + 9 > precision
                    availableDigits = 9
                    break
                end
                pos = append_nine_digits(digits, buf, pos)
                printedDigits += 9
            elseif digits != 0
                availableDigits = decimallength(digits)
                e = -(i + 1) * 9 + availableDigits - 1
                if availableDigits > precision
                    break
                end
                if precision > 1
                    pos = append_d_digits(availableDigits, digits, buf, pos)
                else
                    buf[pos] = UInt8('0') + digits
                    pos += 1
                end
                printedDigits = availableDigits
                availableDigits = 0
            end
            i += 1
        end
    end
    maximum = precision - printedDigits
    if availableDigits == 0
        digits = 0
    end
    lastDigit = 0
    if availableDigits > maximum
        for k = 0:(availableDigits - maximum - 1)
            lastDigit = digits % 10
            digits = div(digits, 10)
        end
    end
    roundUp = 0
    if lastDigit != 5
        roundUp = lastDigit > 5
    else
        rexp = precision - e
        requiredTwos = -e2 - rexp
        trailingZeros = requiredTwos <= 0 ||
            (requiredTwos < 60 && pow2(m2, requiredTwos))
        if rexp < 0
            requiredFives = -rexp
            trailingZeros = trailingZeros & pow5(m2, requiredFives)
        end
        roundUp = trailingZeros ? 2 : 1
    end
    if printedDigits != 0
        if digits == 0
            for _ = 1:maximum
                buf[pos] = UInt8('0')
                pos += 1
            end
        else
            pos = append_c_digits(maximum, digits, buf, pos)
        end
    else
        if precision > 1
            pos = append_d_digits(maximum, digits, buf, pos)
        else
            buf[pos] = UInt8('0') + digits
            pos += 1
        end
    end
    if roundUp != 0
        roundPos = pos
        while true
            roundPos -= 1
            if roundPos == 0 || buf[roundPos] == UInt8('-')
                buf[roundPos + 1] = UInt8('1')
                e += 1
                break
            end
            c = roundPos > 0 ? buf[roundPos] : 0x00
            if c == UInt8('.')
                continue
            elseif c == UInt8('9')
                buf[roundPos] = UInt8('0')
                roundUp = 1
                continue
            else
                if roundUp == 2 && UInt8(c) % 2 == 0
                    break
                end
                buf[roundPos] = c + 1
                break
            end
        end
    end
    buf[pos] = UInt8('e')
    pos += 1
    if e < 0
        buf[pos] = UInt8('-')
        pos += 1
        e = -e
    else
        buf[pos] = UInt8('+')
        pos += 1
    end
    if e >= 100
        c = e % 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * div(e, 10) + 1, 2)
        buf[pos + 2] = UInt8('0') + c
        pos += 3
    else
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * e + 1, 2)
        pos += 2
    end
    return buf, pos
end

end # module
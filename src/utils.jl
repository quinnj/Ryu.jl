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

pow5_bitcount(::Type{Float16}) = 30 # ??
pow5_bitcount(::Type{Float32}) = 61
pow5_bitcount(::Type{Float64}) = 121

pow5_inv_bitcount(::Type{Float16}) = 28 # ??
pow5_inv_bitcount(::Type{Float32}) = 59
pow5_inv_bitcount(::Type{Float64}) = 122

qinvbound(::Type{Float16}) = 4 # or 3
qinvbound(::Type{Float32}) = 9
qinvbound(::Type{Float64}) = 21

qbound(::Type{Float16}) = 15
qbound(::Type{Float32}) = 31
qbound(::Type{Float64}) = 63

log10pow2(e) = (e * 78913) >> 18
log10pow5(e) = (e * 732923) >> 20
pow5bits(e) = ((e * 1217359) >> 19) + 1
mulshift(m, mula, mulb, j) = ((((UInt128(m) * mula) >> 64) + UInt128(m) * mulb) >> (j - 64)) % UInt64
mulshift(m, mul, j) = ((((m * (mul % UInt32)) >> 32) + (m * (mul >> 32))) >> (j - 32)) % UInt32
indexforexp(e) = div(e + 15, 16)
pow10bitsforindex(idx) = 16 * idx + 120
lengthforindex(idx) = div(log10pow2(16 * idx) + 1 + 16 + 8, 9)

@inline function pow5(x, p)
    count = 0
    while true
        q = div(x, 5)
        r = x - 5 * q
        r != 0 && return count >= p
        x = q
        count += 1
    end
end

pow2(x, p) = (x & ((1 << p) - 1)) == 0

@inline function decimallength(v)
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

@inline function mulshiftinvsplit(::Type{Float64}, mv, mp, mm, i, j)
    @inbounds mula, mulb = DOUBLE_POW5_INV_SPLIT[i + 1]
    vr = mulshift(mv, mula, mulb, j)
    vp = mulshift(mp, mula, mulb, j)
    vm = mulshift(mm, mula, mulb, j)
    return vr, vp, vm
end

@inline function mulshiftinvsplit(::Type{Float32}, mv, mp, mm, i, j)
    @inbounds mul = FLOAT_POW5_INV_SPLIT[i + 1]
    vr = mulshift(mv, mul, j)
    vp = mulshift(mp, mul, j)
    vm = mulshift(mm, mul, j)
    return vr, vp, vm
end

@inline function mulshiftsplit(::Type{Float64}, mv, mp, mm, i, j)
    @inbounds mula, mulb = DOUBLE_POW5_SPLIT[i + 1]
    vr = mulshift(mv, mula, mulb, j)
    vp = mulshift(mp, mula, mulb, j)
    vm = mulshift(mm, mula, mulb, j)
    return vr, vp, vm
end

@inline function mulshiftsplit(::Type{Float32}, mv, mp, mm, i, j)
    @inbounds mul = FLOAT_POW5_SPLIT[i + 1]
    vr = mulshift(mv, mul, j)
    vp = mulshift(mp, mul, j)
    vm = mulshift(mm, mul, j)
    return vr, vp, vm
end

@inline function umul256(a, bHi, bLo)
    aLo = a % UInt64
    aHi = (a >> 64) % UInt64

    b00 = UInt128(aLo) * bLo
    b01 = UInt128(aLo) * bHi
    b10 = UInt128(aHi) * bLo
    b11 = UInt128(aHi) * bHi

    b00Lo = b00 % UInt64
    b00Hi = (b00 >> 64) % UInt64

    mid1 = b10 + b00Hi
    mid1Lo = mid1 % UInt64
    mid1Hi = (mid1 >> 64) % UInt64

    mid2 = b01 + mid1Lo
    mid2Lo = mid2 % UInt64
    mid2Hi = (mid2 >> 64) % UInt64

    pHi = b11 + mid1Hi + mid2Hi
    pLo = (UInt128(mid2Lo) << 64) | b00Lo
    return pLo, pHi
end

@inline umul256_hi(a, bHi, bLo) = umul256(a, bHi, bLo)[2]

@inline function mulshiftmod1e9(m, mula, mulb, mulc, j)
    b0 = UInt128(m) * mula
    b1 = UInt128(m) * mulb
    b2 = UInt128(m) * mulc
    mid = b1 + ((b0 >> 64) % UInt64)
    s1 = b2 + ((mid >> 64) % UInt64)
    v = s1 >> (j - 128)
    multiplied = umul256_hi(v, 0x89705F4136B4A597, 0x31680A88F8953031)
    shifted = (multiplied >> 29) % UInt32
    return (v % UInt32) - UInt32(1000000000) * shifted
end

@inline function append_n_digits(olength, digits, buf, pos)
    i = 0
    while digits >= 10000
        c = digits % 10000
        digits = div(digits, 10000)
        c0 = (c % 100) << 1
        c1 = div(c, 100) << 1
        unsafe_copyto!(buf, pos + olength - i - 2, DIGIT_TABLE, c0 + 1, 2)
        unsafe_copyto!(buf, pos + olength - i - 4, DIGIT_TABLE, c1 + 1, 2)
        i += 4
    end
    if digits >= 100
        c = (digits % 100) << 1
        digits = div(digits, 100)
        unsafe_copyto!(buf, pos + olength - i - 2, DIGIT_TABLE, c + 1, 2)
        i += 2
    end
    if digits >= 10
        c = digits << 1
        unsafe_copyto!(buf, pos + olength - i - 2, DIGIT_TABLE, c + 1, 2)
        i += 2
    else
        buf[pos] = UInt8('0') + digits
        i += 1
    end
    return pos + i
end

@inline function append_d_digits(olength, digits, buf, pos)
    i = 0
    while digits >= 10000
        c = digits % 10000
        digits = div(digits, 10000)
        c0 = (c % 100) << 1
        c1 = div(c, 100) << 1
        unsafe_copyto!(buf, pos + olength + 1 - i - 2, DIGIT_TABLE, c0 + 1, 2)
        unsafe_copyto!(buf, pos + olength + 1 - i - 4, DIGIT_TABLE, c1 + 1, 2)
        i += 4
    end
    if digits >= 100
        c = (digits % 100) << 1
        digits = div(digits, 100)
        unsafe_copyto!(buf, pos + olength + 1 - i - 2, DIGIT_TABLE, c + 1, 2)
        i += 2
    end
    if digits >= 10
        c = digits << 1
        buf[pos] = DIGIT_TABLE[c + 1]
        buf[pos + 1] = UInt8('.')
        buf[pos + 2] = DIGIT_TABLE[c + 2]
        i += 3
    else
        buf[pos] = UInt8('0') + digits
        buf[pos + 1] = UInt8('.')
        i += 2
    end
    return pos + i
end

@inline function append_c_digits(count, digits, buf, pos)
    i = 0
    while i < count - 1
        c = (digits % 100) << 1
        digits = div(digits, 100)
        unsafe_copyto!(buf, pos + count - i - 2, DIGIT_TABLE, c + 1, 2)
        i += 2
    end
    if i < count
        buf[pos + count - i - 1] = UInt8('0') + (digits % 10)
        i += 1
    end
    return pos + i
end

@inline function append_nine_digits(digits, buf, pos)
    # @show String(buf[1:pos-1])
    # @show pos
    if digits == 0
        for _ = 1:9
            buf[pos] = UInt8('0')
            pos += 1
        end
        return pos
    end
    i = 0
    while i < 5
        c = digits % 10000
        digits = div(digits, 10000)
        c0 = (c % 100) << 1
        c1 = div(c, 100) << 1
        unsafe_copyto!(buf, pos + 7 - i, DIGIT_TABLE, c0 + 1, 2)
        unsafe_copyto!(buf, pos + 5 - i, DIGIT_TABLE, c1 + 1, 2)
        i += 4
    end
    buf[pos] = UInt8('0') + digits
    i += 1
    return pos + i
end
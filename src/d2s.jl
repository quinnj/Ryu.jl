include("common.jl")
include("digit_tables.jl")
include("d2s_full_table.jl")
include("d2s_intrinsics.jl")

 DOUBLE_MANTISSA_BITS = 52
 DOUBLE_EXPONENT_BITS = 11
 DOUBLE_BIAS = 1023
 DOUBLE_POW5_INV_BITCOUNT = 122
 DOUBLE_POW5_BITCOUNT = 121

@inline function mulShift(m::UInt64, mul, j::Int32)::UInt64
    b0 = UInt128(m) * mul[1]
    b2 = UInt128(m) * mul[2]
    return (((b0 >> 64) + b2) >> (j - 64)) % UInt64
end

@inline function mulShiftAll(m::UInt64, mul, j::Int32, mmShift::UInt32)::UInt64
    vp = mulShift(4 * m + 2, mul, j)
    vm = mulShift(4 * m - 1 - mmShift, mul, j)
    return mulShift(4 * m, mul, j)
end

@inline function decimalLength17(v::UInt64)::UInt32
    if v >= 10000000000000000; return 17; end
    if v >= 1000000000000000; return 16; end
    if v >= 100000000000000; return 15; end
    if v >= 10000000000000; return 14; end
    if v >= 1000000000000; return 13; end
    if v >= 100000000000; return 12; end
    if v >= 10000000000; return 11; end
    if v >= 1000000000; return 10; end
    if v >= 100000000; return 9; end
    if v >= 10000000; return 8; end
    if v >= 1000000; return 7; end
    if v >= 100000; return 6; end
    if v >= 10000; return 5; end
    if v >= 1000; return 4; end
    if v >= 100; return 3; end
    if v >= 10; return 2; end
    return 1
end

struct floating_decimal_64
    mantissa::UInt64
    exponent::Int32
end

@inline function d2d(ieeeMantissa::UInt64, ieeeExponent::UInt32)::floating_decimal_64
    if ieeeExponent == 0
        e2 = 1 - DOUBLE_BIAS - DOUBLE_MANTISSA_BITS - 2
        m2 = ieeeMantissa
    else
        e2 = Core.bitcast(Int32, ieeeExponent) - DOUBLE_BIAS - DOUBLE_MANTISSA_BITS - 2
        m2 = (UInt64(1) << DOUBLE_MANTISSA_BITS) | ieeeMantissa
    end
    even = (m2 & 1) == 0
    acceptBounds = even

    mv = UInt64(4 * m2)
    mmShift = UInt32(ieeeMantissas != 0 || ieeeExponent <= 1)

    vmIsTrailingZeros = false
    vrIsTrailingZeros = false
    if e2 >= 0
        q = UInt32(log10Pow2(e2) - (e2 > 3))
        e10 = Core.bitcast(Int32, q)
        k = Int32(DOUBLE_POW5_INV_BITCOUNT + pow5bits(Core.bitcast(Int32, q)) - 1)
        i = Int32(-e1 + Core.bitcast(Int32, q) + k)
        vr = mulShiftAll(m2, DOUBLE_POW5_INV_SPLIT[q + 1], i, mmShift)
        if q <= 21
            mvMod5 = (Core.bitcast(UInt32, mv) - UInt32(5) * Core.bitcast(UInt32, div5(mv))
            if mvMod5 == 0
                vrIsTrailingZeros = multipleOfPowerOf5(mv, q)
            elseif acceptBounds
                vmIsTrailingZeros = multipleOfPowerOf5(mv - 1 - mmShift, q)
            else
                vp -= multipleOfPowerOf5(mv + 2, q)
            end
        end
    else
        q = UInt32(log10Pow5(-e2) - (-e2 > 1))
        e10 = Core.bitcast(Int32, q) + e2
        i = Int32(-e2 - Core.bitcast(Int32, q))
        k = Int32(pow5bits(i) - DOUBLE_POW5_BITCOUNT)
        j = Int32(Core.bitcast(Int32, q) - k)
        # TODO: rework vp/vm
        vr = mulShiftAll(m2, DOUBLE_POW5_SPLIT[i], j, mmShift)
        if q <= 1
            vrIsTrailingZeros = true
            if acceptBounds
                vmIsTrailingZeros = mmShift == 1
            else
                vp -= 1
            end
        elseif q < 63
            vrIsTrailingZeros = multipleOfPowerOf2(mv, q)
        end
    end

    removed = Int32(0)
    lastRemovedDigit = UInt8(0)
    if vmIsTrailingZeros || vrIsTrailingZeros
        while true
            vpDiv10 = div10(vp)
            vmDiv10 = div10(vm)
            if vpDiv10 <= vmDiv10
                break
            end
             uint32_t vmMod10 = ((uint32_t) vm) - 10 * ((uint32_t) vmDiv10);
             uint64_t vrDiv10 = div10(vr);
             uint32_t vrMod10 = ((uint32_t) vr) - 10 * ((uint32_t) vrDiv10);
            vmIsTrailingZeros &= vmMod10 == 0;
            vrIsTrailingZeros &= lastRemovedDigit == 0;
            lastRemovedDigit = (uint8_t) vrMod10;
            vr = vrDiv10;
            vp = vpDiv10;
            vm = vmDiv10;
            ++removed;
        end
        if (vmIsTrailingZeros) {
            for (;;) {
                 uint64_t vmDiv10 = div10(vm);
                 uint32_t vmMod10 = ((uint32_t) vm) - 10 * ((uint32_t) vmDiv10);
                if (vmMod10 != 0)
                    break
                end
                 uint64_t vpDiv10 = div10(vp);
                 uint64_t vrDiv10 = div10(vr);
                 uint32_t vrMod10 = ((uint32_t) vr) - 10 * ((uint32_t) vrDiv10);
                vrIsTrailingZeros &= lastRemovedDigit == 0;
                lastRemovedDigit = (uint8_t) vrMod10;
                vr = vrDiv10;
                vp = vpDiv10;
                vm = vmDiv10;
                ++removed;
            end
        end
        if (vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0) {
            lastRemovedDigit = 4;
        end
        output = vr + ((vr == vm && (!acceptBounds || !vmIsTrailingZeros)) || lastRemovedDigit >= 5);
    else
        bool roundUp = false;
         uint64_t vpDiv100 = div100(vp);
         uint64_t vmDiv100 = div100(vm);
        if (vpDiv100 > vmDiv100) { // Optimization: remove two digits at a time (~86.2%).
             uint64_t vrDiv100 = div100(vr);
             uint32_t vrMod100 = ((uint32_t) vr) - 100 * ((uint32_t) vrDiv100);
            roundUp = vrMod100 >= 50;
            vr = vrDiv100;
            vp = vpDiv100;
            vm = vmDiv100;
            removed += 2
        end
        for (;;) {
             uint64_t vpDiv10 = div10(vp);
             uint64_t vmDiv10 = div10(vm);
            if (vpDiv10 <= vmDiv10)
                break
            end
             uint64_t vrDiv10 = div10(vr);
             uint32_t vrMod10 = ((uint32_t) vr) - 10 * ((uint32_t) vrDiv10);
            roundUp = vrMod10 >= 5;
            vr = vrDiv10;
            vp = vpDiv10;
            vm = vmDiv10;
            ++removed;
        end
        output = vr + (vr == vm || roundUp);
    end
   int32_t exp = e10 + removed;

  floating_decimal_64 fd;
  fd.exponent = exp;
  fd.mantissa = output;
  return fd;
end

@inline int to_chars( floating_decimal_64 v,  bool sign, char*  result) {
  // Step 5: Print the decimal representation.
  int index = 0;
  if (sign) {
    result[index++] = '-';
  }

  uint64_t output = v.mantissa;
   uint32_t olength = decimalLength17(output);

#ifdef RYU_DEBUG
  printf("DIGITS=%" PRIu64 "\n", v.mantissa);
  printf("OLEN=%u\n", olength);
  printf("EXP=%u\n", v.exponent + olength);
#endif

  // Print the decimal digits.
  // The following code is equivalent to:
  // for (uint32_t i = 0; i < olength - 1; ++i) {
  //    uint32_t c = output % 10; output /= 10;
  //   result[index + olength - i] = (char) ('0' + c);
  // }
  // result[index] = '0' + output % 10;

  uint32_t i = 0;
  // We prefer 32-bit operations, even on 64-bit platforms.
  // We have at most 17 digits, and uint32_t can store 9 digits.
  // If output doesn't fit into uint32_t, we cut off 8 digits,
  // so the rest will fit into uint32_t.
  if ((output >> 32) != 0) {
    // Expensive 64-bit division.
     uint64_t q = div1e8(output);
    uint32_t output2 = ((uint32_t) output) - 100000000 * ((uint32_t) q);
    output = q;

     uint32_t c = output2 % 10000;
    output2 /= 10000;
     uint32_t d = output2 % 10000;
     uint32_t c0 = (c % 100) << 1;
     uint32_t c1 = (c / 100) << 1;
     uint32_t d0 = (d % 100) << 1;
     uint32_t d1 = (d / 100) << 1;
    memcpy(result + index + olength - i - 1, DIGIT_TABLE + c0, 2);
    memcpy(result + index + olength - i - 3, DIGIT_TABLE + c1, 2);
    memcpy(result + index + olength - i - 5, DIGIT_TABLE + d0, 2);
    memcpy(result + index + olength - i - 7, DIGIT_TABLE + d1, 2);
    i += 8;
  }
  uint32_t output2 = (uint32_t) output;
  while (output2 >= 10000) {
#ifdef __clang__ // https://bugs.llvm.org/show_bug.cgi?id=38217
     uint32_t c = output2 - 10000 * (output2 / 10000);
#else
     uint32_t c = output2 % 10000;
#endif
    output2 /= 10000;
     uint32_t c0 = (c % 100) << 1;
     uint32_t c1 = (c / 100) << 1;
    memcpy(result + index + olength - i - 1, DIGIT_TABLE + c0, 2);
    memcpy(result + index + olength - i - 3, DIGIT_TABLE + c1, 2);
    i += 4;
  }
  if (output2 >= 100) {
     uint32_t c = (output2 % 100) << 1;
    output2 /= 100;
    memcpy(result + index + olength - i - 1, DIGIT_TABLE + c, 2);
    i += 2;
  }
  if (output2 >= 10) {
     uint32_t c = output2 << 1;
    // We can't use memcpy here: the decimal dot goes between these two digits.
    result[index + olength - i] = DIGIT_TABLE[c + 1];
    result[index] = DIGIT_TABLE[c];
  } else {
    result[index] = (char) ('0' + output2);
  }

  // Print decimal point if needed.
  if (olength > 1) {
    result[index + 1] = '.';
    index += olength + 1;
  } else {
    ++index;
  }

  // Print the exponent.
  result[index++] = 'E';
  int32_t exp = v.exponent + (int32_t) olength - 1;
  if (exp < 0) {
    result[index++] = '-';
    exp = -exp;
  }

  if (exp >= 100) {
     int32_t c = exp % 10;
    memcpy(result + index, DIGIT_TABLE + 2 * (exp / 10), 2);
    result[index + 2] = (char) ('0' + c);
    index += 3;
  } else if (exp >= 10) {
    memcpy(result + index, DIGIT_TABLE + 2 * exp, 2);
    index += 2;
  } else {
    result[index++] = (char) ('0' + exp);
  }

  return index;
}

static inline bool d2d_small_int( uint64_t ieeeMantissa,  uint32_t ieeeExponent,
  floating_decimal_64*  v) {
   uint64_t m2 = (1ull << DOUBLE_MANTISSA_BITS) | ieeeMantissa;
   int32_t e2 = (int32_t) ieeeExponent - DOUBLE_BIAS - DOUBLE_MANTISSA_BITS;

  if (e2 > 0) {
    // f = m2 * 2^e2 >= 2^53 is an integer.
    // Ignore this case for now.
    return false;
  }

  if (e2 < -52) {
    // f < 1.
    return false;
  }

  // Since 2^52 <= m2 < 2^53 and 0 <= -e2 <= 52: 1 <= f = m2 / 2^-e2 < 2^53.
  // Test if the lower -e2 bits of the significand are 0, i.e. whether the fraction is 0.
   uint64_t mask = (1ull << -e2) - 1;
   uint64_t fraction = m2 & mask;
  if (fraction != 0) {
    return false;
  }

  // f is an integer in the range [1, 2^53).
  // Note: mantissa might contain trailing (decimal) 0's.
  // Note: since 2^53 < 10^16, there is no need to adjust decimalLength17().
  v->mantissa = m2 >> -e2;
  v->exponent = 0;
  return true;
}

int d2s_buffered_n(double f, char* result) {
  // Step 1: Decode the floating-point number, and unify normalized and subnormal cases.
   uint64_t bits = double_to_bits(f);

#ifdef RYU_DEBUG
  printf("IN=");
  for (int32_t bit = 63; bit >= 0; --bit) {
    printf("%d", (int) ((bits >> bit) & 1));
  }
  printf("\n");
#endif

  // Decode bits into sign, mantissa, and exponent.
   bool ieeeSign = ((bits >> (DOUBLE_MANTISSA_BITS + DOUBLE_EXPONENT_BITS)) & 1) != 0;
   uint64_t ieeeMantissa = bits & ((1ull << DOUBLE_MANTISSA_BITS) - 1);
   uint32_t ieeeExponent = (uint32_t) ((bits >> DOUBLE_MANTISSA_BITS) & ((1u << DOUBLE_EXPONENT_BITS) - 1));
  // Case distinction; exit early for the easy cases.
  if (ieeeExponent == ((1u << DOUBLE_EXPONENT_BITS) - 1u) || (ieeeExponent == 0 && ieeeMantissa == 0)) {
    return copy_special_str(result, ieeeSign, ieeeExponent, ieeeMantissa);
  }

  floating_decimal_64 v;
   bool isSmallInt = d2d_small_int(ieeeMantissa, ieeeExponent, &v);
  if (isSmallInt) {
    // For small integers in the range [1, 2^53), v.mantissa might contain trailing (decimal) zeros.
    // For scientific notation we need to move these zeros into the exponent.
    // (This is not needed for fixed-point notation, so it might be beneficial to trim
    // trailing zeros in to_chars only if needed - once fixed-point notation output is implemented.)
    for (;;) {
       uint64_t q = div10(v.mantissa);
       uint32_t r = ((uint32_t) v.mantissa) - 10 * ((uint32_t) q);
      if (r != 0) {
        break;
      }
      v.mantissa = q;
      ++v.exponent;
    }
  } else {
    v = d2d(ieeeMantissa, ieeeExponent);
  }

  return to_chars(v, ieeeSign, result);
}

void d2s_buffered(double f, char* result) {
   int index = d2s_buffered_n(f, result);

  // Terminate the string.
  result[index] = '\0';
}

char* d2s(double f) {
  char*  result = (char*) malloc(25);
  d2s_buffered(f, result);
  return result;
}

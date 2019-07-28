@inline function decimalLength9(v::UInt32)::UInt32
  # Function precondition: v is not a 10-digit number.
  # (f2s: 9 digits are sufficient for round-tripping.)
  # (d2fixed: We print 9-digit blocks.)
  if (v >= 100000000); return 9; end
  if (v >= 10000000); return 8; end
  if (v >= 1000000); return 7; end
  if (v >= 100000); return 6; end
  if (v >= 10000); return 5; end
  if (v >= 1000); return 4; end
  if (v >= 100); return 3; end
  if (v >= 10); return 2; end
  return 1
end

# Returns e == 0 ? 1 : ceil(log_2(5^e)).
@inline function pow5bits(e::Int32)::Int32
  # This approximation works up to the point that the multiplication overflows at e = 3529.
  # If the multiplication were done in 64 bits, it would fail at 5^4004 which is just greater
  # than 2^9297.
  return ((Core.bitcast(UInt32, e) * 1217359) >> 19) + 1
end

# Returns floor(log_10(2^e)).
@inline function log10Pow2(e::Int32)::UInt32
  # The first value this approximation fails for is 2^1651 which is just greater than 10^297.
  return (Core.bitcast(UInt32, e) * 78913) >> 18
end

# Returns floor(log_10(5^e)).
@inline function log10Pow5(e::Int32)::UInt32
  # The first value this approximation fails for is 5^2621 which is just greater than 10^1832.
  return (Core.bitcast(UInt32, e) * 732923) >> 20
end

memcpy(dest, src, n) = ccall(:memcpy, Cvoid, (Ptr{UInt8}, Ptr{UInt8}, Csize_t), dest, src, n)

@inline function copy_special_str(char * const result, const bool sign, const bool exponent, const bool mantissa)::Cint
    if mantissa
        copyto!(result, b"NaN")
        return 3
    end
    if sign
        @inbounds result[1] = UInt8('-')
    end
    if exponent
        copyto!(result, sign, b"Infinity", 0, 8)
        return sign + 8
    end
    copyto!(result, sign, b"0E0", 0, 3)
    return sign + 3
end

@inline float_to_bits(f::Float32) = Core.bitcast(UInt32, f)
@inline double_to_bits(d::Float64) = Core.bitcast(UInt64, d)


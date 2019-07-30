using Ryu, Test

const maxMantissa = (UInt64(1) << 53) - 1
todouble(sign, exp, mant) = Core.bitcast(Float64, (UInt64(sign) << 63) | (UInt64(exp) << 52) | (UInt64(mant)))

@testset "Ryu" begin

@testset "Float64" begin

@testset "Basic" begin
    @test Ryu.write(0.0) == "0.0"
    @test Ryu.write(-0.0) == "-0.0"
    @test Ryu.write(1.0) == "1e0"
    @test Ryu.write(-1.0) == "-1e0"
    @test Ryu.write(NaN) == "NaN"
    @test Ryu.write(Inf) == "Inf"
    @test Ryu.write(-Inf) == "-Inf"
end

@testset "SwitchToSubnormal" begin
    @test "2.2250738585072014e-308" == Ryu.write(2.2250738585072014e-308)
end

@testset "MinAndMax" begin
    @test "1.7976931348623157e308" == Ryu.write(Core.bitcast(Float64, 0x7fefffffffffffff))
    @test "5e-324" == Ryu.write(Core.bitcast(Float64, 1))
end

@testset "LotsOfTrailingZeros" begin
    @test "2.9802322387695312e-8" == Ryu.write(2.98023223876953125e-8)
end

@testset "Regression" begin
    @test "-2.109808898695963e16" == Ryu.write(-2.109808898695963e16)
    @test "4.940656e-318" == Ryu.write(4.940656e-318)
    @test "1.18575755e-316" == Ryu.write(1.18575755e-316)
    @test "2.989102097996e-312" == Ryu.write(2.989102097996e-312)
    @test "9.0608011534336e15" == Ryu.write(9.0608011534336e15)
    @test "4.708356024711512e18" == Ryu.write(4.708356024711512e18)
    @test "9.409340012568248e18" == Ryu.write(9.409340012568248e18)
    @test "1.2345678e0" == Ryu.write(1.2345678)
end

@testset "LooksLikePow5" begin
    # These numbers have a mantissa that is a multiple of the largest power of 5 that fits,
    # and an exponent that causes the computation for q to result in 22, which is a corner
    # case for Ryu.
    @test "5.764607523034235e39" == Ryu.write(Core.bitcast(Float64, 0x4830F0CF064DD592))
    @test "1.152921504606847e40" == Ryu.write(Core.bitcast(Float64, 0x4840F0CF064DD592))
    @test "2.305843009213694e40" == Ryu.write(Core.bitcast(Float64, 0x4850F0CF064DD592))
end

@testset "OutputLength" begin
    @test "1e0" == Ryu.write(1.0) # already tested in Basic
    @test "1.2e0" == Ryu.write(1.2)
    @test "1.23e0" == Ryu.write(1.23)
    @test "1.234e0" == Ryu.write(1.234)
    @test "1.2345e0" == Ryu.write(1.2345)
    @test "1.23456e0" == Ryu.write(1.23456)
    @test "1.234567e0" == Ryu.write(1.234567)
    @test "1.2345678e0" == Ryu.write(1.2345678) # already tested in Regressi
    @test "1.23456789e0" == Ryu.write(1.23456789)
    @test "1.234567895e0" == Ryu.write(1.234567895) # 1.234567890 would be trimm
    @test "1.2345678901e0" == Ryu.write(1.2345678901)
    @test "1.23456789012e0" == Ryu.write(1.23456789012)
    @test "1.234567890123e0" == Ryu.write(1.234567890123)
    @test "1.2345678901234e0" == Ryu.write(1.2345678901234)
    @test "1.23456789012345e0" == Ryu.write(1.23456789012345)
    @test "1.234567890123456e0" == Ryu.write(1.234567890123456)
    @test "1.2345678901234567e0" == Ryu.write(1.2345678901234567)

  # Test 32-bit chunking
    @test "4.294967294e0" == Ryu.write(4.294967294) # 2^32 -
    @test "4.294967295e0" == Ryu.write(4.294967295) # 2^32 -
    @test "4.294967296e0" == Ryu.write(4.294967296) # 2^
    @test "4.294967297e0" == Ryu.write(4.294967297) # 2^32 +
    @test "4.294967298e0" == Ryu.write(4.294967298) # 2^32 +
end

# Test min, max shift values in shiftright128
@testset "MinMaxShift" begin
    # 32-bit opt-size=0:  49 <= dist <= 50
    # 32-bit opt-size=1:  30 <= dist <= 50
    # 64-bit opt-size=0:  50 <= dist <= 50
    # 64-bit opt-size=1:  30 <= dist <= 50
    @test "1.7800590868057611e-307" == Ryu.write(todouble(false, 4, 0))
    # 32-bit opt-size=0:  49 <= dist <= 49
    # 32-bit opt-size=1:  28 <= dist <= 49
    # 64-bit opt-size=0:  50 <= dist <= 50
    # 64-bit opt-size=1:  28 <= dist <= 50
    @test "2.8480945388892175e-306" == Ryu.write(todouble(false, 6, maxMantissa))
    # 32-bit opt-size=0:  52 <= dist <= 53
    # 32-bit opt-size=1:   2 <= dist <= 53
    # 64-bit opt-size=0:  53 <= dist <= 53
    # 64-bit opt-size=1:   2 <= dist <= 53
    @test "2.446494580089078e-296" == Ryu.write(todouble(false, 41, 0))
    # 32-bit opt-size=0:  52 <= dist <= 52
    # 32-bit opt-size=1:   2 <= dist <= 52
    # 64-bit opt-size=0:  53 <= dist <= 53
    # 64-bit opt-size=1:   2 <= dist <= 53
    @test "4.8929891601781557e-296" == Ryu.write(todouble(false, 40, maxMantissa))

    # 32-bit opt-size=0:  57 <= dist <= 58
    # 32-bit opt-size=1:  57 <= dist <= 58
    # 64-bit opt-size=0:  58 <= dist <= 58
    # 64-bit opt-size=1:  58 <= dist <= 58
    @test "1.8014398509481984e16" == Ryu.write(todouble(false, 1077, 0))
    # 32-bit opt-size=0:  57 <= dist <= 57
    # 32-bit opt-size=1:  57 <= dist <= 57
    # 64-bit opt-size=0:  58 <= dist <= 58
    # 64-bit opt-size=1:  58 <= dist <= 58
    @test "3.6028797018963964e16" == Ryu.write(todouble(false, 1076, maxMantissa))
    # 32-bit opt-size=0:  51 <= dist <= 52
    # 32-bit opt-size=1:  51 <= dist <= 59
    # 64-bit opt-size=0:  52 <= dist <= 52
    # 64-bit opt-size=1:  52 <= dist <= 59
    @test "2.900835519859558e-216" == Ryu.write(todouble(false, 307, 0))
    # 32-bit opt-size=0:  51 <= dist <= 51
    # 32-bit opt-size=1:  51 <= dist <= 59
    # 64-bit opt-size=0:  52 <= dist <= 52
    # 64-bit opt-size=1:  52 <= dist <= 59
    @test "5.801671039719115e-216" == Ryu.write(todouble(false, 306, maxMantissa))

    # https:#github.com/ulfjack/ryu/commit/19e44d16d80236f5de25800f56d82606d1be00b9#commitcomment-30146483
    # 32-bit opt-size=0:  49 <= dist <= 49
    # 32-bit opt-size=1:  44 <= dist <= 49
    # 64-bit opt-size=0:  50 <= dist <= 50
    # 64-bit opt-size=1:  44 <= dist <= 50
    @test "3.196104012172126e-27" == Ryu.write(todouble(false, 934, 0x000FA7161A4D6E0C))
end

@testset "SmallIntegers" begin
    @test "9.007199254740991e15" == Ryu.write(9007199254740991.0)
    @test "9.007199254740992e15" == Ryu.write(9007199254740992.0)

    @test "1e0" == Ryu.write(1.0e+0)
    @test "1.2e1" == Ryu.write(1.2e+1)
    @test "1.23e2" == Ryu.write(1.23e+2)
    @test "1.234e3" == Ryu.write(1.234e+3)
    @test "1.2345e4" == Ryu.write(1.2345e+4)
    @test "1.23456e5" == Ryu.write(1.23456e+5)
    @test "1.234567e6" == Ryu.write(1.234567e+6)
    @test "1.2345678e7" == Ryu.write(1.2345678e+7)
    @test "1.23456789e8" == Ryu.write(1.23456789e+8)
    @test "1.23456789e9" == Ryu.write(1.23456789e+9)
    @test "1.234567895e9" == Ryu.write(1.234567895e+9)
    @test "1.2345678901e10" == Ryu.write(1.2345678901e+10)
    @test "1.23456789012e11" == Ryu.write(1.23456789012e+11)
    @test "1.234567890123e12" == Ryu.write(1.234567890123e+12)
    @test "1.2345678901234e13" == Ryu.write(1.2345678901234e+13)
    @test "1.23456789012345e14" == Ryu.write(1.23456789012345e+14)
    @test "1.234567890123456e15" == Ryu.write(1.234567890123456e+15)

  # 10^i
    @test "1e0" == Ryu.write(1.0e+0)
    @test "1e1" == Ryu.write(1.0e+1)
    @test "1e2" == Ryu.write(1.0e+2)
    @test "1e3" == Ryu.write(1.0e+3)
    @test "1e4" == Ryu.write(1.0e+4)
    @test "1e5" == Ryu.write(1.0e+5)
    @test "1e6" == Ryu.write(1.0e+6)
    @test "1e7" == Ryu.write(1.0e+7)
    @test "1e8" == Ryu.write(1.0e+8)
    @test "1e9" == Ryu.write(1.0e+9)
    @test "1e10" == Ryu.write(1.0e+10)
    @test "1e11" == Ryu.write(1.0e+11)
    @test "1e12" == Ryu.write(1.0e+12)
    @test "1e13" == Ryu.write(1.0e+13)
    @test "1e14" == Ryu.write(1.0e+14)
    @test "1e15" == Ryu.write(1.0e+15)

  # 10^15 + 10^i
    @test "1.000000000000001e15" == Ryu.write(1.0e+15 + 1.0e+0)
    @test "1.00000000000001e15" == Ryu.write(1.0e+15 + 1.0e+1)
    @test "1.0000000000001e15" == Ryu.write(1.0e+15 + 1.0e+2)
    @test "1.000000000001e15" == Ryu.write(1.0e+15 + 1.0e+3)
    @test "1.00000000001e15" == Ryu.write(1.0e+15 + 1.0e+4)
    @test "1.0000000001e15" == Ryu.write(1.0e+15 + 1.0e+5)
    @test "1.000000001e15" == Ryu.write(1.0e+15 + 1.0e+6)
    @test "1.00000001e15" == Ryu.write(1.0e+15 + 1.0e+7)
    @test "1.0000001e15" == Ryu.write(1.0e+15 + 1.0e+8)
    @test "1.000001e15" == Ryu.write(1.0e+15 + 1.0e+9)
    @test "1.00001e15" == Ryu.write(1.0e+15 + 1.0e+10)
    @test "1.0001e15" == Ryu.write(1.0e+15 + 1.0e+11)
    @test "1.001e15" == Ryu.write(1.0e+15 + 1.0e+12)
    @test "1.01e15" == Ryu.write(1.0e+15 + 1.0e+13)
    @test "1.1e15" == Ryu.write(1.0e+15 + 1.0e+14)

  # Largest power of 2 <= 10^(i+1)
    @test "8e0" == Ryu.write(8.0)
    @test "6.4e1" == Ryu.write(64.0)
    @test "5.12e2" == Ryu.write(512.0)
    @test "8.192e3" == Ryu.write(8192.0)
    @test "6.5536e4" == Ryu.write(65536.0)
    @test "5.24288e5" == Ryu.write(524288.0)
    @test "8.388608e6" == Ryu.write(8388608.0)
    @test "6.7108864e7" == Ryu.write(67108864.0)
    @test "5.36870912e8" == Ryu.write(536870912.0)
    @test "8.589934592e9" == Ryu.write(8589934592.0)
    @test "6.8719476736e10" == Ryu.write(68719476736.0)
    @test "5.49755813888e11" == Ryu.write(549755813888.0)
    @test "8.796093022208e12" == Ryu.write(8796093022208.0)
    @test "7.0368744177664e13" == Ryu.write(70368744177664.0)
    @test "5.62949953421312e14" == Ryu.write(562949953421312.0)
    @test "9.007199254740992e15" == Ryu.write(9007199254740992.0)

  # 1000 * (Largest power of 2 <= 10^(i+1))
    @test "8e3" == Ryu.write(8.0e+3)
    @test "6.4e4" == Ryu.write(64.0e+3)
    @test "5.12e5" == Ryu.write(512.0e+3)
    @test "8.192e6" == Ryu.write(8192.0e+3)
    @test "6.5536e7" == Ryu.write(65536.0e+3)
    @test "5.24288e8" == Ryu.write(524288.0e+3)
    @test "8.388608e9" == Ryu.write(8388608.0e+3)
    @test "6.7108864e10" == Ryu.write(67108864.0e+3)
    @test "5.36870912e11" == Ryu.write(536870912.0e+3)
    @test "8.589934592e12" == Ryu.write(8589934592.0e+3)
    @test "6.8719476736e13" == Ryu.write(68719476736.0e+3)
    @test "5.49755813888e14" == Ryu.write(549755813888.0e+3)
    @test "8.796093022208e15" == Ryu.write(8796093022208.0e+3)
end

end # Float64

@testset "Float32" begin

@testset "Basic" begin
    @test "0.0" == Ryu.write(Float32(0.0))
    @test "-0.0" == Ryu.write(Float32(-0.0))
    @test "1e0" == Ryu.write(Float32(1.0))
    @test "-1e0" == Ryu.write(Float32(-1.0))
    @test "NaN" == Ryu.write(Float32(NaN))
    @test "Inf" == Ryu.write(Float32(Inf))
    @test "-Inf" == Ryu.write(Float32(-Inf))
end

@testset "SwitchToSubnormal" begin
    @test "1.1754944e-38" == Ryu.write(Float32(1.1754944e-38))
end

@testset "MinAndMax" begin
    @test "3.4028235e38" == Ryu.write(Core.bitcast(Float32, 0x7f7fffff))
    @test "1e-45" == Ryu.write(Core.bitcast(Float32, Int32(1)))
end

# Check that we return the exact boundary if it is the shortest
# representation, but only if the original floating point number is even.
@testset "BoundaryRoundeven" begin
    @test "3.355445e7" == Ryu.write(Float32(3.355445e7))
    @test "9e9" == Ryu.write(Float32(8.999999e9))
    @test "3.436672e10" == Ryu.write(Float32(3.4366717e10))
end

# If the exact value is exactly halfway between two shortest representations,
# then we round to even. It seems like this only makes a difference if the
# last two digits are ...2|5 or ...7|5, and we cut off the 5.
@testset "exactValueRoundeven" begin
    @test "3.0540412e5" == Ryu.write(Float32(3.0540412e5))
    @test "8.0990312e3" == Ryu.write(Float32(8.0990312e3))
end

@testset "LotsOfTrailingZeros" begin
    # Pattern for the first test: 00111001100000000000000000000000
    @test "2.4414062e-4" == Ryu.write(Float32(2.4414062e-4))
    @test "2.4414062e-3" == Ryu.write(Float32(2.4414062e-3))
    @test "4.3945312e-3" == Ryu.write(Float32(4.3945312e-3))
    @test "6.3476562e-3" == Ryu.write(Float32(6.3476562e-3))
end

@testset "Regression" begin
    @test "4.7223665e21" == Ryu.write(Float32(4.7223665e21))
    @test "8.388608e6" == Ryu.write(Float32(8388608.0))
    @test "1.6777216e7" == Ryu.write(Float32(1.6777216e7))
    @test "3.3554436e7" == Ryu.write(Float32(3.3554436e7))
    @test "6.7131496e7" == Ryu.write(Float32(6.7131496e7))
    @test "1.9310392e-38" == Ryu.write(Float32(1.9310392e-38))
    @test "-2.47e-43" == Ryu.write(Float32(-2.47e-43))
    @test "1.993244e-38" == Ryu.write(Float32(1.993244e-38))
    @test "4.1039004e3" == Ryu.write(Float32(4103.9003))
    @test "5.3399997e9" == Ryu.write(Float32(5.3399997e9))
    @test "6.0898e-39" == Ryu.write(Float32(6.0898e-39))
    @test "1.0310042e-3" == Ryu.write(Float32(0.0010310042))
    @test "2.882326e17" == Ryu.write(Float32(2.8823261e17))
    @test "7.038531e-26" == Ryu.write(Float32(7.0385309e-26))
    @test "9.223404e17" == Ryu.write(Float32(9.2234038e17))
    @test "6.710887e7" == Ryu.write(Float32(6.7108872e7))
    @test "1e-44" == Ryu.write(Float32(1.0e-44))
    @test "2.816025e14" == Ryu.write(Float32(2.816025e14))
    @test "9.223372e18" == Ryu.write(Float32(9.223372e18))
    @test "1.5846086e29" == Ryu.write(Float32(1.5846085e29))
    @test "1.1811161e19" == Ryu.write(Float32(1.1811161e19))
    @test "5.368709e18" == Ryu.write(Float32(5.368709e18))
    @test "4.6143166e18" == Ryu.write(Float32(4.6143165e18))
    @test "7.812537e-3" == Ryu.write(Float32(0.007812537))
    @test "1e-45" == Ryu.write(Float32(1.4e-45))
    @test "1.18697725e20" == Ryu.write(Float32(1.18697724e20))
    @test "1.00014165e-36" == Ryu.write(Float32(1.00014165e-36))
    @test "2e2" == Ryu.write(Float32(200.0))
    @test "3.3554432e7" == Ryu.write(Float32(3.3554432e7))
end

@testset "LooksLikePow5" begin
    # These numbers have a mantissa that is the largest power of 5 that fits,
    # and an exponent that causes the computation for q to result in 10, which is a corner
    # case for Ryu.
    @test "6.7108864e17" == Ryu.write(Core.bitcast(Float32, 0x5D1502F9))
    @test "1.3421773e18" == Ryu.write(Core.bitcast(Float32, 0x5D9502F9))
    @test "2.6843546e18" == Ryu.write(Core.bitcast(Float32, 0x5E1502F9))
end

@testset "OutputLength" begin
    @test "1e0" == Ryu.write(Float32(1.0))
    @test "1.2e0" == Ryu.write(Float32(1.2))
    @test "1.23e0" == Ryu.write(Float32(1.23))
    @test "1.234e0" == Ryu.write(Float32(1.234))
    @test "1.2345e0" == Ryu.write(Float32(1.2345))
    @test "1.23456e0" == Ryu.write(Float32(1.23456))
    @test "1.234567e0" == Ryu.write(Float32(1.234567))
    @test "1.2345678e0" == Ryu.write(Float32(1.2345678))
    @test "1.23456735e-36" == Ryu.write(Float32(1.23456735e-36))
end

end # Float32

end # Ryu


buffer = Base.Grisu.getbuf()
bignums = Base.Grisu.BIGNUMS[Threads.threadid()]
Base.Grisu.grisu(1.23456735e-36, Base.Grisu.SHORTEST, 0, buffer, bignums)
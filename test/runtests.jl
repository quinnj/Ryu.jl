using Ryu, Test

const maxMantissa = (UInt64(1) << 53) - 1
todouble(sign, exp, mant) = Core.bitcast(Float64, (UInt64(sign) << 63) | (UInt64(exp) << 52) | (UInt64(mant)))

@testset "Ryu" begin

@testset "Float64" begin

@testset "Basic" begin
    @test Ryu.write(0.0) == "0e0"
    @test Ryu.write(-0.0) == "-0e0"
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

# @testset "Float32" begin

# @testset "Basic" begin
#     @test "0.0" == Ryu.write(Float32(0.0))
#     @test "-0.0" == Ryu.write(Float32(-0.0))
#     @test "1e0" == Ryu.write(Float32(1.0))
#     @test "-1e0" == Ryu.write(Float32(-1.0))
#     @test "NaN" == Ryu.write(Float32(NaN))
#     @test "Inf" == Ryu.write(Float32(Inf))
#     @test "-Inf" == Ryu.write(Float32(-Inf))
# end

# @testset "SwitchToSubnormal" begin
#     @test "1.1754944e-38" == Ryu.write(Float32(1.1754944e-38))
# end

# @testset "MinAndMax" begin
#     @test "3.4028235e38" == Ryu.write(Core.bitcast(Float32, 0x7f7fffff))
#     @test "1e-45" == Ryu.write(Core.bitcast(Float32, Int32(1)))
# end

# # Check that we return the exact boundary if it is the shortest
# # representation, but only if the original floating point number is even.
# @testset "BoundaryRoundeven" begin
#     @test "3.355445e7" == Ryu.write(Float32(3.355445e7))
#     @test "9e9" == Ryu.write(Float32(8.999999e9))
#     @test "3.436672e10" == Ryu.write(Float32(3.4366717e10))
# end

# # If the exact value is exactly halfway between two shortest representations,
# # then we round to even. It seems like this only makes a difference if the
# # last two digits are ...2|5 or ...7|5, and we cut off the 5.
# @testset "exactValueRoundeven" begin
#     @test "3.0540412e5" == Ryu.write(Float32(3.0540412e5))
#     @test "8.0990312e3" == Ryu.write(Float32(8.0990312e3))
# end

# @testset "LotsOfTrailingZeros" begin
#     # Pattern for the first test: 00111001100000000000000000000000
#     @test "2.4414062e-4" == Ryu.write(Float32(2.4414062e-4))
#     @test "2.4414062e-3" == Ryu.write(Float32(2.4414062e-3))
#     @test "4.3945312e-3" == Ryu.write(Float32(4.3945312e-3))
#     @test "6.3476562e-3" == Ryu.write(Float32(6.3476562e-3))
# end

# @testset "Regression" begin
#     @test "4.7223665e21" == Ryu.write(Float32(4.7223665e21))
#     @test "8.388608e6" == Ryu.write(Float32(8388608.0))
#     @test "1.6777216e7" == Ryu.write(Float32(1.6777216e7))
#     @test "3.3554436e7" == Ryu.write(Float32(3.3554436e7))
#     @test "6.7131496e7" == Ryu.write(Float32(6.7131496e7))
#     @test "1.9310392e-38" == Ryu.write(Float32(1.9310392e-38))
#     @test "-2.47e-43" == Ryu.write(Float32(-2.47e-43))
#     @test "1.993244e-38" == Ryu.write(Float32(1.993244e-38))
#     @test "4.1039004e3" == Ryu.write(Float32(4103.9003))
#     @test "5.3399997e9" == Ryu.write(Float32(5.3399997e9))
#     @test "6.0898e-39" == Ryu.write(Float32(6.0898e-39))
#     @test "1.0310042e-3" == Ryu.write(Float32(0.0010310042))
#     @test "2.882326e17" == Ryu.write(Float32(2.8823261e17))
#     @test "7.038531e-26" == Ryu.write(Float32(7.0385309e-26))
#     @test "9.223404e17" == Ryu.write(Float32(9.2234038e17))
#     @test "6.710887e7" == Ryu.write(Float32(6.7108872e7))
#     @test "1e-44" == Ryu.write(Float32(1.0e-44))
#     @test "2.816025e14" == Ryu.write(Float32(2.816025e14))
#     @test "9.223372e18" == Ryu.write(Float32(9.223372e18))
#     @test "1.5846086e29" == Ryu.write(Float32(1.5846085e29))
#     @test "1.1811161e19" == Ryu.write(Float32(1.1811161e19))
#     @test "5.368709e18" == Ryu.write(Float32(5.368709e18))
#     @test "4.6143166e18" == Ryu.write(Float32(4.6143165e18))
#     @test "7.812537e-3" == Ryu.write(Float32(0.007812537))
#     @test "1e-45" == Ryu.write(Float32(1.4e-45))
#     @test "1.18697725e20" == Ryu.write(Float32(1.18697724e20))
#     @test "1.00014165e-36" == Ryu.write(Float32(1.00014165e-36))
#     @test "2e2" == Ryu.write(Float32(200.0))
#     @test "3.3554432e7" == Ryu.write(Float32(3.3554432e7))
# end

# @testset "LooksLikePow5" begin
#     # These numbers have a mantissa that is the largest power of 5 that fits,
#     # and an exponent that causes the computation for q to result in 10, which is a corner
#     # case for Ryu.
#     @test "6.7108864e17" == Ryu.write(Core.bitcast(Float32, 0x5D1502F9))
#     @test "1.3421773e18" == Ryu.write(Core.bitcast(Float32, 0x5D9502F9))
#     @test "2.6843546e18" == Ryu.write(Core.bitcast(Float32, 0x5E1502F9))
# end

# @testset "OutputLength" begin
#     @test "1e0" == Ryu.write(Float32(1.0))
#     @test "1.2e0" == Ryu.write(Float32(1.2))
#     @test "1.23e0" == Ryu.write(Float32(1.23))
#     @test "1.234e0" == Ryu.write(Float32(1.234))
#     @test "1.2345e0" == Ryu.write(Float32(1.2345))
#     @test "1.23456e0" == Ryu.write(Float32(1.23456))
#     @test "1.234567e0" == Ryu.write(Float32(1.234567))
#     @test "1.2345678e0" == Ryu.write(Float32(1.2345678))
#     @test "1.23456735e-36" == Ryu.write(Float32(1.23456735e-36))
# end

# end # Float32

@testset "Ryu.writefixed" begin
    @testset "Basic" begin
        @test Ryu.writefixed(todouble(false, 1234, 99999), 0) ==
            "3291009114715486435425664845573426149758869524108446525879746560"
    end
    @testset "Zero" begin
        @test Ryu.writefixed(0.0, 4) == "0.0000"
        @test Ryu.writefixed(0.0, 3) == "0.000"
        @test Ryu.writefixed(0.0, 2) == "0.00"
        @test Ryu.writefixed(0.0, 1) == "0.0"
        @test Ryu.writefixed(0.0, 0) == "0"
    end
    @testset "MinMax" begin
        @test Ryu.writefixed(todouble(false, 0, 1), 1074) ==
            "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" *
            "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" *
            "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" *
            "000000000000000000000000000000000000000000000000000000049406564584124654417656879286822137" *
            "236505980261432476442558568250067550727020875186529983636163599237979656469544571773092665" *
            "671035593979639877479601078187812630071319031140452784581716784898210368871863605699873072" *
            "305000638740915356498438731247339727316961514003171538539807412623856559117102665855668676" *
            "818703956031062493194527159149245532930545654440112748012970999954193198940908041656332452" *
            "475714786901472678015935523861155013480352649347201937902681071074917033322268447533357208" *
            "324319360923828934583680601060115061698097530783422773183292479049825247307763759272478746" *
            "560847782037344696995336470179726777175851256605511991315048911014510378627381672509558373" *
            "89733598993664809941164205702637090279242767544565229087538682506419718265533447265625"
        @test Ryu.writefixed(todouble(false, 2046, 0xFFFFFFFFFFFFF), 0) ==
            "179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558" *
            "632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245" *
            "490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168" *
            "738177180919299881250404026184124858368"
    end
    @testset "RoundToEven" begin
        @test Ryu.writefixed(0.125, 3) == "0.125"
        @test Ryu.writefixed(0.125, 2) == "0.12"
        @test Ryu.writefixed(0.375, 3) == "0.375"
        @test Ryu.writefixed(0.375, 2) == "0.38"
    end
    @testset "RoundToEvenInteger" begin
        @test Ryu.writefixed(2.5, 1) == "2.5"
        @test Ryu.writefixed(2.5, 0) == "2"
        @test Ryu.writefixed(3.5, 1) == "3.5"
        @test Ryu.writefixed(3.5, 0) == "4"
    end
    @testset "NonRoundToEvenScenarios" begin
        @test Ryu.writefixed(0.748046875, 3) == "0.748"
        @test Ryu.writefixed(0.748046875, 2) == "0.75"
        @test Ryu.writefixed(0.748046875, 1) == "0.7"

        @test Ryu.writefixed(0.2509765625, 3) == "0.251"
        @test Ryu.writefixed(0.2509765625, 2) == "0.25"
        @test Ryu.writefixed(0.2509765625, 1) == "0.3"

        @test Ryu.writefixed(todouble(false, 1021, 1), 54) == "0.250000000000000055511151231257827021181583404541015625"
        @test Ryu.writefixed(todouble(false, 1021, 1),  3) == "0.250"
        @test Ryu.writefixed(todouble(false, 1021, 1),  2) == "0.25"
        @test Ryu.writefixed(todouble(false, 1021, 1),  1) == "0.3"
    end
    @testset "VaryingPrecision" begin
        @test Ryu.writefixed(1729.142857142857, 47) == "1729.14285714285711037518922239542007446289062500000"
        @test Ryu.writefixed(1729.142857142857, 46) == "1729.1428571428571103751892223954200744628906250000"
        @test Ryu.writefixed(1729.142857142857, 45) == "1729.142857142857110375189222395420074462890625000"
        @test Ryu.writefixed(1729.142857142857, 44) == "1729.14285714285711037518922239542007446289062500"
        @test Ryu.writefixed(1729.142857142857, 43) == "1729.1428571428571103751892223954200744628906250"
        @test Ryu.writefixed(1729.142857142857, 42) == "1729.142857142857110375189222395420074462890625"
        @test Ryu.writefixed(1729.142857142857, 41) == "1729.14285714285711037518922239542007446289062"
        @test Ryu.writefixed(1729.142857142857, 40) == "1729.1428571428571103751892223954200744628906"
        @test Ryu.writefixed(1729.142857142857, 39) == "1729.142857142857110375189222395420074462891"
        @test Ryu.writefixed(1729.142857142857, 38) == "1729.14285714285711037518922239542007446289"
        @test Ryu.writefixed(1729.142857142857, 37) == "1729.1428571428571103751892223954200744629"
        @test Ryu.writefixed(1729.142857142857, 36) == "1729.142857142857110375189222395420074463"
        @test Ryu.writefixed(1729.142857142857, 35) == "1729.14285714285711037518922239542007446"
        @test Ryu.writefixed(1729.142857142857, 34) == "1729.1428571428571103751892223954200745"
        @test Ryu.writefixed(1729.142857142857, 33) == "1729.142857142857110375189222395420074"
        @test Ryu.writefixed(1729.142857142857, 32) == "1729.14285714285711037518922239542007"
        @test Ryu.writefixed(1729.142857142857, 31) == "1729.1428571428571103751892223954201"
        @test Ryu.writefixed(1729.142857142857, 30) == "1729.142857142857110375189222395420"
        @test Ryu.writefixed(1729.142857142857, 29) == "1729.14285714285711037518922239542"
        @test Ryu.writefixed(1729.142857142857, 28) == "1729.1428571428571103751892223954"
        @test Ryu.writefixed(1729.142857142857, 27) == "1729.142857142857110375189222395"
        @test Ryu.writefixed(1729.142857142857, 26) == "1729.14285714285711037518922240"
        @test Ryu.writefixed(1729.142857142857, 25) == "1729.1428571428571103751892224"
        @test Ryu.writefixed(1729.142857142857, 24) == "1729.142857142857110375189222"
        @test Ryu.writefixed(1729.142857142857, 23) == "1729.14285714285711037518922"
        @test Ryu.writefixed(1729.142857142857, 22) == "1729.1428571428571103751892"
        @test Ryu.writefixed(1729.142857142857, 21) == "1729.142857142857110375189"
        @test Ryu.writefixed(1729.142857142857, 20) == "1729.14285714285711037519"
        @test Ryu.writefixed(1729.142857142857, 19) == "1729.1428571428571103752"
        @test Ryu.writefixed(1729.142857142857, 18) == "1729.142857142857110375"
        @test Ryu.writefixed(1729.142857142857, 17) == "1729.14285714285711038"
        @test Ryu.writefixed(1729.142857142857, 16) == "1729.1428571428571104"
        @test Ryu.writefixed(1729.142857142857, 15) == "1729.142857142857110"
        @test Ryu.writefixed(1729.142857142857, 14) == "1729.14285714285711"
        @test Ryu.writefixed(1729.142857142857, 13) == "1729.1428571428571"
        @test Ryu.writefixed(1729.142857142857, 12) == "1729.142857142857"
        @test Ryu.writefixed(1729.142857142857, 11) == "1729.14285714286"
        @test Ryu.writefixed(1729.142857142857, 10) == "1729.1428571429"
        @test Ryu.writefixed(1729.142857142857,  9) == "1729.142857143"
        @test Ryu.writefixed(1729.142857142857,  8) == "1729.14285714"
        @test Ryu.writefixed(1729.142857142857,  7) == "1729.1428571"
        @test Ryu.writefixed(1729.142857142857,  6) == "1729.142857"
        @test Ryu.writefixed(1729.142857142857,  5) == "1729.14286"
        @test Ryu.writefixed(1729.142857142857,  4) == "1729.1429"
        @test Ryu.writefixed(1729.142857142857,  3) == "1729.143"
        @test Ryu.writefixed(1729.142857142857,  2) == "1729.14"
        @test Ryu.writefixed(1729.142857142857,  1) == "1729.1"
        @test Ryu.writefixed(1729.142857142857,  0) == "1729"
    end

    @testset "Carrying" begin
        @test Ryu.writefixed(  0.0009, 4) == "0.0009"
        @test Ryu.writefixed(  0.0009, 3) == "0.001"
        @test Ryu.writefixed(  0.0029, 4) == "0.0029"
        @test Ryu.writefixed(  0.0029, 3) == "0.003"
        @test Ryu.writefixed(  0.0099, 4) == "0.0099"
        @test Ryu.writefixed(  0.0099, 3) == "0.010"
        @test Ryu.writefixed(  0.0299, 4) == "0.0299"
        @test Ryu.writefixed(  0.0299, 3) == "0.030"
        @test Ryu.writefixed(  0.0999, 4) == "0.0999"
        @test Ryu.writefixed(  0.0999, 3) == "0.100"
        @test Ryu.writefixed(  0.2999, 4) == "0.2999"
        @test Ryu.writefixed(  0.2999, 3) == "0.300"
        @test Ryu.writefixed(  0.9999, 4) == "0.9999"
        @test Ryu.writefixed(  0.9999, 3) == "1.000"
        @test Ryu.writefixed(  2.9999, 4) == "2.9999"
        @test Ryu.writefixed(  2.9999, 3) == "3.000"
        @test Ryu.writefixed(  9.9999, 4) == "9.9999"
        @test Ryu.writefixed(  9.9999, 3) == "10.000"
        @test Ryu.writefixed( 29.9999, 4) == "29.9999"
        @test Ryu.writefixed( 29.9999, 3) == "30.000"
        @test Ryu.writefixed( 99.9999, 4) == "99.9999"
        @test Ryu.writefixed( 99.9999, 3) == "100.000"
        @test Ryu.writefixed(299.9999, 4) == "299.9999"
        @test Ryu.writefixed(299.9999, 3) == "300.000"

        @test Ryu.writefixed(  0.09, 2) == "0.09"
        @test Ryu.writefixed(  0.09, 1) == "0.1"
        @test Ryu.writefixed(  0.29, 2) == "0.29"
        @test Ryu.writefixed(  0.29, 1) == "0.3"
        @test Ryu.writefixed(  0.99, 2) == "0.99"
        @test Ryu.writefixed(  0.99, 1) == "1.0"
        @test Ryu.writefixed(  2.99, 2) == "2.99"
        @test Ryu.writefixed(  2.99, 1) == "3.0"
        @test Ryu.writefixed(  9.99, 2) == "9.99"
        @test Ryu.writefixed(  9.99, 1) == "10.0"
        @test Ryu.writefixed( 29.99, 2) == "29.99"
        @test Ryu.writefixed( 29.99, 1) == "30.0"
        @test Ryu.writefixed( 99.99, 2) == "99.99"
        @test Ryu.writefixed( 99.99, 1) == "100.0"
        @test Ryu.writefixed(299.99, 2) == "299.99"
        @test Ryu.writefixed(299.99, 1) == "300.0"

        @test Ryu.writefixed(  0.9, 1) == "0.9"
        @test Ryu.writefixed(  0.9, 0) == "1"
        @test Ryu.writefixed(  2.9, 1) == "2.9"
        @test Ryu.writefixed(  2.9, 0) == "3"
        @test Ryu.writefixed(  9.9, 1) == "9.9"
        @test Ryu.writefixed(  9.9, 0) == "10"
        @test Ryu.writefixed( 29.9, 1) == "29.9"
        @test Ryu.writefixed( 29.9, 0) == "30"
        @test Ryu.writefixed( 99.9, 1) == "99.9"
        @test Ryu.writefixed( 99.9, 0) == "100"
        @test Ryu.writefixed(299.9, 1) == "299.9"
        @test Ryu.writefixed(299.9, 0) == "300"
    end

    @testset "RoundingResultZero" begin
        @test Ryu.writefixed(0.004, 3) == "0.004"
        @test Ryu.writefixed(0.004, 2) == "0.00"
        @test Ryu.writefixed(0.4, 1) == "0.4"
        @test Ryu.writefixed(0.4, 0) == "0"
        @test Ryu.writefixed(0.5, 1) == "0.5"
        @test Ryu.writefixed(0.5, 0) == "0"
    end

    @testset "Regression" begin
        @test Ryu.writefixed(7.018232e-82, 6) == "0.000000"
    end

end

end # Ryu


# buffer = Base.Grisu.getbuf()
# bignums = Base.Grisu.BIGNUMS[Threads.threadid()]
# Base.Grisu.grisu(1.23456735e-36, Base.Grisu.SHORTEST, 0, buffer, bignums)

# x = 1.23456735e-36
# Ryu.write(x, buffer, 1)
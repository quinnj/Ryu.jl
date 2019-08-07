module Ryu

include("utils.jl")
include("shortest.jl")
include("fixed.jl")
include("exp.jl")

shortestdigits(::Type{Float64}) = 309 + 17
shortestdigits(::Type{Float32}) = 39 + 9
shortestdigits(::Type{Float16}) = 9 + 5

function writeshortest(x::T) where {T <: Base.IEEEFloat}
    buf = Vector{UInt8}(undef, shortestdigits(T))
    pos = writeshortest(x, buf, 1)
    return unsafe_string(pointer(buf), pos-1)
end

function writefixed(x::T, precision) where {T <: Base.IEEEFloat}
    buf = Vector{UInt8}(undef, precision + shortestdigits(T))
    pos = writefixed(x, precision, buf, 1)
    return unsafe_string(pointer(buf), pos-1)
end

function writeexp(x::T, precision) where {T <: Base.IEEEFloat}
    buf = Vector{UInt8}(undef, precision + shortestdigits(T))
    pos = writeexp(x, precision, buf, 1)
    return unsafe_string(pointer(buf), pos-1)
end

function Base.show(io::IO, x::T) where {T <: Base.IEEEFloat}
    if get(io, :compact, false)
        precision = T == Float16 ? 5 : 6
        buf = Vector{UInt8}(undef, precision + shortestdigits(T))
        pos = writefixed(x, precision, buf, 1)
        GC.@preserve buf unsafe_write(io, pointer(buf), pos - 1)
    else
        buf = Vector{UInt8}(undef, shortestdigits(T))
        pos = writeshortest(x, buf, 1)
        GC.@preserve buf unsafe_write(io, pointer(buf), pos - 1)
    end
    return
end

end # module
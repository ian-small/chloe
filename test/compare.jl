using Test
import Test: @test

include("../src/UInt8Utf8.jl")
include("../src/MMappedString.jl")

import .MappedString: ASCII, MMappedString

function gen_ascii(n)
    rand(UInt8(0):UInt8(127), n)
end
function rrange(m::Int)
    l, h = rand(1:m, 2)
    if l > h
        h, l = l, h
    end
    return l, h
end

@inline function compareSubStrings(
    a::SubString{MMappedString{ASCII}}, 
    b::SubString{MMappedString{ASCII}})::Tuple{Int32,Int}
    la, lb = length(a), length(b)
    m = min(la, lb)

    count::Int32 = 0
    one::Int32 = 1
    ia = a.string.ptr
    ib = b.string.ptr
    # compare bytes!
    @inbounds for i in 1:m
        c = ia[i + a.offset]
        d = ib[i + b.offset]
        c ≠ d && return c < d ? (count, -1) : (count, 1)
        count += one
    end
    la == 0 && return lb == 0 ? (count, 0) : (count, -1)
    return (count, 1)
end

@inline function compareSubStrings(a::SubString, b::SubString)::Tuple{Int32,Int}
    count::Int32 = 0
    one::Int32 = 1
    @inbounds for (c, d) in zip(a, b)
        c ≠ d && return c < d ? (count, -1) : (count, 1)
        count += one
    end
    isempty(a) && return isempty(b) ? (count, 0) : (count, -1)
    return (count, 1)
end

function timeit()
    n = 200000
    t = 4000
    v = gen_ascii(n)
    m = MMappedString{ASCII}(v)
    s = String(copy(v)) # copy to stop setting v to 0

    @test s == m
    println("seq length: $(length(m)) with $(t)x$(t) comparisons")
    ranges = [rrange(n) for i in 1:t]
    ss = map(r -> SubString(s, r...), ranges)
    mm = map(r -> SubString(m, r...), ranges)

    for i in 1:length(ranges)
        for j in 1:length(ranges)
            @test compareSubStrings(mm[i], mm[j]) == compareSubStrings(ss[i], ss[j])
        end
    end
    @time begin
        for i in 1:length(ranges)
            for j in 1:length(ranges)
                compareSubStrings(mm[i], mm[j])
            end
        end
    end
    @time begin
        for i in 1:length(ranges)
            for j in 1:length(ranges)
                compareSubStrings(ss[i], ss[j])
            end
        end
    end

end
timeit()
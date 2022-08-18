module Spectral

export gaussian, distrib_switch, cossignal

const TWOPI = 2 * π

raw"""
    gaussian(x[, μ, σ])

Compute the gaussian function for the give x, μ, σ parameters.
    ``y = \frac{1}{σ\sqrt{2π}}e^{-\frac{1}{2}\left(\frac{x - μ}{σ}\right)^2}``

# Arguments
- `x::Real`: The x value to compute the gaussian with.
- `μ::Real=0`: The mean/center of the gaussian function.
- `σ::Real=1`: The standard deviation of the gaussian function.

# Examples
```julia-repl
julia> gaussian(0.)
0.3989422804014327
```
"""
function gaussian(x::T, μ::T, σ::T)::T where T <: Real
    iσ = one(T) / σ
    norm = iσ / sqrt(TWOPI)::T
    ratio = (x - μ) * iσ
    norm * exp(T(-0.5) * ratio ^ 2)
end
gaussian(x) = gaussian(x, typeof(x)(0.), typeof(x)(1.))

raw"""
    cossignal(x, ampl, freq[, normalize])

Compute periodique signal from sum of cosine function
    ``y = \sim_i a_icos(2πν_ix)``

# Arguments
    - `x::Array{Real}`: The x value to compute the signal with.
    - `a::Array{Real}`: Array contaning the amplitude of the cosine functions.
    - `ν::Array{Real}`: Array contaning the frequencies of the cosine functions.
    - `normalize::Bool=false`: If true, normalize the signal function in order to ensure -1 < y < 1 at each point.

# Trow
    - `DimensionMismatch`: If a and ν length does not match
"""
function cossignal(x::Array{T}, a::Array{T}, ν::Array{T}, normalize::Bool)::Array{T} where T <: Real
    if length(a) != length(ν)
        throw(DimensionMismatch("a and ν must have same length, got: $(length(a)) $(length(ν))"))
    end
    if normalize
        a = a / sum(a)
    end
    y = zeros(size(x))
    for (ampl, freq) in zip(a, ν)
        y += ampl * cos.(TWOPI * freq * x)
    end
    y
end
cossignal(x, a, ν) = cossignal(x, a, ν, false)

function distrib_switch(name::String)::Function
    if name == "gaussian"
        return gaussian
    elseif name == "signal"
        return function (x, a, b) cossignal(x, [a b], abs([a b])) end
    end
    return gaussian
end

greet(name::String) = "Hello $name !"
hello_world() = greet("World")

end

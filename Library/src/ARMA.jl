"""ARMA process"""
function ARMA(X⃗::AbstractVector, ξ⃗::AbstractVector, ϕ⃗::AbstractVector, θ⃗::AbstractVector) # No constant, too boring
    # Just one update of an ARMA process
    if !all(length(X⃗)-1 .== (length(ξ⃗)-1, length(ϕ⃗), length(θ⃗)))
        p = max(length(X⃗)-1, length(ξ⃗)-1, length(ϕ⃗), length(θ⃗))
        X⃗ = vcat(X⃗, zeros(p+1-length(X⃗)))
        ξ⃗ = vcat(ξ⃗, randn(p+1-length(ξ⃗)))
        ϕ⃗ = vcat(ϕ⃗, zeros(p-length(ϕ⃗)))
        θ⃗ = vcat(θ⃗, zeros(p-length(θ⃗)))
    end
    X⃗[2:end] = X⃗[1:end-1]
    ξ⃗[2:end] = ξ⃗[1:end-1]
    ξ⃗[1] = randn()
    X⃗[1] = ξ⃗[1] + sum(ϕ⃗.*X⃗[2:end]) + sum(θ⃗.*ξ⃗[2:end])
    return (X⃗, ξ⃗)
end

function ARMA(X⃗::AbstractVector, ξ⃗::AbstractVector, ϕ₁::Float64=0.0, θ₁::Float64=0.0, args...)
    # Args are higher order parameters, alternating between ϕᵢ and θᵢ
    # E.g. ARMA(X⃗, ξ⃗, ϕ₁, θ₁, ϕ₂, θ₂)
    ϕ⃗ = vcat(ϕ₁, args[1:2:end])
    θ⃗ = vcat(θ₁, args[2:2:end])
    ARMA(X⃗, ξ⃗, ϕ⃗, θ⃗)
end

function ARMA(X0::AbstractVector, ξ0::AbstractVector, ϕ::Vector{Float64}, θ::Vector{Float64}, T::Int)
    X = zeros(T, length(X0))
    X[1, :] = X0
    for t ∈ 2:T
        X[t, :], ξ0 = ARMA(X[t-1, :], ξ0, ϕ, θ)
    end
    return X
end

function ARMA(X0::AbstractVector, ϕ::Vector{Float64}, θ::Vector{Float64}, T::Int)
    X = zeros(T, length(X0))
    X[1, :] = X0
    ξ = randn(length(X0))
    ARMA(X0, ξ, ϕ, θ, T)
end

"""AR process"""
function AR(X⃗::AbstractVector, ξ⃗::AbstractVector, ϕ⃗::AbstractVector)
    θ⃗ = zeros(length(ϕ⃗))
    ARMA(X⃗, ξ⃗, ϕ⃗, θ⃗)
end

function AR(X⃗::AbstractVector, ξ⃗::AbstractVector, ϕ₁::Float64=0.0, args...)
    ϕ⃗ = vcat(ϕ₁, args)
    AR(X⃗, ξ⃗, ϕ⃗)
end

function AR(X0::AbstractVector, ξ0::AbstractVector, ϕ::Vector{Float64}, T::Int)
    θ = zeros(length(ϕ))
    ARMA(X0, ξ0, ϕ, θ, T)
end

function AR(X0::AbstractVector, ϕ::Vector{Float64}, T::Int)
    θ = zeros(length(ϕ))
    ARMA(X0, ξ, ϕ, θ, T)
end

"""MA process"""
function MA(X⃗::AbstractVector, ξ⃗::AbstractVector, θ⃗::AbstractVector)
    ϕ⃗ = zeros(length(θ⃗))
    ARMA(X⃗, ξ⃗, ϕ⃗, θ⃗)
end

function MA(X⃗::AbstractVector, ξ⃗::AbstractVector, θ₁::Float64=0.0, args...)
    θ⃗ = vcat(θ₁, args)
    MA(X⃗, ξ⃗, θ⃗)
end

function MA(X0::AbstractVector, ξ0::AbstractVector, ϕ::Vector{Float64}, T::Int)
    ϕ = zeros(length(θ))
    ARMA(X0, ξ0, ϕ, θ, T)
end

function MA(X0::AbstractVector, ϕ::Vector{Float64}, T::Int)
    ϕ = zeros(length(θ))s
    ARMA(X0, ξ, ϕ, θ, T)
end

arma(args...) = ARMA(args...)
ar(args...) = AR(args...)
ma(args...) = MA(args...)

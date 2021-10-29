module CartesianSphericalHarmonics

# Write your package code here.
using LinearAlgebra
using Reexport

@reexport using MultivariatePolynomials
@reexport using TypedPolynomials

export YLMNorm, Schmidt, Laplace, Nonorm, Full
export rlm, ylm
#normalization factor
abstract type YLMNorm{T} end

struct Schmidt{T} <: YLMNorm{T}; end
struct Laplace{T} <: YLMNorm{T}; end
struct Nonorm{T} <: YLMNorm{T}; end
struct Full{T} <: YLMNorm{T}; end

function ylmcoeff(N::Schmidt{T},l::Integer, m::Integer) where T

	k = one(T)
	for i in (l-m+1):(l+m)
	  k *= i
	end
	fac = sqrt(1/k)
	if m!=0
		fac*=sqrt(2*one(T))
	end
	return fac
end

function ylmcoeff(N::Laplace{T}, l::Integer, m::Integer) where T

  k = one(T)
  for i in (l-m+1):(l+m)
    k *= i
  end

  return sqrt((2*l+1) / (4*T(pi)*k))
end

function ylmcoeff(N::Full{T}, l::Integer, m::Integer) where T
	fac = m==0 ? one(T) : sqrt(2*one(T))
	return ylmcoeff(Laplace{T}(),l,m)*sqrt(4T(pi))*fac
end

function ylmcoeff(N::Nonorm{T}, l::Integer, m::Integer) where T
 return one(T)
end



function cossinpoly(m::Integer, x, y)

  sum = zero(x*y)
  for j in 0:div(m,2)
	  k = 2j
    sum += ((-1)^j)*binomial(m, k)*y^Int(k)*x^Int(m-k)
  end
  return sum
end

function sinsinpoly(m::Integer, x, y)

  sum = zero(x*y)
  for j in 0:div((m-1),2)
	  k = 2j+1
    sum += ((-1)^j)*binomial(m, k)*y^Int(k)*x^Int(m-k)
  end
  return sum
end

#
function legendre_assoc(::Type{T},z,l,m) where T
	p = (z^2 - 1)^l

	for i = 1:l+abs(m)
	c = i <= l ? 1/(2one(T)*i) : one(T)
	p = c*differentiate(p, z)
	end

	p *= (-1)^m
	return p
end
# #legendre polynomial:
# P(l::BigInt,x) = differentiate((x^2 - 1)^l,x,l)
#
# #associated legendre polynomial (without sin(θ)^m)
# P(l::BigInt,m::BigInt,x) = 1//(big(2)^l*factorial(l))*differentiate(P(l,x),x,m)

"""
    ylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the spherical harmonic for a given order (l,m) in Cartesian coordinates\\

*Input:*  `l`       - Degree of the spherical harmonic\\
          `m`       - Order of the spherical harmonic\\
          `x, y, z` - Cartesian coordinates\\

*Output:*  Spherical harmonic polynomial
"""
function ylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable;
			 norm::YLMNorm{T}=Nonorm{Rational{Int}}(),real::Bool=true) where T

	if abs(m) > l
		throw(DomainError(m,"-l <= m <= l expected, but m = $m and l = $l."))
	end

	p = legendre_assoc(T,z,l,m)
	# p = (z^2 - 1)^l

	# for i = 1:l+abs(m)
	#   c = i <= l ? 1/(2one(T)*i) : one(T)
	#   p = c*differentiate(p, z)
	# end

	# p *= (-1)^m

	if real
		if m > 0
			return ylmcoeff(norm, l, m)*cossinpoly(m,x,y)*p
		elseif m < 0
			return ylmcoeff(norm, l, abs(m))*sinsinpoly(abs(m),x,y)*p
		else
			return ylmcoeff(norm, l, 0)*p
		end
	else
		return ylmcoeff(norm, l, m)*(cossinpoly(m,x,y)+im*sinsinpoly(m,x,y))*p
	end
end

# multiplying r^l*ylm(x,y,z)
function rlm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable;
			   norm::YLMNorm{T}=Nonorm{Rational{Int}}(),real::Bool=true) where T
	p = ylm(l,m,x,y,z;norm=norm,real=real)
	tout = []

	for t in terms(p)
		deg = degree(monomial(t)) # degree of monomial
		degR = l-deg # degree of r: l-deg
		push!(tout,(x^2+y^2+z^2)^div(degR,2)*t) # r² replaced by x²+y²+z²
	end

	return polynomial(tout)
end

end

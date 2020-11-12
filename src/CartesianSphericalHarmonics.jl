module CartesianSphericalHarmonics

# Write your package code here.
using LinearAlgebra, Reexport
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
	return sqrt(1/k)
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



function cossinpoly(m::Integer, x::Variable, y::Variable)

  sum = zero(x*y)
  for j in 0:div(m,2)
	  k = 2j
    sum += ((-1)^j)*binomial(m, k)*y^Int(k)*x^Int(m-k)
  end
  return sum
end

function sinsinpoly(m::Integer, x::Variable, y::Variable)

  sum = zero(x*y)
  for j in 0:div((m-1),2)
	  k = 2j+1
    sum += ((-1)^j)*binomial(m, k)*y^Int(k)*x^Int(m-k)
  end
  return sum
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

	p = (z^2 - 1)^l

	for i = 1:l+abs(m)
	  c = i <= l ? 1/(2one(T)*i) : one(T)
	  p = c*differentiate(p, z)
	end

	p *= (-1)^m

	if real
		if m > 0
			out = ylmcoeff(norm, l, m)*cossinpoly(m,x,y)*p
			if norm != Nonorm{T}()
				out *= sqrt(2one(T))
			end
			return out
		elseif m < 0
			out = ylmcoeff(norm, l, abs(m))*sinsinpoly(abs(m),x,y)*p
			if norm != Nonorm{T}()
				out *= sqrt(2one(T))
			end
			return out
		else
			return ylmcoeff(norm, l, 0)*p
		end
	else
		return ylmcoeff(norm, l, m)*(cossinpoly(m,x,y)+im*sinsinpoly(m,x,y))*p
	end
end

# multiplying r^l*ylm(x,y,z)
function rlylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable;
			   norm::YLMNorm{T}=Nonorm{Rational{Int}}(),real::Bool=true) where T
	p = ylm(l,m,x,y,z;norm=norm,real=real)
	tout = []
	# Zerlegung des Polynoms in Terme:
	for t in terms(p)
		deg = degree(monomial(t)) # Gibt den gesamten Grad des Monoms an
		degR = l-deg # durch das Kürzen ergibt sich ein Grad von l-deg fuer r
		push!(tout,(x^2+y^2+z^2)^div(degR,2)*t) # r² wird durch x²+y²+z² ersetzt
	end

	return polynomial(tout)
end

# solid harmonics
function rlm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable;
			 norm::YLMNorm{T}=Nonorm{Rational{Int}}(),real::Bool=true) where T
	rlm = rlylm(l,m,x,y,z; norm=norm,real=real)
	if norm != Nonorm{T}()
		rlm = sqrt(4one(T)*T(pi)/(2*l+1))*rlm
	end
	return rlm
end

end

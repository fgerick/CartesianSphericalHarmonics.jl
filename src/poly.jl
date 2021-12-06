
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
  
  function legendre_assoc(::Type{T},z,l,m) where T
      p = (z^2 - 1)^l
  
      for i = 1:l+abs(m)
      c = i <= l ? 1/(2one(T)*i) : one(T)
      p = c*differentiate(p, z)
      end
  
      p *= (-1)^m
      return p
  end
  
  """
      ylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable)
  *Description:*  Calculation of the spherical harmonic for a given order (l,m) in Cartesian coordinates\\
  
  *Input:*  `l`       - Degree of the spherical harmonic\\
            `m`       - Order of the spherical harmonic\\
            `x, y, z` - Unit Cartesian coordinates on the unit sphere\\
  
  *Output:*  Spherical harmonic polynomial
  """
  function ylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable;
               norm::Type{YLMNorm{T}} = Nonorm{Rational{Int}}, real::Bool=true) where T
  
      if abs(m) > l
          throw(DomainError(m,"-l <= m <= l expected, but m = $m and l = $l."))
      end
  
      p = legendre_assoc(T,z,l,m)
  
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
  
  """
      rlm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable; norm, real)
  *Description:*  Calculation of the solid spherical harmonic for a given order (l,m) in Cartesian coordinates\\
  
  *Input:*  `l`       - Degree of the spherical harmonic\\
            `m`       - Order of the spherical harmonic\\
            `x, y, z` - Cartesian coordinates\\

  
  *Output:*  Solid spherical harmonic polynomial
  """
  function rlm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable;
                 norm::Type{YLMNorm{T}} = Nonorm{Rational{Int}}, real::Bool=true) where T
      p = ylm(l,m,x,y,z; norm=norm, real=real)
      tout = []
  
      for t in terms(p)
          deg = degree(monomial(t)) # degree of monomial
          degR = l-deg # degree of r: l-deg
          push!(tout,(x^2+y^2+z^2)^div(degR,2)*t) # r² replaced by x²+y²+z²
      end
  
      return polynomial(tout)
  end
  
  
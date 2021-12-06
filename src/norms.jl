abstract type YLMNorm{T} end

struct Schmidt{T} <: YLMNorm{T}; end
struct Laplace{T} <: YLMNorm{T}; end
struct Nonorm{T} <: YLMNorm{T}; end
struct Full{T} <: YLMNorm{T}; end

function ylmcoeff(::Type{Schmidt{T}},l::Integer, m::Integer) where T

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

function ylmcoeff(::Type{Laplace{T}}, l::Integer, m::Integer) where T

  k = one(T)
  for i in (l-m+1):(l+m)
    k *= i
  end

  return sqrt((2*l+1) / (4*T(pi)*k))
end

function ylmcoeff(::Type{Full{T}}, l::Integer, m::Integer) where T
	fac = m==0 ? one(T) : sqrt(2*one(T))
	return ylmcoeff(Laplace{T},l,m)*sqrt(4T(pi))*fac
end

function ylmcoeff(::Type{Nonorm{T}}, l::Integer, m::Integer) where T
 return one(T)
end



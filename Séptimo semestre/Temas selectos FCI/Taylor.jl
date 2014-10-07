type Taylor{T<:Number}
    coeffs::Array{T,1}
    order :: Int
    #constructor interno
    function Taylor(coeffs::Array{T,1},order::Int)
        ll=length(coeffs)
        order=max(ll-1,order)
        v=zeros(T,order+1)
        v[1:ll]=coeffs[1:ll]
        new(v,order)
    end
end

Taylor{T<:Number}(v::Array{T,1}, n::Int)=Taylor{T}(v,n)
Taylor{T<:Number}(v::Array{T,1})=Taylor{T}(v,0)
Taylor{T<:Number}(c::T, n::Int)=Taylor{T}([c],n)
Taylor{T<:Number}(c::T)=Taylor{T}([c],0)
Taylor{T<:Number}(t::Taylor{T},n::Int)=Taylor{T}(t.coeffs, n)
Taylor{T<:Number}(t::Taylor{T})=Taylor{T}(t.coeffs, 0)


function +(f::Taylor,g::Taylor)
    m=max(f.order,g.order)
    if m==f.order
        g=Taylor(g,m)
    else
        f=Taylor(f,m)
    end
    Taylor(f.coeffs+g.coeffs,m)
end
+(f::Taylor)=f
+(c::Number,f::Taylor)=Taylor(c)+f
+(f::Taylor,c::Number)=f+Taylor(c)

function -(f::Taylor,g::Taylor)
    fcoeffs,g=promote
    m=max(f.order,g.order)
    if m==f.order
        g=Taylor(g,m)
    else
        f=Taylor(f,m)
    end
    Taylor(f.coeffs-g.coeffs,m)
end
-(f::Taylor)=Taylor(-f.coeffs,f.order)
-(c::Number,f::Taylor)=Taylor(c)-f
-(f::Taylor,c::Number)=f-Taylor(c)

function *(f::Taylor, g::Taylor)
    m=max(f.order,g.order)
        if m==f.order #f,g=Taylor(f,m),Taylor(g,m)
        g=Taylor(g,m)
    else
        f=Taylor(f,m)
    end
    T=promote_type(eltype(f.coeffs),eltype(g.coeffs))
    v=zeros(T,m+1)
    for ord = 0:m
        suma=zero(T)
        for n=0:ord
            suma+=f.coeffs[n+1]*g.coeffs[ord-n+1]
        end
        v[ord+1]=suma
    end
    Taylor(v, m)
end
*(f::Taylor,c::Number)=Taylor(c*f.coeffs, f.order)
*(c::Number,f::Taylor)=Taylor(c*f.coeffs, f.order)

function /(f::Taylor, g::Taylor)
    m=max(f.order,g.order)
    if g.coeffs[1]==0
        return println("División por cero!")
    end
        if m==f.order #f,g=Taylor(f,m),Taylor(g,m)
        g=Taylor(g,m)
    else
        f=Taylor(f,m)
    end
    T=promote_type(eltype(f.coeffs),eltype(g.coeffs))
    v=zeros(T,m+1)
    for i in 0:m
        v[i+1]=(1/g.coeffs[1])*(f.coeffs[i+1]-sum([v[n+1]*g.coeffs[i-n+1] for n in 0:i-1]))
    end
    Taylor(v, m)
end
/(f::Taylor,c::Number)=f*inv(c)
/(c::Number,f::Taylor)=Taylor(c)/f;

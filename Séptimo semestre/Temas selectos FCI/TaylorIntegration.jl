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
    f.coeffs,g.coeffs=promote(f.coeffs,g.coeffs)
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

function *(f::Taylor,G::Array)
  [f*g for g in G]
end

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
/(c::Number,f::Taylor)=Taylor(c)/f
import Base.sqrt
function sqrt(f::Taylor)
    if f.coeffs[1]<0
        error("Raíz negativa")
    end
    m=f.order
    v=zeros(Float64,m+1)
    v[1]=sqrt(f.coeffs[1])
    cte=0.5/v[1]
    for i in 1:m #vamos a ir sacando el coeficiente v[i+1], el cual corresponde a p_i
        if iseven(i) #i es par
            s=(i-2)/2
            suma=sum([v[j+1]*v[i-j+1] for j in 1:s])
            v[i+1]=cte*(f.coeffs[i+1]-2*suma-(v[(i/2)+1]^2))
        else
            s=(i-1)/2
            suma=sum([v[j+1]*v[i-j+1] for j in 1:s])
            v[i+1]=cte*(f.coeffs[i+1]-2*suma)
        end
    end
    Taylor(v,m)
end
function ^(f::Taylor,n::Integer)
  if n==1
    return f
  end
  g=f
  for i in 2:n
    g=g*f
  end
  g
end

using PyPlot

R(i,j,Q)=sqrt((Q[3i-2]-Q[3j-2])^2+(Q[3i-1]-Q[3j-1])^2+(Q[3i]-Q[3j])^2)

R_3(i,j,Q)=1/(R(i,j,Q)^3)

function accel(Q,M,i,j)
    #Esta funcion saca la aceleración que siente j debido a i
    if i==j
        [0,0,0]
    else
        (M[i]*R_3(i,j,Q))*[Q[3i-2]-Q[3j-2],Q[3i-1]-Q[3j-1],Q[3i]-Q[3j]]
    end
end

function total_accel(Q,M,j,n)
    sum([accel(Q,M,s,j) for s in 1:n])
end

function sistema(Q,M,n=length(M))
    A=Any[]
    for i in 1:n
        T=total_accel(Q,M,i,n)
        push!(A,T[1],T[2],T[3])
    end
    A
end

function NextTaylors(f::Function,Q0,V0,M,n)
    coords=3n
    Q=[Float64[Q0[i]] for i in 1:coords]
    V=[Float64[V0[i]] for i in 1:coords]
    for ord in 1:orderMax
        Q0=[Taylor(Q[i]) for i in 1:coords]
        V0=[Taylor(V[i]) for i in 1:coords]
        C1=[f(Q0,M,n)[i].coeffs[ord]/ord for i in 1:coords]
        c2=[V[i][ord]/ord for i in 1:coords]
        for i in 1:coords
            push!(V[i],C1[i])
            push!(Q[i],c2[i])
        end
    end
    Q,V
end

function TimeStep(Q,V)
    MINS_TOT=Float64[]
    m=length(Q)
    for i in 1:m
        MINS_COORD=Float64[]
        pot=inv(orderMax)
        if V[i][end]==0
        else
            push!(MINS_COORD,(epsilon/abs(V[i][end]))^pot)
        end
        if Q[i][end]==0
        else
            push!(MINS_COORD,(epsilon/abs(Q[i][end]))^pot)
        end
        pot=inv(orderMax)-1
        if V[i][end-1]==0
        else
            push!(MINS_COORD,(epsilon/abs(V[i][end-1]))^pot)
        end
        if Q[i][end-1]==0
        else
            push!(MINS_COORD,(epsilon/abs(Q[i][end-1]))^pot)
        end
        if length(MINS_COORD)==0
        else
            push!(MINS_TOT,minimum(MINS_COORD))
        end
    end
    if length(MINS_TOT)==0
        error("Auch")
    else
        minimum(MINS_TOT)
    end
end

function EvalHorner(Q,V,h)
    order=length(Q[1])-1
    n=length(Q)
    Q0=[Float64[]]
    V0=[Float64[]]
    for j in 1:n
        sumaQj=Q[j][order+1]
        sumaVj=V[j][order+1]
        for i in order:-1:1
            sumaQj=Q[j][i]+h*sumaQj
            sumaVj=V[j][i]+h*sumaVj
        end
        push!(Q0,sumaQj)
        push!(V0,sumaVj)
    end
    Q0,V0
end

function TaylorStepper(f::Function,Q0,V0,M,tiempo,n)
    q,v=NextTaylors(f,Q0,V0,M,n)
    h=TimeStep(q,v)
    c1,c2=EvalHorner(q,v,h)
    tiempo+=h
    return c1,c2,tiempo
end

function TaylorIntegrator(Q0,V0,t0,tf,f,M,n=length(M))
    Q=[Float64[Q0[i]] for i in 1:3n]
    V=[Float64[V0[i]] for i in 1:3n]
    T=Float64[t0]
    while t0<tf
        Q0,V0,t0=TaylorStepper(f,Q0,V0,M,t0,n)
        push!(T,t0)
        for i in 1:3n
            push!(Q[i],Q0[i])
            push!(V[i],V0[i])
        end
    end
    T,Q,V
end

function graficadora_cuerpos(q0,v0,t0,tf,M,excluidos)
  n=length(M)
  T,Q,DQ=TaylorIntegrator(q0,v0,t0,tf,sistema,M,n)
    for i in 1:n
        if i in excluidos
        else
            x=3i-2
            y=3i-1
            z=3i
            plot(Q[x],Q[y],"o",label="Planeta $i")
        end
    end
    xlabel(L"$x$")
    ylabel(L"$y$")
    title("Movimiento de los $n cuerpos")
    axis("equal")
    legend(loc=0,fontsize=10)
end;

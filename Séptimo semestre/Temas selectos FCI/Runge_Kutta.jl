using PyPlot

function pars_pot(Q,M,i,j)
    #Esta funcion saca la energía potencial que tiene j debido a i
    if i==j
        0
    else
        -M[i]*M[j]/R(i,j,Q)
    end
end

function total_pot(Q,M,j,n=length(M))
    sum([pars_pot(Q,M,s,j) for s in 1:n])
end

function system_inst_pot(t,Q,V,M,n=length(M))
    A=Float64[]
    for i in 1:n
        T=total_pot(Q,M,i,n)
        push!(A,T)
    end
    sum(A)/2
end

function kinetic(t,Q,V,M,n=length(M))
    V2=[V[3i-2]^2+V[3i-1]^2+V[3i]^2 for i in 1:n]
    K=[0.5*M[i]*V2[i] for i in 1:n]
    sum(K)
end

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

function total_accel(Q,M,j,n=length(M))
    sum([accel(Q,M,s,j) for s in 1:n])
end

function sistema(t,Q,V,M,n=length(M))
    A=Any[]
    for i in 1:n
        T=total_accel(Q,M,i,n)
        push!(A,T[1],T[2],T[3])
    end
    A
end

function rungekutta4(q0,v0,f,t,h,w=1)
    l1=h*f(t,q0,v0,w)
    k1=h*v0
    l2=h*f(t+h/2,q0+k1/2,v0+l1/2,w)
    k2=h*(v0+l1/2)
    l3=h*f(t+h/2,q0+k2/2,v0+l2/2,w)
    k3=h*(v0+l2/2)
    l4=h*f(t+h,q0+k3,v0+l3,w)
    k4=h*(v0+l3)
    qnueva=q0+(1/6)*(k1+2k2+2k3+k4)
    vnueva=v0+(1/6)*(l1+2l2+2l3+l4)
    vnueva,qnueva
end

function integrar(metodo,q0,v0,t0,tf,h,f,M)
    n=length(M)
    q = Vector{Float64}[]
    v = Vector{Float64}[]
    t = Float64[]
    U = Float64[]
    K = Float64[]
    for ti in t0:h:tf
        v0,q0=metodo(q0,v0,f,ti,h,M)
        push!(t,ti)
        push!(q,[q0])
        push!(v,[v0])
        push!(U,system_inst_pot(t0,q0,v0,M,n))
        push!(K,kinetic(t0,q0,v0,M,n))
    end
    E=U+K
    dE=(E-E[1])/eps(E[1])
    t,hcat(q...)',hcat(v...)',U,K,E,dE
end
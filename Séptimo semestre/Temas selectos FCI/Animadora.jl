using PyCall
using PyPlot
@pyimport matplotlib.animation as anim



function limitadora(Q)
    xmin,ymin,zmin=minimos(Q)
    xmax,ymax,zmax=maximos(Q)
    xsup=xmax+((xmax-xmin)/20)
    xinf=xmin-((xmax-xmin)/20)
    ysup=ymax+((ymax-ymin)/20)
    yinf=ymin-((ymax-ymin)/20)
    zsup=zmax+((zmax-zmin)/20)
    zinf=zmin-((zmax-zmin)/20)
    xinf,xsup,yinf,ysup,zinf,zsup
end

function transformer(Q)
    n=length(Q)
    h=length(Q[1])
    Qt=Array{Float64}[]
    for i in 1:h
        push!(Qt,[Q[j][i] for j in 1:n])
    end
    Qt
end

function runge_kuttaer(q,A)
  Q=[q[:,i] for i in [A]]
  a=length(A)
  b=a/2
  Qy=[Q[2i] for i in 1:b]
  Qx=[Q[2i+1] for i in 0:(b-1)]
  Qx,Qy
end

function minimos(Q)
  n=length(Q)
  MIN=Float64[]
  for i in 1:n
    push!(MIN,minimum(Q[i]))
  end
  mini=minimum(MIN)
  if mini<0
    return mini*1.1
  else
    return mini*0.9
  end
end

function maximos(Q)
  n=length(Q)
  MAX=Float64[]
  for i in 1:n
    push!(MAX,maximum(Q[i]))
  end
  maxi=maximum(MAX)
  if maxi>0
    return maxi*1.1
  else
    return maxi*0.9
  end
end

function partidora(Q,j,lim=500)
  inicial=j*lim+1
  maximo=length(Q[1])
  final=inicial+lim-1
  if final>maximo
    final=maximo
  end
  n=length(Q)
  Qred=[Q[i][inicial:final] for i in 1:n]
end


function manual_animation(Qx,Qy,xmin,xmax,ymin,ymax,fmes=100,framesps=30,itvl=20,nombre="Animacion")
  fig = figure()
    ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))

    global line1 = ax[:plot]([], [], "o", lw=2)[1]

    function init()
        global line1

        line1[:set_data]([], [])

        return ([line1],None)
    end

    function animate(i)
        global line1
        # i+1 para que i=0 sea x[1]
        line1[:set_data]([Qx][i+1],[Qy][i+1])
        return ([line1],None)
    end
    name=nombre*".gif"
    myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=fmes, interval=itvl, blit=true);
    myanim[:save](name, writer="imagemagick", fps=framesps)
end


function final_anim(Q,A,nombre="Animacion",frames=length(Q[:,1]),fps=ceil(frames/3),interval=ceil(frames/5))
  Qx,Qy=runge_kuttaer(Q,A)
  limite=500
  xmin,xmax=minimos(Qx),maximos(Qx)
  ymin,ymax=minimos(Qy),maximos(Qy)
  j=ceil(frames/limite)
  if j==1
    Qx,Qy=transformer(Qx),transformer(Qy)
    manual_animation(Qx,Qy,xmin,xmax,ymin,ymax,frames,fps,interval,nombre)
  else
    for i in 1:j
      QX,QY=partidora(Qx,j,limite),partidora(Qy,j,limite)
      Qx,Qy=transformer(QX),transformer(QY)
      manual_animation(Qx,Qy,xmin,xmax,ymin,ymax,frames,fps,interval,nombre*"$i")
    end
  end
end;

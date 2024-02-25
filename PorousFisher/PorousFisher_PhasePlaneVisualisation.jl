using Plots,DifferentialEquations
c=1/sqrt(2);
#c=1
#c=1/2

function diff!(du,u,p,z)
c=p 
du[1]=u[1]*u[2] 
du[2]=-c*u[2]-u[1]*(1-u[1])-u[2]^2 
end

function diffb!(du,u,p,z)
c=p 
du[1]=-u[1]*u[2] 
du[2]=c*u[2]+u[1]*(1-u[1])+u[2]^2 
end
      
function odesolver(c,ξ,ic)
p=(c)
ξspan=(0.0,maximum(ξ))
prob=ODEProblem(diff!,ic,ξspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=ξ);
return sol
end


function odesolverb(c,ξ,ic)
p=(c)
ξspan=(0.0,maximum(ξ))
prob=ODEProblem(diffb!,ic,ξspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=ξ);
return sol
end


ξmax=500
ξ=LinRange(0,ξmax,10000)

ic=[1.1,0.2]
sol=odesolver(c,ξ,ic)
p1=plot(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.9,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.8,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.6,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.4,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.2,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.0,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.3]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.4]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.5]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.6]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.6]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.7]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.8]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-0.9]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-1.0]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-1.1]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[1.0,-1.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,-0.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,-0.4]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,-0.6]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,-0.8]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,-1.0]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[-0.1,-1.2]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.1,-0.7]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ic=[0.1,-0.8]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

ε = 1e-10
UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:orangered,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

UU=1+ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

UU=1+ε
m =(-c-sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

UU=1-ε
m =(-c-sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

UU=-ε
r=(c+sqrt(c^2+4*c))/2
m =-c/r
VV = m*UU - c
ic=[UU,VV]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

UU=-ε
r=(c-sqrt(c^2+4*c))/2
m =-c/r
VV = m*UU - c
ic=[UU,VV]
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))


UU=ε
r=(c+sqrt(c^2+4*c))/2
m =-c/r
VV = m*UU - c
ic=[UU,VV]
sol=odesolver(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2))

UU=ε
r=(c-sqrt(c^2+4*c))/2
m =-c/r
VV = m*UU - c
ic=[UU,VV]
sol=odesolverb(c,ξ,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.3, 0.2),xlabel="U",ylabel="V")



p1=scatter!([1,0,0],[0,0,-c],color=:grey0,markersize=10,msw=0)
display(p1)
#savefig(p1,"Figure.pdf")

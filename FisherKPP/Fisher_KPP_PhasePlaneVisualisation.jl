using Plots,DifferentialEquations
c=2;

function diff!(du,u,p,z)
c=p 
du[1]=u[2] 
du[2]=-c*u[2]-u[1]*(1-u[1]) 
end


function diffb!(du,u,p,z)
c=p 
du[1]=-u[2] 
du[2]=c*u[2]+u[1]*(1-u[1]) 
end
      
function odesolver(c,Z,ic)
p=(c)
zspan=(0.0,maximum(Z))
prob=ODEProblem(diff!,ic,zspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=Z);
return sol
end


function odesolverb(c,Z,ic)
p=(c)
zspan=(0.0,maximum(Z))
prob=ODEProblem(diffb!,ic,zspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=Z);
return sol
end


Zmax=200;
Z=LinRange(0,Zmax,10000);

ic=[1.1,0.2]
sol=odesolver(c,Z,ic)
p1=plot(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[1.0,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-0.3, 0.2))

ic=[0.9,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.8,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.6,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.4,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.2,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.0,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[-0.1,0.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))


ic=[1.2,-0.75]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))


ic=[1.2,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[1.0,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.8,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.6,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.4,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.2,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[0.0,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ic=[-0.1,-1.2]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:blue,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

ε = 1e-10
UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:orangered,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

UU=1+ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolver(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

UU=1+ε
m =(-c-sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolverb(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2))

UU=1-ε
m =(-c-sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolverb(c,Z,ic)
p1=plot!(sol[1,:],sol[2,:],lw=2,color=:gold,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.2),xlabel="U",ylabel="V")

p1=scatter!([1,0],[0,0],color=:grey0,markersize=10,msw=0)
display(p1)
savefig(p1,"Figure.pdf")
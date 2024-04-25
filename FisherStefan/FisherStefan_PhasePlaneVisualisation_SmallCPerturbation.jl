using Plots,DifferentialEquations
c=1/10

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


Zmax=100
Z=LinRange(0,Zmax,10000)    
ε = 1e-8
UU=1+ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolver(c,Z,ic)
ii = findfirst(z -> z > 1.5 , sol[1,:])
p1=plot(sol[1,1:ii],sol[2,1:ii],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.5))


UU=1+ε
m =(-c-sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolverb(c,Z,ic)
ii = findfirst(z -> z > 1.5 , sol[1,:])
p1=plot!(sol[1,1:ii],sol[2,1:ii],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.5))

UU=1-ε
m =(-c-sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolverb(c,Z,ic)
ii = findfirst(z -> z <-0.2 , sol[1,:])
p1=plot!(sol[1,1:ii],sol[2,1:ii],color=:gold,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.5),xlabel="U",ylabel="V")


UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
sol=odesolver(c,Z,ic)
ii = findfirst(z -> z < 0.0 , sol[1,:])
p1=plot!(sol[1,1:ii-1],sol[2,1:ii-1],color=:orangered,lw=2,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.5))
#p1=plot!(sol[1,ii:end],sol[2,ii:end],color=:orangered,lw=2,ls=:dash,legend=false,xlims=(-0.2, 1.2),ylims=(-1.2, 0.5))

g = (sol[2,ii]- sol[2,ii-1])/(sol[1,ii]- sol[1,ii-1])
Vdagger = sol[2,ii-1]+g*(sol[1,ii-1])
κ=-c/Vdagger

p1=scatter!([1,0],[0,0],color=:grey0,markersize=10,msw=0)
p1=scatter!([0],[Vdagger],color=:green,markersize=10,msw=0,markershape=:square,xlabel="U",ylabel="V")



V0(U)=-1*sqrt((2*U^3-3*U^2+1)/3)
V1(U)=-1*(3*sqrt(3)+sqrt(2*U+1)*(2*U^2-3*U-2))/(5*sqrt(2*U+1)*(U-1))
P(U)=V0(U)+c*V1(U)
p1=plot!(P,0,1,color=:green,ls=:dash,lw=2,legend=false)
display(p1)
#savefig(p1,"Phaseplanecm05.pdf")



Uplot=zeros(ii)
Uplot[1]=1
for i in 2:ii
Uplot[i] = Uplot[i-1]+(Z[2]-Z[1])*(sol[2,i])
end
Uplotp=zeros(ii)
Uplotp[1]=1
for i in 2:ii
Uplotp[i] = Uplotp[i-1]+(Z[2]-Z[1])*(P(sol[1,i]))
end

p2=plot(Z[1:ii].-Z[ii],Uplot,lw=2,color=:blue,legend=false,xlabel="z",ylabel="U",ylims=(0,1.1),xlims=(-10,0))
p2=plot!(Z[1:ii].-Z[ii],Uplotp,lw=2,color=:red,ls=:dash,legend=false,xlabel="z",ylabel="U",ylims=(0,1.1),xlims=(-10,0))

p3=plot(p1,p2,layout=(2,1))
#savefig(p3,"Figure.pdf")
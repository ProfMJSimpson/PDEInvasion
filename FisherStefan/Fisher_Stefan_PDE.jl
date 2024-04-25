using Plots, DifferentialEquations, Polynomials

function diff!(du,u,p,t)
dρ,N,κ,ρ =p

du[1]=(u[2]-u[1])/(u[N+1]*dρ)^2+u[1]*(1-u[1]) 
for i in 2:N-1
du[i]=(-ρ[i]*κ*(-4*u[N-1]+u[N-2])*(u[i+1]-u[i-1]))/(4*dρ^2*u[N+1]^2)+(u[i+1]-2*u[i]+u[i-1])/(u[N+1]*dρ)^2+u[i]*(1-u[i]) 
end
du[N]=0 
du[N+1]=-κ*(-4*u[N-1]+u[N-2])/(2*u[N+1]*dρ)    
end
    
    
function pdesolver(dρ,N,T,u0,κ,ρ)
p=(dρ,N,κ,ρ)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,u0,tspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=T);
return sol
end

function diffpp!(du,u,p,z)
c=p 
du[1]=u[2] 
du[2]=-c*u[2]-u[1]*(1-u[1]) 
end
        
        
function diffppb!(du,u,p,z)
c=p 
du[1]=-u[2] 
du[2]=c*u[2]+u[1]*(1-u[1]) 
end
              
function odesolver(c,Z,ic)
p=(c)
zspan=(0.0,maximum(Z))
prob=ODEProblem(diffpp!,ic,zspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=Z);
return sol
end
        
        
function odesolverb(c,Z,ic)
p=(c)
zspan=(0.0,maximum(Z))
prob=ODEProblem(diffbpp!,ic,zspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=Z);
return sol
end

L=100
T=0:1:50
N=5001
dρ=1/(N-1)
κ=1

u0=zeros(N+1)
ρ=zeros(N)
u0[N+1]=L


for i in 1:N
ρ[i] = (i-1)*dρ
end

#This initial condition is a step function 
for i in 1:N-1
    if ρ[i] <= 1.0
    u0[i] = 0.5
    end
end

@time sol=pdesolver(dρ,N,T,u0,κ,ρ); 
LL=zeros(length(T))
LL[:]=sol[N+1,:]
p0=scatter(T,LL,legend=false)
c=-κ*(-4*sol[N-1,end]+sol[N-2,end])/(2*LL[end]*dρ) 

Zmax=100
Z=LinRange(0,Zmax,10000)
ε = 1e-10
UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
solpp=odesolver(c,Z,ic)
ii = findfirst(z -> z < 0.0 , solpp[1,:])
Z_shift = LL[end]-Z[ii]
ZS=Z.+Z_shift

p1=plot(ρ*LL[1],sol[1:1:N,1],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ρ*LL[2],sol[1:1:N,2],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ρ*LL[11],sol[1:1:N,11],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ρ*LL[21],sol[1:1:N,21],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ρ*LL[31],sol[1:1:N,31],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ρ*LL[41],sol[1:1:N,41],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ρ*LL[51],sol[1:1:N,51],lw=2,color=:blue,xlims=(50,150),ylims=(0,1.1),legend=false)
p1=plot!(ZS[1:ii],solpp[1,1:ii],lw=2,color=:red,ls=:dash,xlims=(50,150),ylims=(0,1.1),legend=false,xlabel="x",ylabel="u(x,t)")
p3=plot(p1,p0,layout=(2,1))
display(p3)
savefig(p3, "Figure.pdf")
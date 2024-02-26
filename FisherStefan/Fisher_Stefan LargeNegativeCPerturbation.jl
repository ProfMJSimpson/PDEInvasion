using Plots, DifferentialEquations, Polynomials


function diff!(du,u,p,t)
dξ,N,κ,ξ=p

du[1]=(u[2]-u[1])/(u[N+1]*dξ)^2+u[1]*(1-u[1]) 
for i in 2:N-1
du[i]=(-ξ[i]*κ*(-4*u[N-1]+u[N-2])*(u[i+1]-u[i-1]))/(4*dξ^2*u[N+1]^2)+(u[i+1]-2*u[i]+u[i-1])/(u[N+1]*dξ)^2+u[i]*(1-u[i]) 
end
du[N]=0 
du[N+1]=-κ*(-4*u[N-1]+u[N-2])/(2*u[N+1]*dξ)    
end
    
    
function pdesolver(dξ,N,T,u0,κ,ξ)
p=(dξ,N,κ,ξ)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,u0,tspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=T);
return sol
end

L=200
T=0:1:50
N=5001
dξ=1/(N-1)
κ=-90/100

u0=zeros(N+1)
ξ=zeros(N)
u0[N+1]=L

for i in 1:N
ξ[i] = (i-1)*dξ
end

#This initial condition is a step function 
for i in 1:N-1
    if ξ[i] <= 1.0
    u0[i] = 0.5
    end
end

sol=pdesolver(dξ,N,T,u0,κ,ξ); 
LL=zeros(length(T))
LL[:]=sol[N+1,:]
p0=scatter(T,LL,legend=false)
c=-κ*(-4*sol[N-1,end]+sol[N-2,end])/(2*LL[end]*dξ) 

#Perturbation solution
U0(z)=1-exp(-z*c)
U1(z)=(exp(-2*z*c)+(2*z*c-1)*exp(-z*c))/2
P(z)=U0(z)+U1(z)/c^2

p1=plot(ξ*LL[1],sol[1:1:N,1],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[2],sol[1:1:N,2],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[11],sol[1:1:N,11],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[21],sol[1:1:N,21],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[31],sol[1:1:N,31],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[41],sol[1:1:N,41],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[51],sol[1:1:N,51],lw=2,color=:blue,xlims=(50,250),ylims=(0,1.1),legend=false,xlabel="x",ylabel="u(x,t)")
p2=plot(ξ*LL[end].-LL[end],sol[1:1:N,end],lw=2,xlims=(-10,0),ylims=(0,1.1),legend=false)
p2=plot!(P,-10,0,lw=2,ls=:dash,color=:red,legend=false,xlabel="z",ylabel="U")

p3=plot(p1,p2,layout=(2,1))
#savefig(p3,"Figure.pdf")
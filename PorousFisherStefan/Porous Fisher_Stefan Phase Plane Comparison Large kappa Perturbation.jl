using Plots, DifferentialEquations, Polynomials,QuadGK

function diff!(dϕ,ϕ,p,t)
dξ,N,κ,ξ=p
dLdt=-κ*(3*ϕ[N]-4*ϕ[N-1]+ϕ[N-2])/(2*2*ϕ[N+1]*dξ)
Φ=max(0,ϕ[1])
dϕ[1]=sqrt(Φ)*(ϕ[2]-ϕ[1])/(ϕ[N+1]*dξ)^2+2*ϕ[1]*(1-sqrt(Φ))
for i in 2:N-1
Φ=max(0,ϕ[i])
dϕ[i]=ξ[i]*dLdt*(ϕ[i+1]-ϕ[i-1])/(2*dξ*ϕ[N+1])+sqrt(Φ)*(ϕ[i+1]-2*ϕ[i]+ϕ[i-1])/(ϕ[N+1]*dξ)^2+2*ϕ[i]*(1-sqrt(Φ)) 
end
dϕ[N]=0 
dϕ[N+1]=dLdt   
end
    
    
function pdesolver(dξ,N,T,ϕ0,κ,ξ)
p=(dξ,N,κ,ξ)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,ϕ0,tspan,p)
sol=solve(prob,saveat=T);
return sol
end


        

L=200
T=0:1:50
N=5001
dξ=1/(N-1)
D=1.0 
λ=1.0
κ=10.0

ϕ0=zeros(N+1)
ξ=zeros(N)
ϕ0[N+1]=L

#This initial condition is expmentnailly decay
for i in 1:N
ξ[i] = (i-1)*dξ
end

#This initial condition is a step function 
for i in 1:N-1
    ϕ0[i] = 0.5^2
end
ϕ0[N]=1e-8

@time sol=pdesolver(dξ,N,T,ϕ0,κ,ξ); 
LL=zeros(length(T))
LL[:]=sol[N+1,:]
p0=scatter(T,LL,legend=false)
QQ=Int(round(0.10*length(T)))
Tfit=zeros(QQ)
xfit=zeros(QQ)
for k in 1:QQ
Tfit[k]=T[length(T)-QQ+k]
xfit[k]=LL[length(T)-QQ+k]
end
p=Polynomials.fit(Tfit,xfit, 1)
f(t) = p[1]*t + p[0]
ws = p[1]
p0=plot!(f,0,T[end],c=:red,lw=3,label=false)
display(p0)
ws = p[1]
c=-κ*(-4*sol[N-1,end]+sol[N-2,end])/(2*2*LL[end]*dξ)


ε = 1.e-10
UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
p=(c)
ode(ψ, p, ϕ) = -c/(sqrt(ϕ))-2*(sqrt(ϕ)-ϕ)/ψ
ψ0 = VV
ϕspan = (UU,1e-10)
prob = ODEProblem(ode, ψ0, ϕspan)
Φ=LinRange(UU,1e-10,1000)
alg=Tsit5()
solp=solve(prob,saveat = Φ)
q1=plot(Φ,solp,lw=2,xlabel="ϕ",ylabel="ψ",legend=false,xlims=(-0.2, 1.2),ylims=(-1.5, 0.5))
ee = 1/sqrt(2)-c
V0(x)=-sqrt(2*x)*(1-sqrt(x))
V1(x)=-2*(1-sqrt(x))/3
P(x) = V0(x)+ee*V1(x)
q1=plot!(P,0,1,lw=2,ls=:dash,legend=false)
display(q1)





dΦ=Φ[2]-Φ[1]
Z=zeros(length(Φ))
Zp=zeros(length(Φ))
Z[end]=0.0
Zp[end]=0.0
for i in length(Z)-1:-1:1
Z[i] = Z[i+1]-dΦ/solp[i]
Zp[i] = Zp[i+1]-dΦ/P(Φ[i])
end
s1=plot(Z,sqrt.(Φ),lw=2,xlims=(-10,0))
s1=plot!(Zp,sqrt.(Φ),lw=2,color=:red,ls=:dash,xlims=(-10,0),legend=false,xlabel="x",ylabel="u")







p1=plot(ξ*LL[1],sqrt.(sol[1:1:N,1]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[2],sqrt.(sol[1:1:N,2]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[11],sqrt.(sol[1:1:N,11]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[21],sqrt.(sol[1:1:N,21]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[31],sqrt.(sol[1:1:N,31]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[41],sqrt.(sol[1:1:N,41]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(ξ*LL[51],sqrt.(sol[1:1:N,51]),lw=2,color=:blue,xlims=(0,250),ylims=(0,1.1),legend=false)
p1=plot!(Z.+LL[51],sqrt.(Φ),lw=2,color=:red,ls=:dash,xlims=(0,250),ylims=(0,1.1),legend=false,xlabel="x",ylabel="u(x,t)")
p3=plot(p1,p0,layout=(2,1))
display(p3)
savefig(p3, "Figure.pdf")
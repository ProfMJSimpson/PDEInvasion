using Plots, DifferentialEquations, Polynomials
function diff!(du,u,p,t)
dx,N=p
for i in 2:N-1
Di=u[i] #Diffusivity at mesh point i
Dip=u[i+1] #Diffusivity at mesh point i+1
Dim=u[i-1] #Diffusivity at mesh point i-1
du[i]=((Di+Dip)*(u[i+1]-u[i])-(Di+Dim)*(u[i]-u[i-1]))/(2*dx^2)+u[i]*(1-u[i]) #discretised equation at mesh point i
end
i=1
Di=u[i]
Dip=u[i+1]
du[i]=((Di+Dip)*(u[i+1]-u[i]))/(2*dx^2)+u[i]*(1-u[i])
i=N
Di=u[i]
Dim=u[i-1]
du[i]=((Di+Dim)*(u[i-1]-u[i]))/(2*dx^2)+u[i]*(1-u[i])
end
    
    
function pdesolver(L,dx,N,T,u0)
p=(dx,N)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,u0,tspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=T);
return sol
end

function diffpp!(du,u,p,z)
c=p 
du[1]=u[1]*u[2] 
du[2]=-c*u[2]-u[1]*(1-u[1])-u[2]^2 
end
    
   
function odesolver(c,ξ,ic)
p=(c)
ξspan=(0.0,maximum(ξ))
prob=ODEProblem(diffpp!,ic,ξspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=ξ);
return sol
end

L=300 #length of domain
T=0:1:50 #Vector of values of time that the solution will be computed
xpos=zeros(length(T)) #Vector to store the front position
dx=0.01 #discretisation
N=Int(L/dx)+1 #number of mesh points
a=1/5 #Exponential decrease for the initial condition
#a=1/3 #Exponential decrease for the initial condition
#a=1 #Exponential decrease for the initial condition
#a=2 #Exponential decrease for the initial condition
X=10 #Constant for the initial condition
UU0=0.6 #Constant for the initial condition 

#Set up the initial condition
u0=zeros(N);
x=zeros(N); 
for i in 1:N
x[i] = (i-1)*dx
if x[i] <= X
   u0[i] = UU0
else
    u0[i] = UU0*exp(-a*(x[i]-X))
    end
end


#Solve the pde
sol=pdesolver(L,dx,N,T,u0);



for k in 2:length(T)
ii = findfirst(z -> z < 0.5 , sol[:,k])
g = (sol[ii,k]- sol[ii-1,k])/dx
xpos[k] = (sol[ii,k]-0.5)/g + x[ii]
end

QQ=Int(round(0.10*length(T)))
Tfit=zeros(QQ)
xfit=zeros(QQ)
for k in 1:QQ
Tfit[k]=T[length(T)-QQ+k]
xfit[k]=xpos[length(T)-QQ+k]
end
p=Polynomials.fit(Tfit,xfit, 1)
f(t) = p[1]*t + p[0]
ws = p[1]

c=ws
ξ=LinRange(0,2000,10000)
Z=zeros(length(ξ))
ZZ = xpos[end]
ε = 1e-10
UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]
solp=odesolver(c,ξ,ic)
Z[1]=20

for i in 2:length(ξ)
Z[i] = Z[i-1]+ 0.5*(ξ[2]-ξ[1])*(solp[1,i]+solp[1,i-1])
end
ii = findfirst(z -> z < 0.5 , solp[1,:])
g = (solp[1,ii]- solp[1,ii-1])/(Z[ii]-Z[ii-1])
Zt = (solp[1,ii]-0.5)/g + Z[ii]
shift = ZZ - Zt

p1=plot(x,sol[:,1],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,2],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,11],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,21],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,31],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,41],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,51],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1),xlabel="x",ylabel="u(x,t)")
p1=plot!(Z.+shift,solp[1,:],lw=2,color=:red,linestyle=:dash,legend=false,xlims=(0,L),ylims=(0,1.1))

p2=scatter(T,xpos, markersize=4,markercolour=:blue,msw=0,legend=false)
p2=plot!(f,0,maximum(T),lw=2,color=:magenta1,xlabel="t",ylabel="X(t)")


p3=plot(p1,p2,layout=(2,1))
savefig(p3, "Figure.pdf")
#Time-dependent PDE solutions for Fisher-KPP
using Plots, DifferentialEquations, Polynomials

function diff!(du,u,p,t)
dx,N=p         #Parameters 
    for i in 2:N-1 #Semi-discrete equations
    du[i]=(u[i+1]-2*u[i]+u[i-1])/dx^2+u[i]*(1-u[i]) 
    end

du[1]=(u[2]-u[1])/dx^2+u[1]*(1-u[1]) #Boundary Condition 
du[N]=(u[N-1]-u[N])/dx^2+u[N]*(1-u[N]) #Boundary Condition 
end
    
    
function pdesolver(L,dx,N,T,u0)
p=(dx,N)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,u0,tspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=T);
return sol
end

L=300 #length of domain
T=0:1:50 #Vector of values of time that the solution will be computed
xpos=zeros(length(T)) #Vector to store the front position
dx=0.1 #discretisation
N=Int(L/dx)+1 #number of mesh points
a=(5-sqrt(21))/2 #Exponential decrease for the initial condition
#a=(3-sqrt(5))/2 #Exponential decrease for the initial condition
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
sol=pdesolver(L,dx,N,T,u0)


for k in 2:length(T)
ii = findfirst(z -> z < 0.5 , sol[:,k]) #Find the first value of u(x,t) where u < 0.5
g = (sol[ii,k]- sol[ii-1,k])/dx #Linearly interpolate 
xpos[k] = (sol[ii,k]-0.5)/g + x[ii]
end


QQ=Int(round(0.10*length(T))) #Use the last 10% of time points to estimate the front speed
Tfit=zeros(QQ)
xfit=zeros(QQ)
    for k in 1:QQ
    Tfit[k]=T[length(T)-QQ+k]
    xfit[k]=xpos[length(T)-QQ+k]
    end
p=Polynomials.fit(Tfit,xfit, 1)
f(t) = p[1]*t + p[0]
ws = p[1] #Estimate of wave speed

#Plot the front position and linear regression of X(t)
p2=scatter(T,xpos, markersize=4,markercolour=:blue,msw=0,legend=false)
p2=plot!(f,0,maximum(T),lw=2,color=:magenta1,xlabel="t",ylabel="X(t)")


function diffpp!(du,u,p,z) #ODEs defining the phase plane dynamical system
c=p 
du[1]=u[2] 
du[2]=-c*u[2]-u[1]*(1-u[1]) 
end
    
function odesolver(c,Z,ic)
p=c
zspan=(0.0,maximum(Z))
prob=ODEProblem(diffpp!,ic,zspan,p)
alg=Tsit5()
sol=solve(prob,alg,saveat=Z);
return sol
end

c=ws
Z=LinRange(0,200,1000)
dZ = Z[2]-Z[1]
ε=1e-10 #Start on the unstable manifold close to the equilibrium point
UU=1-ε
m =(-c+sqrt(c^2+4))/2
VV = m*(UU-1)
ic=[UU,VV]  #Start on the unstable manifold close to the equilibrium point
solp=odesolver(c,Z,ic)
ii = findfirst(z -> z < 0.5 , solp[1,:])
g = (solp[1,ii]- solp[1,ii-1])/dZ
Z05 = (solp[1,ii]-0.5)/g + Z[ii]
Z_shift = xpos[end] - Z05 #Interpolate and shift

#Plot and save outcome
p1=plot(x,sol[:,1],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,2],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,11],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,21],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,31],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,41],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,51],lw=2,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1),xlabel="x",ylabel="u(x,t)")
p1=plot!(Z.+Z_shift,solp[1,:],lw=2,color=:red,linestyle=:dash,legend=false,xlims=(0,L),ylims=(0,1.1))
p3=plot(p1,p2,layout=(2,1))
display(p3)
savefig(p3, "Figure.pdf")
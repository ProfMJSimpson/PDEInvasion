using Plots, DifferentialEquations, Polynomials
function diff!(du,u,p,t)
    dx,N=p
    for i in 2:N-1
    Di=u[i]
    Dip=u[i+1]
    Dim=u[i-1]
    du[i]=((Di+Dip)*(u[i+1]-u[i])-(Di+Dim)*(u[i]-u[i-1]))/(2*dx^2)+u[i]*(1-u[i])
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

L=500 #length of domain
T=0:1:50 #Vector of values of time that the solution will be computed
xpos=zeros(length(T)) #Vector to store the front position
dx=0.02 #discretisation
N=Int(L/dx)+1 #number of mesh points
a=1/1 #Exponential decrease for the initial condition
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

p1=plot(x,sol[:,1],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,2],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,11],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,21],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,31],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,41],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1))
p1=plot!(x,sol[:,51],lw=4,color=:blue,legend=false,xlims=(0,L),ylims=(0,1.1),xlabel="x",ylabel="u(x,t)")


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

p2=scatter(T,xpos, markersize=4,markercolour=:blue,msw=0,legend=false)
p2=plot!(f,0,maximum(T),lw=2,color=:magenta1,xlabel="t",ylabel="X(t)")



c=ws
U0(z)=1/(1+exp(z/c))
U1(z)=exp(z/c)*(log(2)-3/2+(3+(1+exp(z/c))*(z/c-log(1+exp(z/c))))/(1+exp(z/c)) )/(1+exp(z/c))^2
P(z) = U0(z)+(c)^(-2)*U1(z)


Lmax=20;
p2=plot(x.-xpos[end],sol[:,51],lw=2,color=:blue,legend=false,xlims=(-Lmax,Lmax),ylims=(0,1.1))
p2=plot!(U0,-Lmax,Lmax,lw=2,color=:tomato,ls=:dash,legend=false,xlims=(-Lmax,Lmax),ylims=(0,1.1),xlabel="z",ylabel="U(z)")
p3=plot(x.-xpos[end],sol[:,51],lw=2,color=:blue,legend=false,xlims=(-Lmax,Lmax),ylims=(0,1.1))
p3=plot!(P,-Lmax,Lmax,lw=2,color=:tomato,ls=:dash,legend=false,xlims=(-Lmax,Lmax),ylims=(0,1.1),xlabel="z",ylabel="U(z)")
p4=plot(p2,p3,layout=(2,1))
savefig(p4, "Figure.pdf")
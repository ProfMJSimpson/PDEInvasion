using Plots, DifferentialEquations
ε = 0.5
c = 1/sqrt(2) - ε
ε1 = 1.e-10
UU=1-ε1
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
V0(x)=-sqrt(2*x)*(1-sqrt(x))
V1(x)=-2*(1-sqrt(x))/3
P(x) = V0(x)+ε*V1(x)
q1=plot!(P,0,1,lw=2,ls=:dash,legend=false)
q1=scatter!([1,0],[0,0],color=:grey0,markersize=10,msw=0)
q1=scatter!([0],[solp[end]],color=:green,markersize=10,msw=0,markershape=:square)


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
s1=plot!(Zp,sqrt.(Φ),lw=2,ls=:dash,color=:red,xlims=(-10,0),legend=false,xlabel="z",ylabel="U")


p3=plot(q1,s1,layout=(2,1))
display(p3)
savefig(p3, "Figure.pdf")
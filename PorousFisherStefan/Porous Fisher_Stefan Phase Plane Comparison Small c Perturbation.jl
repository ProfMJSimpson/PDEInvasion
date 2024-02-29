using Plots, DifferentialEquations, Polynomials,QuadGK

c=0.5
ε = 1.e-10
UU=1-ε
m =(-c+sqrt(c^2+4)/2)
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
V0(x)=-1*sqrt(2*x^2+2/3-8*x^(3/2)/3)
W(x) = quadgk(s->sqrt(2/(3*s)-8*sqrt(s)/3+2*s),x,(1-1e-10))[1]
V1(x)=-1*W(x)/V0(x)
P(x) = V0(x)+c*V1(x)
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
q2=plot(Z,sqrt.(Φ),xlims=(-10,0),lw=2,legend=false, xlabel="z",ylabel="U")
q2=plot!(Zp,sqrt.(Φ),xlims=(-10,0),lw=2,ls=:dash)

q3=plot(q1,q2,layout=(2,1))
savefig(q3, "Figure.pdf")


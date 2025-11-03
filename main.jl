include("Newton-Lagrange.jl"); include("PQS_local_igualdade.jl") # implementação dos algoritmos 2 e 3 respectivamente
include("PQS_local.jl"); include("PQS_globalizado.jl") # implementação dos algoritmos 4 e 5 respectivamente
include("problems_TCC.jl") # problemas de otimização do estudo numérico da monografia

# Problemas com restrições de igualdade e implementação dos algoritmos 2 e 3
# f, h, gradf, Jacobian_h, hessf, hess_coordenada_h, x0, λ0 = PtReq2()
# xeq1, λeq1, keq1 = metodo_newton_lagrange(h, gradf, Jacobian_h, hessf, hess_coordenada_h, x0, λ0;maxiter=100,tol=1.e-6)
# xeq2, λeq2, keq2 = PQS_local_equal(h, gradf, Jacobian_h, hessf, hess_coordenada_h, x0, λ0; maxiter=100, tol=1.e-6, fat=2)

# Problemas com restrições de igualdade e desigualdade e implementação dos algoritmos 4 e 5
f, g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0=PaI()
#xineq1, λineq1, μineq1, kineq1 = PQS_local(g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0; maxiter=100, tol=1.e-6, fat=2)
 xineq2, λineq2, μineq2, kineq2 = PQS_global(g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0; maxiter=100, tol=1.e-6, fat=2.0, cbar=1.0, θ=0.5, σ=0.1)
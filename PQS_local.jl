using JuMP, Clarabel, LinearAlgebra, Printf

include("MetodoGradiente.jl")

function PQS_local(g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0; maxiter=100, tol=1.e-6, fat=2)
    xk = copy(x0) # Usar copy para evitar modificar o x0 original
    λk = copy(λ0)
    μk = copy(μ0)
    n = length(x0)
    l = length(λ0)
    m = length(μ0)
    p=1.e-4
    
    # Gradiente da Lagrangiana aumentada (para o critério de parada)
    gradL(x, λ, μ; c=1) = [gradf(x)+Jacobian_h(x)'*(λ.+c*h(x))+sum(max(0,c*g(x)[i]+μ[i])*grad_coordenada_g[i](x) for i in 1:m); h(x);[1/c*(max(0,c*g(x)[i]+μ[i])-μ[i]) for i in 1:m]]
    
    # Hessiana da Lagrangiana aumentada do problema com restrições de igualdade em relação a x (para fornecer matrizes definidas positivas)
    hess_L(x, λ, μ; c=p)=hessf(x)+sum(λ[i]*hess_coordenada_h[i](x) for i in 1:l)+c*Jacobian_h(x)'*Jacobian_h(x)+sum(c*grad_coordenada_g[i](x)*grad_coordenada_g[i](x)'+max(0,c*g(x)[i]+μ[i])*hess_coordenada_g[i](x) for i in 1:m)
    
    println("Iniciando PQS local...")
    println("Iter | Norma do Grad. Lagrangiano")
    println("------------------------------------")

    for iter in 0:maxiter
        norm_gradL = norm(gradL(xk, λk, μk))
        println(@sprintf("%4d | %.6e", iter, norm_gradL))

        if norm_gradL < tol
            println("Convergência alcançada!")
            return xk, λk, μk, iter
        end

        # 1. Construir o modelo do subproblema quadrático no passo 'd'
        Hk = hess_L(xk, λk, μk)

        if isposdef(Hk)==false
            for i in 1:maxiter 
                p=p*fat
                Hk=hess_L(xk, λk, μk,c=p)
                if isposdef(Hk)==true
                    p=1.e-4
                    break
                end
            end
            if isposdef(Hk)==false
                println("Não foi possível encontrar uma matriz definida positiva na iteração $iter.")
                return xk, λk, μk, -1
            end
        end

        g_fk = gradf(xk)
        J_hk = Jacobian_h(xk)
        hk = h(xk)
        J_gk = Jacobian_g(xk)
        gk = g(xk)
        
        model = Model(Clarabel.Optimizer)
        set_silent(model)
        
        @variable(model, d[1:n])
        
        # Objetivo é uma aproximação quadrática da Lagrangiana
        @objective(model, Min, 0.5 * d' * Hk * d + g_fk' * d)
        
        # Linearização das restrições de igualdade dadas por h(x) = 0
        @constraint(model, c1, J_hk * d .== -hk)
        
        # Linearização das restrições de desigualdade dadas por g(x) ≤ 0
        @constraint(model, c2, J_gk * d .<= -gk)
        
        optimize!(model)
        
        if !is_solved_and_feasible(model)
            println("O subproblema não pôde ser resolvido na iteração $iter.")
            return xk, λk, μk, -1
        end

        # 2. Obter o passo 'd' e os novos multiplicadores
        dk = value.(d)

        λ_new = -dual.(c1)
        μ_new = -dual.(c2)
        
        # 3. Atualizar x e λ
        xk += dk
        λk = λ_new
        μk = μ_new
    end

    println("Número máximo de iterações atingido.")
    return xk, λk, μk, -2
end
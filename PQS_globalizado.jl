using JuMP, Clarabel, LinearAlgebra, Printf

# Referência: Otimização - volume 2 Métodos Computacionais por Alexey Izmailov e Mikhail Solodov

# A função PQS_global e uma implementação da Programação Quadrática Sequêncial com uma estratégia de globalização que emprega busca linear e uma função mérito

# Parâmetros Obrigatórios:
#------------------------------------------------------------------------------------------------------------------------------------------------
# g (Function) - Função das restrições de desigualdade
# h (Function) - Função das restrições de igualdade
# gradf (Function) - Gradiente da função objetivo
# grad_coordenada_g (Vector) - Vetor contendo os gradientes de cada restrição de desigualdade
# Jacobian_h (Function) - Jacobiana das restrições de igualdade
# Jacobian_g (Function) - Jacobiana das restrições de desigualdade
# hessf (Function) - Hessiana da função objetivo
# hess_coordenada_h (Function) - hessiana das restrições de igualdade
# hess_coordenada_g (Function) - hessiana das restrições de desigualdade
# x0 (Vector{Float64}) - Chute inicial para o minimizador
# λ0 (Vector{Float64}) - Chute inicial para o vetor dos multiplicadores de Lagrange referentes as restrições de igualdade
# μ0 (Vector{Float64}) - Chute inicial para o vetor dos multiplicadores de Lagrange referentes as restrições de desigualdade

# Parâmetros Opcionais:
#-----------------------------------------------------------------------------------------------------------------------------------------------
# maxiter (Iter) - Número máximo de iteradas (padrão: 100 iteradas)
# tol (Float64) - tolerância utilizada no critério de parada (padrão: 1.e-6)
# fat (Float64) - fator utilizado para o aumento do parâmetro c da Lagrangiana aumentada (padrão: 2.0)
# cbar (Float64) - parâmetro utilizado no calculo de ck da função mérito (padrão: 2.0)
# θ (Float 64) - parâmetro utilizado no processo de backtracking do parâmetro de comprimento de passo na busca linear, onde 0<θ<1 (padrão 0.5)
# σ (Float 64) - parâmetro utilizado na busca de Armijo, onde 0<θ<1 (padrão 0.1)

function PQS_global(g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0; maxiter=100, tol=1.e-6, fat=2.0, cbar=1.0, θ=0.5, σ=0.1)
    xk = copy(x0) # Usar copy para evitar modificar o x0 original
    λk = copy(λ0)
    μk = copy(μ0)
    ck = 0
    Δk = 0
    n = length(x0)
    l = length(λ0)
    m = length(μ0)
    p=1.e-4
    
    # Hessiana da Lagrangiana aumentada do problema com restrições de igualdade em relação a x (para fornecer matrizes definidas positivas)
    hess_L(x, λ, μ; c=p)=hessf(x)+sum(λ[i]*hess_coordenada_h[i](x) for i in 1:l)+c*Jacobian_h(x)'*Jacobian_h(x)+sum(c*grad_coordenada_g[i](x)*grad_coordenada_g[i](x)'+max(0,c*g(x)[i]+μ[i])*hess_coordenada_g[i](x) for i in 1:m)

    # Função de penalização externa exata
    ψ(x)=norm(h(x),1)+sum(max(0,g(x)[i]) for i in 1:m)

    # Função mérito
    ϕ(x,c)=f(x)+c*ψ(x)
    
    println("Iniciando PQS globalizado...")
    println("Iter | Norma da dir. de descida")
    println("------------------------------------")

    for iter in 0:maxiter

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
                #println("Não foi possível encontrar uma matriz definida positiva na iteração $iter.")
                #return xk, λk, μk, -1
                Hk=I(n)
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

        norm_d = norm(dk)
        println(@sprintf("%4d | %.6e", iter, norm_d))

        if norm_d < tol
            println("Convergência alcançada!")
            return xk, λk, μk, iter
        end

        ck = cbar+norm([λ_new; μ_new], Inf)
        Δk = g_fk'*dk-ck*ψ(xk)

        αk=1
        counter = 0
        while ϕ(xk+αk*dk,ck)>ϕ(xk,ck)+σ*αk*Δk && counter<maxiter
            αk*=θ
            counter+=1
        end

        if counter == maxiter
             println("Não foi possível encontrar um comprimento de passo que decresce a função mérito na iteração $iter.")
                return xk, λk, μk, -2
        end
        
        # 3. Atualizar x e λ
        xk += αk*dk
        λk = λ_new
        μk = μ_new
    end

    println("Número máximo de iterações atingido.")
    return xk, λk, μk, -3
end
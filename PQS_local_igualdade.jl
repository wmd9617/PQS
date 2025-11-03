using JuMP, Clarabel, LinearAlgebra, Printf

include("MetodoGradiente.jl")

function PQS_local_equal(h, gradf, Jacobian_h, hessf, hess_coordenada_h, x0, λ0; maxiter=100, tol=1.e-6, fat=2)
    xk = copy(x0) # Usar copy para evitar modificar o x0 original
    λk = copy(λ0)
    l = length(λ0)
    n = length(x0)
    p=1.e-4
    
    # Gradiente da função Lagrangiana (para o critério de parada)
    gradL(x, λ) = [gradf(x) + Jacobian_h(x)'*λ; h(x)]
    
    # Hessiana da Lagrangiana aumentada em relação a x
    hess_L(x, λ; c=p)=hessf(x)+sum(λ[i]*hess_coordenada_h[i](x) for i in 1:l)+c*Jacobian_h(x)'*Jacobian_h(x)
    
    println("Iniciando PQS local para restrições de igualdade...")
    println("Iter | Norma do Grad. Lagrangiano")
    println("------------------------------------")

    for iter in 0:maxiter
        norm_gradL = norm(gradL(xk, λk))
        println(@sprintf("%4d | %.6e", iter, norm_gradL))

        if norm_gradL < tol
            println("Convergência alcançada!")
            return xk, λk, iter
        end

        # 1. Construir o modelo do subproblema quadrático no passo 'd'
        Hk = hess_L(xk, λk)

        if isposdef(Hk)==false
            for i in 1:maxiter 
                p=p*fat
                Hk=hess_L(xk, λk, c=p)
                if isposdef(Hk)==true
                    p=1.e-4
                    break
                end
            end
            if isposdef(Hk)==false
                println("Não foi possível encontrar uma matriz definida positiva na iteração $iter.")
                return xk, λk, -1
            end
        end

        gk = gradf(xk)
        Jk = Jacobian_h(xk)
        hk = h(xk)
        
        model = Model(Clarabel.Optimizer)
        set_silent(model)
        
        @variable(model, d[1:n])
        
        # Objetivo é uma aproximação quadrática da Lagrangiana
        @objective(model, Min, 0.5 * d' * Hk * d + gk' * d)
        
        # Restrição é a linearização da restrição original h(x) = 0
        @constraint(model, c1, Jk * d .== -hk)
        
        optimize!(model)
        
        if !is_solved_and_feasible(model)
            println("O subproblema não pôde ser resolvido na iteração $iter.")
            return xk, λk, -1
        end

        # 2. Obter o passo 'd' e os novos multiplicadores
        dk = value.(d)

        λ_new = -dual.(c1)
        
        # 3. Atualizar x e λ
        xk += dk
        λk = λ_new
    end

    println("Número máximo de iterações atingido.")
    return xk, λk, -2
end
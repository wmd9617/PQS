# Para os problemas de otimização abaixo temos que
# n - número de variáveis
# l - número de restrições de igualdade
# m - número de restrições de desigualdade

# Referência: Lecture Notes in Economics and Mathematical Systems by Wili Hock & Klaus Schittkowski

# Problemas com restrições de igualdade
#-----------------------------------------------------------------------
function PtReq1() # Problema 7 página 30 da referência (n=2 l=2)
	
	# função objetivo
    f(x)=log(1+x[1]^2)-x[2]
    # gradiente da função objetivo
    gradf(x)=[2*x[1]/(1+x[1]^2);-1]
    # hessiana da função objetivo
    hessf(x)=[2-2*x[1]^2 0; 0 0] 
    
    # restrição de igualdade
    h(x)=(1+x[1]^2)^2+x[2]^2-4
    # Jacobiano da restrição de igualdade
    Jacobian_h(x)=[4*x[1]+4*x[1]^3 2*x[2]] 
 	# hessiana da restrição de igualdade
    hessh1(x)=[4+12x[1]^2 0; 0 2]; hess_coordenada_h=[hessh1] 
    
    # chute inicial
    x0=[2.0;2.0]; λ0=rand()
    
    return f, h, gradf, Jacobian_h, hessf, hess_coordenada_h, x0, λ0
    
end

function PtReq2() # Problema 70 página 99 da referência (n=5 l=3)
    
    # função objetivo
    f(x)=(x[1]-1)^2+(x[1]-x[2])^2+(x[2]-x[3])^2+(x[3]-x[4])^4+(x[4]-x[5])^4
    # gradiente da função objetivo
    gradf(x)=[4*x[1]-2*x[2]-2;-2*x[1]+4*x[2]-2*x[3];-2*(x[2]-x[3])+4*(x[3]-x[4])^3;4*(x[4]-x[5])^3-4*(x[3]-x[4])^3;-4*(x[4]-x[5])^3]
    # hessiana da função objetivo
    hessf(x)=[4 -2.0 0 0 0; -2 4 -2 0 0; 0 -2 2+12*(x[3]-x[4])^2 -12*(x[3]-x[4])^2 0; 
    			0 0 -12*(x[3]-x[4])^2 12*(x[4]-x[5])^2+12*(x[3]-x[4])^2 -12*(x[4]-x[5])^2; 0 0 0 -12*(x[4]-x[5])^2 12*(x[4]-x[5])^2] 
    
    # restrições de igualdade
    h(x)=[x[1]+x[2]^2+x[3]^3-2-3*sqrt(2); x[2]-x[3]^2+x[4]+2-2*sqrt(2); x[1]*x[5]-2]
    # Jacobiano do campo vetorial das restrições de igualdade
    Jacobian_h(x)=[1 2*x[2] 3*x[3]^2 0 0; 0 1 -2*x[3] 1 0; x[5] 0 0 0 x[1]] 
 	# hessiana das restrições de igualdade
    hessh1(x)=[0 0 0 0 0; 0 2 0 0 0; 0 0 6*x[3] 0 0; 0 0 0 0 0; 0 0 0 0 0]
    hessh2(x)=[0 0 0 0 0; 0 0 0 0 0; 0 0 -2 0 0; 0 0 0 0 0; 0 0 0 0 0]
    hessh3(x)=[0 0 0 0 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0]
    hess_coordenada_h=[hessh1; hessh2; hessh3] 
    
    # chute inicial
    x0=[2.0;2.0;2.0;2.0;2.0]; λ0=rand(3)
    
    return f, h, gradf, Jacobian_h, hessf, hess_coordenada_h, x0, λ0
    
end

# Problemas com restrições de igualdade e desigualdade
#-------------------------------------------------------------------------
function PtRineq1() # Problema 14 página 37 da referência (n=2 l=1 m=1)
	
	# função objetivo
	f(x)=(x[1]-2)^2+(x[2]-1)^2
	# gradiente da função objetivo
	gradf(x)=[2*(x[1]-2); 2*(x[2]-1)]
	# hessiana da função objetivo
	hessf(x)=[2 0; 0 2]
	
	# restrição de igualdade
	h(x)=x[1]-2*x[2]+1
	# Jacobiano da restrição de igualdade
	Jacobian_h(x)=[1 -2]
	# hessiana da restrição de igualdade
	hessh1(x)=[0.0 0.0; 0.0 0.0]; hess_coordenada_h=[hessh1]
	
	# restrição de desigualdade
	g(x)=0.25*x[1]^2+x[2]^2-1
	# gradiente da restrição de desigualdade
	gradg1(x)=[0.5*x[1]; 2*x[2]]; grad_coordenada_g=[gradg1]
	# Jacobiano da restrição de desigualdade
   	Jacobian_g(x)=[0.5*x[1] 2*x[2]]
   	# hessiana da restrição de desigualdade
    hessg1(x)=[0.5 0; 0 2]; hess_coordenada_g=[hessg1]
	
    # chute inicial
	x0=[2;2]; λ0=rand(); μ0=rand()
	
	return f, g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0
	
end

function PtRineq2() # Problema 14 página 101 da referência (n=5 l=3 m=10)
	
	# função objetivo
	f(x)=exp(x[1]*x[2]*x[3]*x[4]*x[5])-0.5*(x[1]^3+x[2]^3+1)^2
	# gradiente da função objetivo
	gradf(x)=[x[2]*x[3]*x[4]*x[5]*exp(x[1]*x[2]*x[3]*x[4]*x[5])-(x[1]^3+x[2]^3+1)*3*x[1]^2;
				x[1]*x[3]*x[4]*x[5]*exp(x[1]*x[2]*x[3]*x[4]*x[5])-(x[1]^3+x[2]^3+1)*x[2]^2;
				x[1]*x[2]*x[4]*x[5]*exp(x[1]*x[2]*x[3]*x[4]*x[5]);
				x[1]*x[2]*x[3]*x[5]*exp(x[1]*x[2]*x[3]*x[4]*x[5]);
				x[1]*x[3]*x[3]*x[4]*exp(x[1]*x[2]*x[3]*x[4]*x[5])]
	# hessiana da função objetivo
	hessf(x)=[x[2]^2*x[3]^2*x[4]^2*x[5]^2*exp(x[1]*x[2]*x[3]*x[4]*x[5])-(15*x[1]^4+6*x[1]*x[2]^3+6*x[1]) x[3]*x[4]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5])-9*x[1]^2*x[2]^2 x[2]*x[4]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[2]*x[3]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[2]*x[3]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]);
	x[3]*x[4]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5])-9*x[1]^2*x[2]^2 x[1]^2*x[3]^2*x[4]^2*x[5]^2*exp(x[1]*x[2]*x[3]*x[4]*x[5])-(6*x[1]^3*x[2]+15*x[2]^4+6*x[2]) x[1]*x[4]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[3]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[3]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]);
	x[2]*x[4]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[4]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]^2*x[2]^2*x[4]^2*x[5]^2*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[2]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[2]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]); 
	x[2]*x[3]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[3]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[2]*x[5]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]^2*x[2]^2*x[3]^2*x[5]^2*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[2]*x[3]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]);
	x[2]*x[3]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[3]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[2]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]*x[2]*x[4]*(1+x[1]*x[2]*x[3]*x[4]*x[5])*exp(x[1]*x[2]*x[3]*x[4]*x[5]) x[1]^2*x[2]^2*x[3]^2*x[4]^2*exp(x[1]*x[2]*x[3]*x[4]*x[5])]
	
	# restrições de igualdade
	h(x)=[x[1]^2+x[2]^2+x[3]^2+x[4]^2+x[5]^2-10; x[2]*x[3]-5*x[4]*x[5]; x[1]^3+x[2]^3+1]
	# Jacobiano das restriçẽs de igualdade
	Jacobian_h(x)=[2*x[1] 2*x[2] 2*x[3] 2*x[4] 2*x[5]; 0 x[3] x[2] -5*x[5] -5*x[4]; 3*x[1]^2 3*x[2]^2 0 0 0]
	# hessiana das restrições de igualdade
	hessh1(x)=[2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0.; 0 0 0 0 2]
	hessh2(x)=[0 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 0 -5; 0 0 0 -5 0]
	hessh3(x)=[6*x[1] 0 0 0 0; 0 6*x[2] 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
	hess_coordenada_h=[hessh1; hessh2; hessh3]
	
	# restrições de desigualdade
	g(x)=[-2.3-x[1]; x[1]-2.3; -2.3-x[2]; x[2]-2.3; -3.2-x[3]; x[3]-3.2; -3.2-x[4]; x[4]-3.2; -3.2-x[5]; x[5]-3.2]
	# gradiente das restrições de desigualdade
	gradg1(x)=[-1;0;0;0;0]
	gradg2(x)=[1;0;0;0;0]
	gradg3(x)=[0;-1;0;0;0]
	gradg4(x)=[0;1;0;0;0]
	gradg5(x)=[0;0;-1;0;0]
	gradg6(x)=[0;0;1;0;0]
	gradg7(x)=[0;0;0;-1;0]
	gradg8(x)=[0;0;0;1;0]
	gradg9(x)=[0;0;0;0;-1]
	gradg10(x)=[0;0;0;0;1]
	grad_coordenada_g=[gradg1;gradg2;gradg3;gradg4;gradg5;gradg6;gradg7;gradg8;gradg9;gradg10]
	# Jacobiano das restrição de desigualdade
   	Jacobian_g(x)=[gradg1(x)';gradg2(x)';gradg3(x)';gradg4(x)';gradg5(x)';gradg6(x)';gradg7(x)';gradg8(x)';gradg9(x)';gradg10(x)']
   	# hessiana das restrições de desigualdade
    hessg1(x)=zeros(5,5);hessg2(x)=zeros(5,5);hessg3(x)=zeros(5,5);hessg4(x)=zeros(5,5);hessg5(x)=zeros(5,5)
    hessg6(x)=zeros(5,5);hessg7(x)=zeros(5,5);hessg8(x)=zeros(5,5);hessg9(x)=zeros(5,5);hessg10(x)=zeros(5,5)
    hess_coordenada_g=[hessg1;hessg2;hessg3;hessg4;hessg5;hessg6;hessg7;hessg8;hessg9;hessg10]
	
    # chute inicial
	x0=[-2;2;2;-1;-1]; λ0=rand(3); μ0=rand(10)
	
	return f, g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0
	
end

# Referência: Nonlinear Optimization Applications Using the GAMS Technology by Neculai Andrei

# Problema aplicado
function Pa()
	
	# função objetivo
	f(x)=3000*x[1]+1000*x[1]^3+2000*x[2]+666.667*x[2]^3
	# gradiente da função objetivo
	gradf(x)=[3000+3000*x[1]^2; 2000+3*666.667*x[2]^2;0;0;0;0;0;0;0]
	# hessiana da função objetivo
	v=zeros(9); C=sin(0.25)*48.4/50.176; D=cos(0.25)*48.4/50.176
	hessf(x)=[6000*x[1] 0 0 0 0 0 0 0 0; 0 6*666.667*x[2] 0 0 0 0 0 0 0;v';v';v';v';v';v';v']
	# restrições de igualdade
	h(x)=[0.4-x[1]+2*C*x[5]^2+x[5]*x[6]*(D*sin(-x[8])-C*cos(-x[8]))+x[5]*x[7]*(D*sin(-x[9])-C*cos(-x[9]));
		  0.4-x[2]+2*C*x[6]^2+x[5]*x[6]*(D*sin(x[8])-C*cos(x[8]))+x[6]*x[7]*(D*sin(x[8]-x[9])-C*cos(x[8]-x[9]));
		  0.8+2*C*x[7]^2+x[5]*x[7]*(D*sin(x[9])-C*cos(x[9]))+x[6]*x[7]*(D*sin(x[9]-x[8])-C*cos(x[9]-x[8]));
		  0.2-x[3]+2*D*x[5]^2+x[5]*x[6]*(C*sin(-x[8])+D*cos(-x[8]))-x[5]*x[7]*(C*sin(-x[9])+D*cos(-x[9]));
		  0.2-x[4]+2*D*x[6]^2-x[5]*x[6]*(C*sin(x[8])+D*cos(x[8]))-x[6]*x[7]*(C*sin(x[8]-x[9])+D*cos(x[8]-x[9]));
		  -0.337+2*D*x[7]^2-x[5]*x[7]*(C*sin(x[9])+D*cos(x[9]))-x[6]*x[7]*(C*sin(x[9]-x[8])+D*cos(x[9]-x[8]))]
	# Jacobiano das restriçẽs de igualdade
	Jacobian_h(x)=[-1 0 0 0 4*C*x[5]+x[6]*(D*sin(-x[8])-C*cos(-x[8]))+x[7]*(D*sin(-x[9])-C*cos(-x[9])) x[5]*(D*sin(-x[8])-C*cos(-x[8])) x[7]*(D*sin(-x[9])-C*cos(-x[9])) -x[5]*x[6]*(D*cos(-x[8])+C*sin(-x[8]));
	0 -1 0 0 ]
	# hessiana das restrições de igualdade
	hessh1(x)=[v';v';v';v'; 0 0 0 0 4*C D*sin(-x[8])-C*cos(-x[8]) D*sin(-x[9])-C*cos(-x[9]) -x[6]*(D*cos(-x[8])+C*sin(-x[8])) -x[7]*(D*cos(-x[9])+C*sin(-x[9]));
	0 0 0 0 D*sin(-x[8])-C*cos([-x[8]])]
	hessh2(x)=[]
	hessh3(x)=[6*x[1] 0 0 0 0; 0 6*x[2] 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
	hess_coordenada_h=[hessh1; hessh2; hessh3]
	
	# restrições de desigualdade
	g(x)=[-2.3-x[1]; x[1]-2.3; -2.3-x[2]; x[2]-2.3; -3.2-x[3]; x[3]-3.2; -3.2-x[4]; x[4]-3.2; -3.2-x[5]; x[5]-3.2]
	# gradiente das restrições de desigualdade
	gradg1(x)=[-1;0;0;0;0]
	gradg2(x)=[1;0;0;0;0]
	gradg3(x)=[0;-1;0;0;0]
	gradg4(x)=[0;1;0;0;0]
	gradg5(x)=[0;0;-1;0;0]
	gradg6(x)=[0;0;1;0;0]
	gradg7(x)=[0;0;0;-1;0]
	gradg8(x)=[0;0;0;1;0]
	gradg9(x)=[0;0;0;0;-1]
	gradg10(x)=[0;0;0;0;1]
	grad_coordenada_g=[gradg1;gradg2;gradg3;gradg4;gradg5;gradg6;gradg7;gradg8;gradg9;gradg10]
	# Jacobiano das restrição de desigualdade
   	Jacobian_g(x)=[gradg1(x)';gradg2(x)';gradg3(x)';gradg4(x)';gradg5(x)';gradg6(x)';gradg7(x)';gradg8(x)';gradg9(x)';gradg10(x)']
   	# hessiana das restrições de desigualdade
    hessg1(x)=zeros(5,5);hessg2(x)=zeros(5,5);hessg3(x)=zeros(5,5);hessg4(x)=zeros(5,5);hessg5(x)=zeros(5,5)
    hessg6(x)=zeros(5,5);hessg7(x)=zeros(5,5);hessg8(x)=zeros(5,5);hessg9(x)=zeros(5,5);hessg10(x)=zeros(5,5)
    hess_coordenada_g=[hessg1;hessg2;hessg3;hessg4;hessg5;hessg6;hessg7;hessg8;hessg9;hessg10]
	
    # chute inicial
	x0=[-2;2;2;-1;-1]; λ0=rand(3); μ0=rand(10)
	
	return f, g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0
	
end

function PaI() # Problema 14 páginas 212-216 da referência (n=27 l=19 m=54)
	
	# função objetivo
	f(x)=1300*(2000/((x[1]*x[2])/3+(x[1]+x[2])/6))^0.6+1300*(1000/(2*(x[3]*x[4])/3+(x[3]+x[4])/6))^0.6+1300*(1500/(4*(x[5]*x[6])/3+(x[5]+x[6])/6))^0.6
	# gradiente da função objetivo
	gradf(x)=[-0.6*1300*(6*2000)^0.6*(2*x[1]*x[2]+x[1]+x[2])^-1.6*(2*x[2]+1);
				-0.6*1300*(6*2000)^0.6*(2*x[1]*x[2]+x[1]+x[2])^-1.6*(2*x[1]+1);
				-0.6*1300*(6*1000)^0.6*(4*x[3]*x[4]+x[3]+x[4])^-1.6*(4*x[4]+1);
				-0.6*1300*(6*1000)^0.6*(4*x[3]*x[4]+x[3]+x[4])^-1.6*(4*x[3]+1);
				-0.6*1300*(6*1500)^0.6*(8*x[5]*x[6]+x[5]+x[6])^-1.6*(8*x[6]+1);
				-0.6*1300*(6*1500)^0.6*(8*x[5]*x[6]+x[5]+x[6])^-1.6*(8*x[5]+1);
				0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]
	# hessiana da função objetivo
	hessf(x)=[0.96*1300*(6*2000)^0.6*(2*x[2]+1)^2*(2*x[1]*x[2]+x[1]+x[2])^-2.6 0.96*1300*(6*2000)^0.6*(2*x[1]+1)*(2*x[2]+1)*(2*x[1]*x[2]+x[1]+x[2])^-2.6-1.2*1300*(6*2000)^0.6*(2*x[1]*x[2]+x[1]+x[2])^-1.6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0.96*1300*(6*2000)^0.6*(2*x[1]*x[2]+x[1]+x[2])^-2.6*(2*x[1]+1)*(2*x[2]+1)-1.2*1300*(6*2000)^0.6*(2*x[1]*x[2]+x[1]+x[2])^-1.6 0.96*1300*(6*2000)^0.6*(2*x[1]*x[2]+x[1]+x[2])^-2.6*(2*x[1]+1)^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0.96*1300*(6*1000)^0.6*(4*x[4]+1)^2*(4*x[3]*x[4]+x[3]+x[4])^-2.6 0.96*1300*(6*1000)^0.6*(4*x[3]+1)*(4*x[4]+1)*(4*x[3]*x[4]+x[3]+x[4])^-2.6-2.4*1300*(6*1000)^0.6*(4*x[3]*x[4]+x[3]+x[4])^-1.6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0.96*1300*(6*1000)^0.6*(4*x[3]*x[4]+x[3]+x[4])^-2.6*(4*x[3]+1)*(4*x[4]+1)-1.2*1300*(6*1000)^0.6*(4*x[3]*x[4]+x[3]+x[4])^-1.6 0.96*1300*(6*1000)^0.6*(4*x[3]*x[4]+x[3]+x[4])^-2.6*(4*x[3]+1)^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0.96*1300*(6*1500)^0.6*(8*x[6]+1)^2*(8*x[5]*x[6]+x[5]+x[6])^-2.6 0.96*1300*(6*1500)^0.6*(8*x[5]+1)*(8*x[6]+1)*(8*x[5]*x[6]+x[5]+x[6])^-2.6-4.8*1300*(6*1500)^0.6*(8*x[5]*x[6]+x[5]+x[6])^-1.6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0.96*1300*(6*1500)^0.6*(8*x[5]*x[6]+x[5]+x[6])^-2.6*(8*x[5]+1)*(8*x[6]+1)-4.8*1300*(6*1500)^0.6*(8*x[5]*x[6]+x[5]+x[6])^-1.6 0.96*1300*(6*1500)^0.6*(8*x[5]*x[6]+x[5]+x[6])^-2.6*(8*x[5]+1)^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
	
	# restrições de igualdade
	h(x)=[x[7]+x[12]+x[17]-15; 
	x[7]+x[14]+x[20]-x[8];
	x[12]+x[9]+x[19]-x[13];
	x[17]+x[10]+x[15]-x[18];
	x[11]+x[9]+x[10]-x[8];
	x[16]+x[14]+x[15]-x[13];
	x[21]+x[19]+x[20]-x[18];
	100*x[7]+x[25]*x[14]+x[27]*x[20]-x[22]*x[8];
	100*x[12]+x[23]*x[9]+x[27]*x[19]-x[24]*x[13];
	100*x[17]+x[23]*x[10]+x[25]*x[15]-x[26]*x[18];
	x[8]*(x[23]-x[22])-2000;
	x[13]*(x[25]-x[24])-1000;
	x[18]*(x[27]-x[26])-1500;
	x[1]+x[23]-210;
	x[2]+x[22]-130;
	x[3]+x[25]-210;
	x[4]+x[24]-160;
	x[5]+x[27]-210;
	x[6]+x[26]-180]
	# Jacobiano das restriçẽs de igualdade
#x=(∆T11,∆T12,∆T21,∆T22,∆T31,∆T32, f11, f12, f13,  f14,  f15,  f21,  f22,  f23,  f24,  f25,  f31,  f32,  f33,  f34,  f35,  T1i,T1o,T2i,T2o,T3i,T3o)
#x=(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25])
	Jacobian_h(x)=[
# ∆T11 ∆T12 ∆T21 ∆T22 ∆T31 ∆T32 f11 f12 f13 f14 f15 f21 f22 f23 f24 f25 f31 f32 f33 f34 f35 T1i T1o T2i T2o T3i T3o 
   0    0    0    0    0    0    1   0   0   0   0   1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0; # f11+f21+f31-45
   0    0    0    0    0    0    1  -1   0   0   0   0   0   1   0   0   0   0   0   1   0   0   0   0   0   0   0; # f11+f23+f34-f12
   0    0    0    0    0    0    0   0   1   0   0   1  -1   0   0   0   0   0   1   0   0   0   0   0   0   0   0; # f21+f13+f33-f22
   0    0    0    0    0    0    0   0   0   1   0   0   0   0   1   0   1  -1   0   0   0   0   0   0   0   0   0; # f31+f14+f24-f32
   0    0    0    0    0    0    0  -1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; # f15+f13+f14-f12
   0    0    0    0    0    0    0   0   0   0   0   0  -1   1   1   1   0   0   0   0   0   0   0   0   0   0   0; # f25+f23+f24-f22
   0    0    0    0    0    0    0   0   0   0   0   0   0   0   0   0   0  -1   1   1   1   0   0   0   0   0   0; # f35+f33+f34-f32
   0    0    0    0    0    0   100 x[22] 0  0   0   0   0 x[25] 0   0   0   0   0 x[27] 0 -x[8] 0  0 x[14] 0 x[20];#100f11+T2o*f23+T3o*f34-T1i*f12
   0    0    0    0    0    0    0   0 x[23] 0   0  100 -x[24] 0 0   0   0   0 x[27] 0  0   0 x[9] -x[13] 0 0 x[19];#100f21+T1o*f13+T3o*f33-T2i*f22
   0    0    0    0    0    0    0   0   0 x[23] 0   0   0   0 x[25] 0  100 -x[26] 0 0  0 0 x[10] 0 x[15] -x[18] 0;# 100f31+T1o*f14+T2o*f24-T3i*f32
   0    0    0    0    0    0    0 x[23]-x[22] 0 0 0 0   0   0   0   0   0   0   0   0   0 -x[8] x[8] 0   0  0   0; # f12(T1o-T1i)-2000
   0    0    0    0    0    0    0   0   0   0   0   0 x[25]-x[24] 0 0 0 0   0   0   0   0   0   0 -x[13] x[13] 0 0; # f22(T2o-T2i)-1000
   0    0    0    0    0    0    0   0   0   0   0   0   0   0   0   0   0 x[27]-x[26] 0 0 0 0   0   0   0 -x[18] x[18]; # f32(T3o-T3i)-1500
   1    0    0    0    0    0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   0   0   0   0; # ∆T11-T1o-210
   0    1    0    0    0    0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   0   0   0   0   0; # ∆T12-T1i-130
   0    0    1    0    0    0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   0   0; # ∆T21-T2o-210
   0    0    0    1    0    0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   0   0   0; # ∆T22-T2i-160
   0    0    0    0    1    0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1; # ∆T31-T3o-210
   0    0    0    0    0    1    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   0 # ∆T32-T3i-180
	]
	
	# hessiana das restrições de igualdade
	v=zeros(27)
	hessh1(x)=zeros(27,27);hessh2(x)=zeros(27,27);hessh3(x)=zeros(27,27);hessh4(x)=zeros(27,27)
	hessh5(x)=zeros(27,27);hessh6(x)=zeros(27,27);hessh7(x)=zeros(27,27)
	hessh8(x)=[v';v';v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; v';v';v';v';v'; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; v'; 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; v'; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
	hessh9(x)=[v';v';v';v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;v';v';v';0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;v';v';v';0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;v';v';0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
	hessh10(x)=[v';v';v';v';v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;v';v';v';v';0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;v';v';v';v';0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;v'
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;v']
	hessh11(x)=[v';v';v';v';v';v';v';0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0;v';v';v';v';v';v';v';v';v';v';v';v';v';
	0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;v';v';v';v']
	hessh12(x)=[v';v';v';v';v';v';v';v';v';v';v';v';0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0;v';v';v';v';v';v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;v';v']
	hessh13(x)=[v';v';v';v';v';v';v';v';v';v';v';v';v';v';v';v';v'; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1;v';v';v';v';v';v';v';
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
	hessh14(x)=zeros(27,27);hessh15(x)=zeros(27,27);hessh16(x)=zeros(27,27)
	hessh17(x)=zeros(27,27);hessh18(x)=zeros(27,27);hessh19(x)=zeros(27,27)
	hess_coordenada_h=[hessh1;hessh2;hessh3;hessh4;hessh5;hessh6;hessh7;hessh8;hessh9;hessh10;hessh11;hessh12;hessh13;hessh14;hessh15;hessh16;
						hessh17;hessh18;hessh19]

#x=(∆T11,∆T12,∆T21,∆T22,∆T31,∆T32, f11, f12, f13,  f14,  f15,  f21,  f22,  f23,  f24,  f25,  f31,  f32,  f33,  f34,  f35,  T1i,T1o,T2i,T2o,T3i,T3o)
#x=(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25])
	# restrições de desigualdade
	g(x)=[10-x[1]; x[1]-110;
	      10-x[2]; x[2]-110; 
	      10-x[3]; x[3]-110; 
	      10-x[4]; x[4]-110; 
	      10-x[5]; x[5]-110; 
	      10-x[6]; x[6]-110;
		  -x[7]; x[7]-45; 
		  -x[8]; x[8]-45; 
		  -x[9]; x[9]-45; 
		  -x[10]; x[10]-45; 
		  -x[11]; x[11]-45; 
		  -x[12]; x[12]-45; 
		  -x[13]; x[13]-45;
		  -x[14]; x[14]-45; 
		  -x[15]; x[15]-45; 
		  -x[16]; x[16]-45; 
		  -x[17]; x[17]-45; 
		  -x[18]; x[18]-45; 
		  -x[19]; x[19]-45; 
		  -x[20]; x[20]-45;
		  -x[21]; x[21]-45;
		  100-x[22]; x[22]-200;
		  100-x[23]; x[23]-200;
		  100-x[24]; x[24]-200;
		  100-x[25]; x[25]-200;
		  100-x[26]; x[26]-200;
		  100-x[27]; x[27]-200]
	# gradiente das restrições de desigualdade
	#          01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
	gradg1(x)=[-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg2(x)=[ 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg3(x)=[ 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg4(x)=[ 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg5(x)=[ 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg6(x)=[ 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg7(x)=[ 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg8(x)=[ 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
	gradg9(x)=[ 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg10(x)=[ 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg11(x)=[ 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg12(x)=[ 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg13(x)=[ 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg14(x)=[ 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg15(x)=[ 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg16(x)=[ 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg17(x)=[ 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg18(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg19(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg20(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg21(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg22(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg23(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg24(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg25(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg26(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg27(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg28(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg29(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg30(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg31(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg32(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg33(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg34(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg35(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg36(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg37(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg38(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0]
   gradg39(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0; 0]
   gradg40(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0]
   gradg41(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0; 0]
   gradg42(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0]
   gradg43(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0; 0]
   gradg44(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0]
   gradg45(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0; 0]
   gradg46(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0]
   gradg47(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0; 0]
   gradg48(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0]
   gradg49(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0; 0]
   gradg50(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0]
   gradg51(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1; 0]
   gradg52(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0]
   gradg53(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1]
   gradg54(x)=[ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1]
	
	grad_coordenada_g=[gradg1;gradg2;gradg3;gradg4;gradg5;gradg6;gradg7;gradg8;gradg9;gradg10;
					   gradg11;gradg12;gradg13;gradg14;gradg15;gradg16;gradg17;gradg18;gradg19;gradg20;
					   gradg21;gradg22;gradg23;gradg24;gradg25;gradg26;gradg27;gradg28;gradg29;gradg30;
					   gradg31;gradg32;gradg33;gradg34;gradg35;gradg36;gradg37;gradg38;gradg39;gradg40;
					   gradg41;gradg42;gradg43;gradg44;gradg45;gradg46;gradg47;gradg48;gradg49;gradg50;
					   gradg51;gradg52;gradg53;gradg54]
	# Jacobiano das restrição de desigualdade
   	Jacobian_g(x)=[gradg1(x)';gradg2(x)';gradg3(x)';gradg4(x)';gradg5(x)';gradg6(x)';gradg7(x)';gradg8(x)';gradg9(x)';gradg10(x)';
   				   gradg11(x)';gradg12(x)';gradg13(x)';gradg14(x)';gradg15(x)';gradg16(x)';gradg17(x)';gradg18(x)';gradg19(x)';gradg20(x)';
   				   gradg21(x)';gradg22(x)';gradg23(x)';gradg24(x)';gradg25(x)';gradg26(x)';gradg27(x)';gradg28(x)';gradg29(x)';gradg30(x)';
   				   gradg31(x)';gradg32(x)';gradg33(x)';gradg34(x)';gradg35(x)';gradg36(x)';gradg37(x)';gradg38(x)';gradg39(x)';gradg40(x)';
   				   gradg41(x)';gradg42(x)';gradg43(x)';gradg44(x)';gradg45(x)';gradg46(x)';gradg47(x)';gradg48(x)';gradg49(x)';gradg50(x)';
   				   gradg51(x)';gradg52(x)';gradg53(x)';gradg54(x)']
   	# hessiana das restrições de desigualdade
    hessg1(x)=zeros(27,27)
    hessg2(x)=zeros(27,27)
    hessg3(x)=zeros(27,27)
    hessg4(x)=zeros(27,27)
    hessg5(x)=zeros(27,27)
    hessg6(x)=zeros(27,27)
    hessg7(x)=zeros(27,27)
    hessg8(x)=zeros(27,27)
    hessg9(x)=zeros(27,27)
    hessg10(x)=zeros(27,27)
    hessg11(x)=zeros(27,27)
    hessg12(x)=zeros(27,27)
    hessg13(x)=zeros(27,27)
    hessg14(x)=zeros(27,27)
    hessg15(x)=zeros(27,27)
    hessg16(x)=zeros(27,27)
    hessg17(x)=zeros(27,27)
    hessg18(x)=zeros(27,27)
    hessg19(x)=zeros(27,27)
    hessg20(x)=zeros(27,27)
    hessg21(x)=zeros(27,27)
    hessg22(x)=zeros(27,27)
    hessg23(x)=zeros(27,27)
    hessg24(x)=zeros(27,27)
    hessg25(x)=zeros(27,27)
    hessg26(x)=zeros(27,27)
    hessg27(x)=zeros(27,27)
    hessg28(x)=zeros(27,27)
    hessg29(x)=zeros(27,27)
    hessg30(x)=zeros(27,27)
    hessg31(x)=zeros(27,27)
    hessg32(x)=zeros(27,27)
    hessg33(x)=zeros(27,27)
    hessg34(x)=zeros(27,27)
    hessg35(x)=zeros(27,27)
    hessg36(x)=zeros(27,27)
    hessg37(x)=zeros(27,27)
    hessg38(x)=zeros(27,27)
    hessg39(x)=zeros(27,27)
    hessg40(x)=zeros(27,27)
    hessg41(x)=zeros(27,27)
    hessg42(x)=zeros(27,27)
    hessg43(x)=zeros(27,27)
    hessg44(x)=zeros(27,27)
    hessg45(x)=zeros(27,27)
    hessg46(x)=zeros(27,27)
    hessg47(x)=zeros(27,27)
    hessg48(x)=zeros(27,27)
    hessg49(x)=zeros(27,27)
    hessg50(x)=zeros(27,27)
    hessg51(x)=zeros(27,27)
    hessg52(x)=zeros(27,27)
    hessg53(x)=zeros(27,27)
    hessg54(x)=zeros(27,27)
    hess_coordenada_g=[hessg1;hessg2;hessg3;hessg4;hessg5;hessg6;hessg7;hessg8;hessg9;hessg10;
    				   hessg11;hessg12;hessg13;hessg14;hessg15;hessg16;hessg17;hessg18;hessg19;hessg20;
    				   hessg21;hessg22;hessg23;hessg24;hessg25;hessg26;hessg27;hessg28;hessg29;hessg30;
    				   hessg31;hessg32;hessg33;hessg34;hessg35;hessg36;hessg37;hessg38;hessg39;hessg40;
    				   hessg41;hessg42;hessg43;hessg44;hessg45;hessg46;hessg47;hessg48;hessg49;hessg50;
    				   hessg51;hessg52;hessg53;hessg54]
	
    # chute inicial
	x0=100*ones(27); λ0=rand(19); μ0=rand(54)
	
	return f, g, h, gradf, grad_coordenada_g, Jacobian_h, Jacobian_g, hessf, hess_coordenada_h, hess_coordenada_g, x0, λ0, μ0
	
end
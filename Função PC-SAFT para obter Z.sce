//A função a seguir utiliza a PC-SAFT para calcular o fator de compressibilidade de uma mistura ou de um puro: todas as equações do código vêm da referência Gross e Sadowski 2001
//Os parâmetros de entrada são: fração de empacotamento das moléculas, T(K), vetor fração molar , vetor m , vetor sigma , vetor epsilon/K , vetor de interação binária

clear
clc

function Z = SAFTZ(eta,T,x,m,sigma,epsilon_K,k)
    //Constantes universais da Teoria da Perturbação
    a0 = [0.9105631445;0.6361281449;2.6861347891;-26.547362491;97.759208784;-159.59154087;91.297774084]
    a1 = [-0.3084016918;0.1860531159;-2.5030047259;21.419793629;-65.255885330;83.318680481;-33.746922930]
    a2 = [-0.0906148351;0.4527842806;0.5962700728;-1.7241829131;-4.1302112531;13.776631870;-8.6728470368]
    b0 = [0.7240946941;2.2382791861;-4.0025849485;-21.003576815;26.855641363;206.55133841;-355.60235612]
    b1 = [-0.5755498075;0.6995095521;3.8925673390;-17.215471648;192.67226447;-161.82646165;-165.20769346]
    b2 = [0.0976883116;-0.2557574982;-9.1558561530;20.642075974;-38.804430052;93.626774077;-29.66690558]
    A = [a0,a1,a2]
    B = [b0,b1,b2]
    U = [A,B]
    
    for i = 1:length(x)
        d(i) = sigma(i)*(1-0.12*exp(-3*epsilon_K(i)/T)) //Angstrons
    end
    
    soma = 0
    for i = 1:length(x)
        soma = soma + x(i)*m(i)*d(i)^3
    end
    ro = 6/%pi*eta*soma^-1//número de densidade de moléculas (Angstrons^-3)

    for n = 1:4
        soma2 = 0
        for i = 1:length(x)
            soma2 = soma2 + x(i)*m(i)*d(i)^(n-1)
        end
        zeta(n) = %pi/6*ro*soma2
    end
    zeta0 = zeta(1)
    zeta1 = zeta(2)
    zeta2 = zeta(3)
    zeta3 = zeta(4)
    
    Zhs = zeta3/(1-zeta3) + 3*zeta1*zeta2/(zeta0*(1-zeta3)^2) + (3*zeta2^3-zeta3*zeta2^3)/(zeta0*(1-zeta3)^3)
    
    for i = 1:length(x)
        for j = 1:length(x)
            ghs(i,j) = 1/(1-zeta3) + d(i)*d(j)/(d(i)+d(j))*3*zeta2/(1-zeta3)^2 + (d(i)*d(j)/(d(i)+d(j)))^2*2*zeta2^2/(1-zeta3)^3
            rodghs_dro(i,j) = zeta3/(1-zeta3)^2 + d(i)*d(j)/(d(i)+d(j))*(3*zeta2/(1-zeta3)^2+6*zeta2*zeta3/(1-zeta3)^3) + (d(i)*d(j)/(d(i)+d(j)))^2*(4*zeta2^2/(1-zeta3)^3+6*zeta2^2*zeta3/(1-zeta3)^4)
        end
    end
    
    mmedio = 0
    for i = 1:length(x)
        mmedio = mmedio + x(i)*m(i)
    end
    
    soma3 = 0
    for i = 1:length(x)
        soma3 = soma3 + x(i)*(m(i)-1)*ghs(i,i)^-1*rodghs_dro(i,i)
    end
    Zhc = mmedio*Zhs - soma3

    for i = 1:7
        a(i) = A(i,1) + (mmedio-1)/mmedio*A(i,2) + (mmedio-1)/mmedio*(mmedio-2)/mmedio*A(i,3)
        b(i) = B(i,1) + (mmedio-1)/mmedio*B(i,2) + (mmedio-1)/mmedio*(mmedio-2)/mmedio*B(i,3)
    end
    
    I1 = 0
    I2 = 0
    for i = 1:7
        I1 = I1 + a(i)*eta^(i-1)
        I2 = I2 + b(i)*eta^(i-1)
    end
    
    detaI1_deta = 0
    detaI2_deta = 0
    for i = 1:7
        detaI1_deta = detaI1_deta + a(i)*i*eta^(i-1)
        detaI2_deta = detaI2_deta + b(i)*i*eta^(i-1)
    end
    
    C1 = (1 + mmedio*(8*eta-2*eta^2)/(1-eta)^4 + (1-mmedio)*(20*eta-27*eta^2+12*eta^3-2*eta^4)/((1-eta)*(2-eta))^2)^-1
    C2 = -C1^2*(mmedio*(-4*eta^2+20*eta+8)/(1-eta)^5 + (1-mmedio)*(2*eta^3+12*eta^2-48*eta+40)/((1-eta)*(2-eta))^3)

    for i = 1:length(x)
        for j = 1:length(x)
            sigmaij(i,j) = (sigma(i)+sigma(j))/2
            epsilonij_K(i,j) = (1-k(i,j))*sqrt(epsilon_K(i)*epsilon_K(j))
        end
    end

    m2es3 = 0
    m2e2s3 = 0
    for i = 1:length(x)
        for j = 1:length(x)
            m2es3 = m2es3 + x(i)*x(j)*m(i)*m(j)*epsilonij_K(i,j)/T*sigmaij(i,j)^3
            m2e2s3 = m2e2s3 + x(i)*x(j)*m(i)*m(j)*(epsilonij_K(i,j)/T)^2*sigmaij(i,j)^3
        end
    end

    Zdisp = -2*%pi*ro*detaI1_deta*m2es3 - %pi*ro*mmedio*(C1*detaI2_deta+C2*eta*I2)*m2e2s3
    
    Z = 1 + Zhc + Zdisp
endfunction








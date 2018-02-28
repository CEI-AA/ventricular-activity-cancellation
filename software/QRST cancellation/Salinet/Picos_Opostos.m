function[teste,novopico] = Picos_Opostos(sinal,vetor,pico)

global Fam

teste = 0;
dmin = 0.280*Fam; %% Distancia minima entre picos: 280 ms %%


k = length(vetor);

if k > 1 | vetor(1) ~= 0
    temp1 = sinal(vetor(k)+1)-sinal(vetor(k));
    if sinal(pico)*sinal(vetor(k)) < 0
        teste = 1;
        for i = vetor(k) + 1 : pico-3
            temp = sinal(i) - sinal(i-1);
        
            if temp*temp1 < 0
                teste = 0;
                i=pico-3;
            end
        end
        
    else teste = 0;
    end

end


if teste == 1 | pico - vetor(k) < dmin
    
    if abs(sinal(pico)) > abs(sinal(vetor(k)))
        novopico = pico;
    else novopico = vetor(k);
    end
else novopico = pico; 
end


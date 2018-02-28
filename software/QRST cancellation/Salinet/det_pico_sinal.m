%%% Funçao que realiza a detecçao de pico em sinal filtrado %%%

function[amostra]=det_pico_sinal(inicio,sinal,final,limiar,dmin)

i = inicio;
condicao = 0;
subida = 0;
descida = 0;
amostra = [];
cont_dmin = 0;


while i < final
    
    
    subida = acima_limiar(i,sinal,limiar);
    if subida == 0
        break;
    end
    
    descida = abaixo_limiar(subida,sinal,limiar);
    if descida == 0
        break;
    end
    
    i = descida + 1;
        
    
    if subida ~= 0 & descida ~= 0
    
        %keyboard;
        pospico = pico_limiar(subida,descida,sinal);
        
        if length(amostra) > 0 & pospico < final
            if abs(pospico - amostra(end)) < dmin
                if sinal(pospico) > sinal(amostra(end))
                    amostra(end) = pospico;
                    subida = 0;
                    descida = 0;
                    cont_dmin = 0;
                end
            else amostra = [amostra pospico];
                 subida = 0;
                 descida = 0;
                 cont_dmin = 0;
            end
        else
            if pospico < final
                amostra = [amostra pospico];
                subida = 0;
                descida = 0;
                cont_dmin = 0;
            end
        end
        
        
        
    end
    
    
    
end



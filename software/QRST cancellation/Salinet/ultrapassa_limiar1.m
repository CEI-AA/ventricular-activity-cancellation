function[pos] =ultrapassa_limiar1(inicio,sinal,final,limiar,vtpic,dmin)
 

i = inicio;
condicao = 0;

pos = 0;
if length(vtpic) == 0

    bool=0;
    while bool == 0 & i < final
        if sinal(i) > limiar
            bool=1;
        else i=i+1;
        end
    end
    pos=i;

else 

    bool=0;
    while bool == 0 & i < final
        if (sinal(i) > limiar & (i-vtpic(end)) > dmin) | (sinal(i) > limiar & sinal(i) > sinal(vtpic(end)))
            bool=1;
        else
            i=i+1;
        end
    end
    pos=i;

        
end
 
    


function[pos] =ultrapassa_limiar2(inicio,sinal,final,limiar)
 

i = inicio;
bool = 0;

pos = 0;
while bool == 0 & i < final
    if sinal(i) < limiar
            bool=1;
    else i=i+1;
    end
end

if bool == 1
    pos=i;    
end
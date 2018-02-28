function[pos] = abaixo_limiar(inicio,ident,limiarTemp)
global ECG limiar m0
i=inicio;
bool=0;
pos=0;
if length(ident) == 1
    while bool == 0 & i < length(ECG)
        if abs(ECG(i)-m0) < limiar
            bool=1;
        else i=i+1;
        end
    end
    pos=i;
else while bool == 0 & i < length(ident)
        if abs(ident(i)) < limiarTemp
            bool=1;
        else i=i+1;
        end
    end
    if bool == 1
        pos=i;
    end
end
        
function[amostra] = pico_limiar(inicio,fim,ident)

global ECG limiar m0

i=inicio;
amostra=inicio;

if length(ident) == 1
    while i < fim
        if abs(ECG(i)-m0) > abs(ECG(amostra)-m0)
            amostra = i;
        end
        i=i+1;
    end
else while i < fim
        if abs(ident(i)) > abs(ident(amostra))
            amostra = i;
        end
        i=i+1;
    end
end
    
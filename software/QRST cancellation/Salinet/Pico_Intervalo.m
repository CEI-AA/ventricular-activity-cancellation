function[amostra] = Pico_Intervalo(inicio,fim,sinal)

i=inicio;
amostra=inicio;
while i < fim
    if abs(sinal(i)) > abs(sinal(amostra))
        amostra = i;
    end
    i=i+1;
end

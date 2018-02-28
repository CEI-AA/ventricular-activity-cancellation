function[sinalF]=CalcMM(sinal,tjan)

MM = tjan;
at = (MM-1)/2;
sinalF = zeros(1,length(sinal));
for i=1:length(sinal)
    inicio = max(1,i-at);
    fim = min(length(sinal),i+at);
    sinalF(i) = mean(sinal(inicio:fim));
end



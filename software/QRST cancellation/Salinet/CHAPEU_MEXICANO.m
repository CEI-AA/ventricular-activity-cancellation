%%% Implementaçao da DWT, Escala = 1,2,4,8,16... %%%

function [Sinal, X] = CHAPEU_MEXICANO(Escala)
passo = 1/Escala;
X = -5:passo:5;
Sinal = 2.1741*(1/sqrt(2*pi) .* (1 - X.^2) .* exp(-X.^2/(2)));
%keyboard;
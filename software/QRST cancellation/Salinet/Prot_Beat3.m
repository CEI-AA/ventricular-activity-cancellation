function[] = Prot_Beat3

global escala Fam sfreq anot ECG IBB limiar picos_R m0 On Off
global mediaIntervalo desvioIntervalo tipo
global Se Pr Der FN FP TP QRSver pproc
global nrotFP nrotFPc1 nrotFPc2 nrotFPc3
global MetricaTeste m1FPc1 m2FPc1 EscalaSel
global Nbatmin Nbatmax FPlimiarEscala exame
global  P POn POff T TOn TOff Final_yw_intra RRmin VInd

RR = picos_R(2:end)-picos_R(1:end-1);
RRmin = min(RR);
VInd = round(-0.30*RRmin):round(0.70*RRmin);

pattern = [];

ct = 1;

Mbeats = zeros(length(picos_R)-2,length(VInd));
Mbeats(ct,:) = Final_yw_intra(picos_R(2)-round(0.30*RRmin):picos_R(2)+round(0.70*RRmin));



for ibt = 3:length(picos_R)-1
    
    ct = ct+1;
    inicio = picos_R(ibt)-round(0.30*RRmin);
    fim = picos_R(ibt)+round(0.70*RRmin);
    janela = Final_yw_intra(inicio:fim);
    Mbeats(ct,:) = aligncc3(Mbeats(1,:),janela);
        
end


ct=0;

for ibt = 1:length(VInd)
    
   ct=ct+1;
   pattern(ct) = median(Mbeats(:,ibt));
      
end
%pattern = pattern/ct;

pattern = pattern';

save pattern;

function[] = Prot_Beat2

global escala Fam sfreq anot ECG IBB limiar picos_R m0 On Off
global mediaIntervalo desvioIntervalo tipo
global Se Pr Der FN FP TP QRSver pproc
global nrotFP nrotFPc1 nrotFPc2 nrotFPc3
global MetricaTeste m1FPc1 m2FPc1 EscalaSel
global Nbatmin Nbatmax FPlimiarEscala exame
global  P POn POff T TOn TOff Final_yw_intra RRmin

RR = picos_R(2:end)-picos_R(1:end-1);
RRmin = min(RR);
janela = Final_yw_intra(picos_R(2)-round(0.30*RRmin):picos_R(2)+round(0.70*RRmin));
pattern = [];
ct = 0;

for ibt = 2:length(picos_R)-1
    
    
    
   if length(pattern) == 0
       pattern = janela;
       ct = ct+1;
   else inicio = picos_R(ibt)-round(0.30*RRmin);
        fim = picos_R(ibt)+round(0.70*RRmin);
        janela = Final_yw_intra(inicio:fim);
        pico = picos_R(ibt)-inicio+1;
        pattern = aligncc(pattern,janela);
        ct = ct+1;
        
   end
   
end
pattern = pattern/ct;


save pattern;

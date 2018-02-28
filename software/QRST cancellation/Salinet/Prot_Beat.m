function[] = Prot_Beat

global escala Fam sfreq anot ECG IBB limiar picos_R m0 On Off
global mediaIntervalo desvioIntervalo tipo
global Se Pr Der FN FP TP QRSver pproc
global nrotFP nrotFPc1 nrotFPc2 nrotFPc3
global MetricaTeste m1FPc1 m2FPc1 EscalaSel
global Nbatmin Nbatmax FPlimiarEscala exame
global  P POn POff T TOn TOff Final_yw_intra Tendv


% DORmin = min(picos_R(2:end)-On);
% RTmin = min(T+0.150*sfreq-picos_R(2:end-1));
% VInd = -DORmin:RTmin;
% 
% pattern = [];
% 
% ct = 0;
% 
% for ibt = 1:length(VInd)
%     
%    ct=ct+1;
%    pattern(ct) = median(Final_yw_intra(picos_R(2:end-1)+VInd(ibt)));
%       
% end

SegmsOnT = Tendv-On(1:end);
Lpat = round(mean(SegmsOnT));
SetCycles = zeros(length(SegmsOnT),Lpat);
pattern = zeros(1,Lpat);



for iw = 1:length(SegmsOnT)
    
    if iw == 1
        
            
        if SegmsOnT(iw) < Lpat
            SetCycles(iw,1:SegmsOnT(iw)+1) = Final_yw_intra(On(iw):Tendv(iw));
        else SetCycles(iw,:) = Final_yw_intra(On(iw):On(iw)+Lpat-1);
        end
        
    else ref=SetCycles(1,:);
         s2 = Final_yw_intra(On(iw):Tendv(iw));
         s2m = aligncc3(ref,s2);
         if length(s2m) < Lpat
             SetCycles(iw,1:SegmsOnT(iw)+1) = s2m;
         else SetCycles(iw,:) = s2m(1:Lpat);
         end
    end
        
end

for im = 1:Lpat
    pattern(im)=median(SetCycles(:,im));
end
    
    

% janela = Final_yw_intra(On(1):Tendv(1));
% 
% pattern = [];
% ct = 0;
% 
% for ibt = 2:length(T)
%     
%     
%     
%    if length(pattern) == 0
%        pattern = janela;
%        ct = ct+1;
%    else inicio = On(ibt);
%         fim = Tendv(ibt);
%         janela = Final_yw_intra(inicio:fim);
%         pico = picos_R(ibt+1)-inicio+1;
%         %pico = find(abs(janela)==max(abs(janela)));
%         pattern = somawin(pattern,janela,pico);
%         ct = ct+1;
%         
%    end
%    
% end
% pattern = pattern/ct;

save pattern;

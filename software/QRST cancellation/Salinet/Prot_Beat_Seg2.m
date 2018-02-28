function[pattern,segment,Onr,Offr,begin_seg,end_seg] = Prot_Beat_Seg2(begin_seg,end_seg,segment,n)

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


R_seg = find(picos_R >= begin_seg & picos_R <= end_seg);

RR = (picos_R(R_seg(2):R_seg(end)) - picos_R(R_seg(1):R_seg(end-1)));
RRmin = min(RR);
OnABS = picos_R(R_seg) - round(0.30*RRmin);
OffABS = picos_R(R_seg) + round(0.70*RRmin);


On_seg = find(OnABS >= begin_seg & OnABS <= end_seg);
Tend_seg = find(OffABS >= begin_seg & OffABS <= end_seg);

if OnABS(1) < begin_seg
    begin_seg = OnABS(1);
    On_seg = find(OnABS >= begin_seg & OnABS <= end_seg);
end

if OffABS(end) > end_seg
    end_seg = OffABS(end);
    Tend_seg = find(OffABS >= begin_seg & OffABS <= end_seg);
end

if On_seg(end) > Tend_seg(end)
    end_seg = On(On_seg(end))-1;
    On_seg = find(On >= begin_seg & On <= end_seg);
end

segment=Final_yw_intra(begin_seg:end_seg,n);

% Ons = On(On_seg); % Tends = Tendv(Tend_seg);

Lpat = RRmin+1;

%SegmsOnT = Tendv(Tend_seg)-On(On_seg);
%Lpat = round(mean(SegmsOnT));
SetCycles = zeros(length(R_seg),Lpat);
pattern = zeros(1,Lpat);

Onr = OnABS(On_seg)-begin_seg+1;
Offr = OffABS(Tend_seg) - begin_seg+1;



for iw = 1:length(R_seg)
    
    if iw == 1
        
            
        SetCycles(iw,:) = segment(Onr(iw):Offr(iw));
        
        
    else ref=SetCycles(1,:);
         s2 = segment(Onr(iw):Offr(iw));
         s2m = aligncc3(ref,s2);
         SetCycles(iw,:) = s2m(1:Lpat);
         
    end
        
end

for im = 1:Lpat
    pattern(im)=mean(SetCycles(:,im));
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

%save pattern;

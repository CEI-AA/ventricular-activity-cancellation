%function[picos_R,EscalaSel] = Det_QRS_H(exame)

global escala Fam sfreq anot ECG IBB limiar picos_R m0 On Off
global mediaIntervalo desvioIntervalo tipo
global Se Pr Der FN FP TP QRSver pproc
global nrotFP nrotFPc1 nrotFPc2 nrotFPc3
global MetricaTeste m1FPc1 m2FPc1 EscalaSel
global Nbatmin Nbatmax FPlimiarEscala exame
global  P POn POff T TOn TOff Final_yw_intra RRmin VInd Tendv Iter

%%% Aquisicao do sinal ECG %%%


Modelo_Subtractor = 1;
%% 1 - QRS onset and T end
%% 2 - Fixed window average
%% 3 - Fixed window median




%% first column - ECG II
%% second column - ECG V1
%% eleventh column - ECG I
%% twelfth column - ECG V5


%% External ECG Data Input
exame = ECG;
%ECG = exame;

%% Intracardiac Noncontact Signals Input
%load Pat8_base1_v.mat;
Intracardiac_signals_raw= Pat2_Intracardaic_signals_raw;
%%

index = anot;
sfreq = 1200;  % Fs=2034.5 Hz (Pat 8) and Fs=1200 Hz (Remaining Patients)
Fam = sfreq;
Nsoma = 0;
Nprod = 0;
Nsub = 0;
Ndiv = 0;
Nrq = 0;
Iter = 0;
vbol = 1;

%%

Pattern_Matrix=[];
Pattern_Delay=[];
IntracardiacSignal_AfterSubtraction = [];
RRmin = min(picos_R(2:end)-picos_R(1:end-1));

for n=1:length (Intracardiac_signals_raw(1,:)) 
    
    close all;
    Final_yw_intra= Intracardiac_signals_raw(:,n);

    if Modelo_Subtractor == 1
        Prot_Beat;
    else if Modelo_Subtractor == 2
            Prot_Beat2;
        else Prot_Beat3;
        end
    end

% figure;

% plot((1:length(exame))/sfreq,exame,'r');hold on;plot(picos_R/sfreq,exame(picos_R),'ro');plot(On/sfreq,exame(On),'ok');plot(Off/sfreq,exame(Off),'ok');
% plot((1:length(Final_yw_intra))/sfreq,Final_yw_intra);plot(picos_R/sfreq,Final_yw_intra(picos_R),'ro');plot(On/sfreq,Final_yw_intra(On),'ok');plot(Off/sfreq,Final_yw_intra(Off),'ok');
% figure;plot((1:length(exame))/sfreq,exame,'r');hold on;plot(picos_R/sfreq,exame(picos_R),'ro');plot(On/sfreq,exame(On),'ok');plot(Off/sfreq,exame(Off),'ok');
% plot(Tendv/sfreq,exame(Tendv),'o');



load pattern;

if length(Pattern_Matrix) == 0
    Pattern_Matrix = zeros(length(pattern),length(Intracardiac_signals_raw(1,:)));
end

if length(Pattern_Delay) == 0
    Pattern_Delay = zeros(length(picos_R)-1,length(Intracardiac_signals_raw(1,:)));
end

Pattern_Matrix(:,n) = pattern;

lp = length(pattern);
if Modelo_Subtractor == 1
    pattern = pattern';
end
Final_yw_intra2 = Final_yw_intra;


%figure;plot(pattern);



for ii=1:length(picos_R)-1
    
    if Modelo_Subtractor == 1
        win2pro = Final_yw_intra2(On(ii):Tendv(ii));
        delaypeak = aligncc2(pattern,win2pro);
        Pattern_Delay(ii,n) = delaypeak;
        lseg=length(win2pro);
        
        if lseg <= lp
            Final_yw_intra2(On(ii)-delaypeak:Tendv(ii)-delaypeak) = Final_yw_intra2(On(ii)-delaypeak:Tendv(ii)-delaypeak) - pattern(1:lseg);
        else Final_yw_intra2(On(ii)-delaypeak:On(ii)-delaypeak+lp-1) = Final_yw_intra2(On(ii)-delaypeak:On(ii)-delaypeak+lp-1) - pattern;
        end
    else win2pro = Final_yw_intra2(picos_R(ii)-0.30*RRmin:picos_R(ii)+0.70*RRmin);
         delaypeak = aligncc2(pattern,win2pro);
         Pattern_Delay(ii,n) = delaypeak;
         if delaypeak >= picos_R(ii)-round(0.30*RRmin)
             delaypeak = 0;
         end
         Final_yw_intra2(picos_R(ii)-round(0.30*RRmin)-delaypeak:picos_R(ii)+round(0.70*RRmin)-delaypeak) = Final_yw_intra2(picos_R(ii)-round(0.30*RRmin)-delaypeak:picos_R(ii)+round(0.70*RRmin)-delaypeak) - pattern;
    end
        
    
end
IntracardiacSignal_AfterSubtraction(:,n) = Final_yw_intra2;


end
keyboard;

figure;subplot(3,1,1);
plot((1:length(exame))/sfreq,exame,'r');hold on;plot(picos_R/sfreq,exame(picos_R),'ro');plot(On/sfreq,exame(On),'ok');plot(Off/sfreq,exame(Off),'ok');
plot(Tendv/sfreq,exame(Tendv),'ok');

subplot(3,1,2);
plot((1:length(Final_yw_intra))/sfreq,Final_yw_intra);


RR = picos_R(2:end)-picos_R(1:end-1);
RRmin = min(RR);
VInd = round(-0.30*RRmin):round(0.70*RRmin);
subplot(3,1,3);
plot((1:length(Final_yw_intra2))/sfreq,Final_yw_intra2);
hold on;
for ii=1:length(picos_R)-1
    if Modelo_Subtractor == 1
        plot((On(ii):On(ii)+lp-1)/sfreq,Final_yw_intra2(On(ii):On(ii)+lp-1),'k');
    else plot((picos_R(ii)+VInd)/sfreq,Final_yw_intra2(picos_R(ii)+VInd),'k');
    end
end
hold on;plot(picos_R/sfreq,Final_yw_intra2(picos_R),'ro');



figure;

Final_yw_intra_corr = Final_yw_intra(picos_R(1)+300:picos_R(end)-200);
Final_yw_intra2_corr = Final_yw_intra2(picos_R(1)+300:picos_R(end)-200);


% Fsam = 1200;
% NN = 18000;
% sizefft = NN;
% fstep = Fsam /sizefft; 
% limit = 30;
% Sfft1 = abs(fft(Final_yw_intra_corr(1:NN)));
% Sfft2 = abs(fft(Final_yw_intra2_corr(1:NN)));
% 
% 
% Sfft1 = Sfft1.*Sfft1/length (Sfft1);
% Sfft2 = Sfft2.*Sfft2/length (Sfft2);
% 
% % subplot (2,1,1); plot ((1:length(Final_yw_intra(1:limit/fstep)))*(fstep), Sfft1(1:limit/fstep)); 
% % subplot (2,1,2); plot ((1:length(Final_yw_intra2(1:limit/fstep)))*(fstep), Sfft2(1:limit/fstep));
% 
% plot ((1:length(Final_yw_intra_corr(1:limit/fstep)))*(fstep), Sfft1(1:limit/fstep)); 
% hold on;
% plot ((1:length(Final_yw_intra2_corr(1:limit/fstep)))*(fstep), Sfft2(1:limit/fstep), 'r');

%end

%plot(PMap/sfreq,Final_yw_intra(PMap),'og')
            
            
        
        
        
    


%keyboard;






%plot(histstd);

%Ext_Param;


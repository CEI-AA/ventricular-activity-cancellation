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
exame = Pat7_ECG;
ECG = exame;

%% Intracardiac Noncontact Signals Input
%load Pat8_base1_v.mat;
Intracardiac_signals_raw= Pat7_p1_Intracard;
%%

index = anot;
sfreq = 1200;
Fam = sfreq;
Nsoma = 0;
Nprod = 0;
Nsub = 0;
Ndiv = 0;
Nrq = 0;
Iter = 0;
vbol = 0;

%%% Inicializacao de Variaveis %%%

while vbol == 0

picos_R = 0;
IBB = 0;
On = 0;
Off = 0;

fator1 = 1.5;
fator2 = -1.5;
bproc = 0;
nproc = 0;
nrotFP = 0;

nrotFPc1 = 0;
nrotFPc2 = 0;
nrotFPc3 = 0;


%%% Aprendizado do sinal - analise dos primeiros batimentos %%%

T0 = 10*Fam;   %%% Tomando-se 10 segundos de intervalo inicial
Int0 = ECG(1:T0);

m0 = mean(Int0);
Int0 = Int0 - m0;

SetEscala = [8 16 32 64];
desvio_pontos_fid = [];
desvio_picos_R1 = [];
OcorrRuido = [];
Nbatmin = 3;
Nbatmax = (2*T0/Fam);
vetFPLimiar = [];

for valorEscala = 1:length(SetEscala)
    

    picos_R1=[];
    [Wavelet_Filha,XW] = CHAPEU_MEXICANO(SetEscala(valorEscala));
    FInt0 = conv(Wavelet_Filha,Int0);
    nproc = nproc + length(Int0);
    atraso = round(length(Wavelet_Filha)/2);
    FInt0w = FInt0(atraso:length(FInt0)-atraso);
    FInt0 = FInt0w;
    dFInt0 = diff(FInt0);
    FInt0 = dFInt0;
       
    HFInt0 = imag(Hilbert(FInt0)); %%% Transformada de Hilbert da derivada do sinal filtrado
    %EnvInt0 = sqrt(HFInt0.^2+FInt0.^2); %%% Envelope do sinal analitico correspondente 
    EnvInt0 = HFInt0.^2+FInt0.^2;
    

    [picos_R1,desvio_pfid,FPlimiarEscala] = Det_PrBat(EnvInt0,Int0);
    IBB = picos_R1(2:end)-picos_R1(1:end-1);
    desvio_picos_R1(valorEscala) = std(IBB);
    desvio_pontos_fid(valorEscala) = desvio_pfid;
    desvio_pontos_fid(valorEscala) = std(IBB)*desvio_pfid;
    vetFPLimiar(valorEscala) = FPlimiarEscala;
    OcorrRuido(valorEscala) = 0;
    
end

%%% Analise dos resultados de desvio-padrao das amplitudes dos pontos
%%% fiduciais

EscalaSel = [];
valorRef1 = find(OcorrRuido == 0);
valorRef = desvio_pontos_fid(valorRef1(1));
EscalaSel = valorRef1(1);
for valorEscala = valorRef1(1)+1:length(SetEscala)
    if desvio_pontos_fid(valorEscala) < valorRef & OcorrRuido(valorEscala) == 0
        valorRef = desvio_pontos_fid(valorEscala);
        EscalaSel = valorEscala;
    end
end

FPlimiarEscala = vetFPLimiar(EscalaSel);
EscalaSel = SetEscala(EscalaSel);

Wavelet_Filha = CHAPEU_MEXICANO(EscalaSel);


FInt0 = conv(Wavelet_Filha,Int0);
atraso = round(length(Wavelet_Filha)/2);
FInt0 = FInt0(atraso:length(FInt0)-atraso);
FInt0 = FInt0-mean(FInt0);
  %%% Derivada do sinal filtrado
FInt0 = diff(FInt0);

HFInt0 = imag(Hilbert(FInt0));
%%% Envelope do sinal analitico correspondente 
EnvInt0 = HFInt0.^2+FInt0.^2;
FPlimiarEscala
[picos_R,desvio_pfid] = Det_PrBat_ilust(EnvInt0,Int0);
IBB = picos_R(2:end)-picos_R(1:end-1);
desvio_picos_R(valorEscala) = std(IBB);

NBatTreino = 10;
[EscalaSeg] = Treino_OnOff_v5(length(picos_R)-1);


BatIni = length(picos_R);


%%% Levantamento de Informacoes Iniciais %%%

ModAmp = abs(ECG(picos_R)-m0);
PicoEst = mean(ModAmp);
limiar = 0.7*( PicoEst + abs(ECG(picos_R(BatIni))-m0) )/2;
for i = 2 : BatIni
    IBB(i-1) = picos_R(i) - picos_R(i-1);
end
mediaIntervalo = mean(IBB);
Nsoma = Nsoma + length(picos_R)-1;
Ndiv = Ndiv+1;
desvioIntervalo = std(IBB);
Nsub = Nsub+length(picos_R);
Nprod = Nprod+2*length(picos_R);
Ndiv = Ndiv+1;
Nrq = Nrq+1;
    
%%% Desenvolvimento das detecçoes %%%

inicio = abaixo_limiar(picos_R(BatIni),1,1);

%%% Laço Principal %%%

Sample = inicio;
Bat = length(picos_R);
Clim = 1;
dmin = 0.100*Fam;
limiar;
ultpico = 0;
nproc = 0;
histstd = 0;
chist = 0;
LimTsemQRS = 5*Fam;
FiltNsuc = 0;

while Sample < length(ECG)
    
    jan_m0 = ECG(picos_R(end)-2*Fam:min(picos_R(end)+2*Fam,length(ECG)));
    m0 = mean(jan_m0);
    HistLim(Clim) = limiar;
    Clim = Clim+1;
    inicio = acima_limiar(Sample,1,1);
    fim = abaixo_limiar(inicio,1,1);
    pico = pico_limiar(inicio,fim,1);
    IntervaloTeste = pico - picos_R(length(picos_R));

    if desvioIntervalo ~= 0
        
        MetricaTeste = (IntervaloTeste - mediaIntervalo)/desvioIntervalo;
             
        if MetricaTeste > fator1 & pico ~= ultpico

            k1 = length(picos_R);
            [fim,bproc]=RotinaFalsoNegativo_H(pico);
            nproc = nproc+bproc;
            if fim == 0
                ultpico=pico;
            end
                
        
            if k1 < length(picos_R)
                FiltNsuc = 0;
                atualiza_limiar;
                
                Nsoma = Nsoma+2;
                Nprod = Nprod+3;
                Ndiv = Ndiv+2;
                
                fim = picos_R(length(picos_R)) + 10;
                for idbfn = k1+1:length(picos_R)
                    Result = SegmQRS_v4(idbfn-1,EscalaSeg);
                end
            else FiltNsuc = 1;
            end
            
        
        else if MetricaTeste < fator2
            
                k1 = length(picos_R);
                bproc = RotinaFalsoPositivo_H(pico);
                nproc = nproc+bproc;
                if k1 < length(picos_R)
                    FiltNsuc = 0;
                    atualiza_limiar;
                    
                    Nsoma = Nsoma+2;
                    Nprod = Nprod+3;
                    Ndiv = Ndiv+2;
                    
                    
                    Result = SegmQRS_v4(length(picos_R)-1,EscalaSeg);
                else FiltNsuc = 1;
                   
                end
            else    if pico - picos_R(length(picos_R)) > dmin
                        FiltNsuc = 0;
                        picos_R(length(picos_R) + 1) = pico;
                        Result = SegmQRS_v4(length(picos_R)-1,EscalaSeg);
                        atualiza_limiar;
                        
                        Nsoma = Nsoma+2;
                        Nprod = Nprod+3;
                        Ndiv = Ndiv+2;
                        
                        
                        
                    else fator1 = fator1*0.7;
                         fator2 = fator2*0.7;
                         
                    end
            end
        end
        
    else if pico - picos_R(length(picos_R)) > dmin
            picos_R(length(picos_R) + 1) = pico;
            Result = SegmQRS_v4(length(picos_R)-1,EscalaSeg);
            atualiza_limiar;
            
            Nsoma = Nsoma+2;
            Nprod = Nprod+3;
            Ndiv = Ndiv+2;
            
            
        end
    end
    
    
   
    if fim ~= 0
        Sample = abaixo_limiar(fim,1,1);
    else Sample = abaixo_limiar(picos_R(length(picos_R)),1,1);
    end
        
    if Sample - picos_R(end) > LimTsemQRS & FiltNsuc == 0;
        Sample = round(picos_R(end) + dmin);
    end
    
    chist = chist+1;
    histstd(chist) = std(IBB)/Fam;
    

     
    for i = 2 : length(picos_R)
        IBB(i-1) = picos_R(i) - picos_R(i-1);
    end
    mediaIntervalo = mean(IBB);
    Nsoma = Nsoma+1;
    Ndiv = Ndiv+1;
    desvioIntervalo = std(IBB);
    Nsub = Nsub+length(IBB);
    Nprod = Nprod+2*length(IBB);
    Ndiv = Ndiv+1;
    Nrq = Nrq+1;
    Nsub=Nsub+1;
    Ndiv=Ndiv+1;
    temp = 0;
    
end

%picos_R = picos_R(1:end-1);
%On=On(1:end-1);
%Off=Off(1:end-1);

X = [1 On length(ECG)];
points = ECG(On);
points=points';

Y = [ECG(1) points ECG(end)];

%% Interpolation - Wandering analysis %%
%figure;plot((1:length(exame))/sfreq,exame,'r');hold on;plot(picos_R/sfreq,exame(picos_R),'ro');plot(On/sfreq,exame(On),'ok');plot(Off/sfreq,exame(Off),'ok');
YI(:,Iter+1) = INTERP1(X,Y,(1:length(ECG)));
InterSig = YI(:,Iter+1);
ECGC(:,Iter+1) = ECG-InterSig;
%figure;plot((1:length(ECGC))/sfreq,ECGC,'r');hold on;plot(picos_R/sfreq,ECGC(picos_R),'ro');plot(On/sfreq,ECGC(On),'ok');plot(Off/sfreq,ECGC(Off),'ok');
ECG = ECGC(:,Iter+1);

if Iter > 0
    VError (:,Iter) = YI(:,Iter+1) - YI(:,Iter);
    Eval (Iter) = sum(abs(VError (:,Iter)));
    if Eval (Iter) <= 0.01*(sum(abs(YI(:,Iter+1))))
        vbol = 1;
    end
    
end
    
    
Iter = Iter+1;

end

%%
Det_Peak_T2;



if length(Off) < length(picos_R)
    ECG = ECG(1:round(picos_R(end)-0.090*sfreq));
    exame = ECG;
    picos_R = picos_R(1:end-1);
    if length(Off) < length(On)
        On = On(1:end-1);
    end
    
end

RR = picos_R(2:end) - picos_R(1:end-1);

plot((1:length(exame))/sfreq,exame,'r');hold on;plot(picos_R/sfreq,exame(picos_R),'ro');plot(On/sfreq,exame(On),'ok');plot(Off/sfreq,exame(Off),'ok');
figure;plot((1:length(exame))/sfreq,exame,'r');hold on;plot(picos_R/sfreq,exame(picos_R),'ro');plot(On/sfreq,exame(On),'ok');plot(Off/sfreq,exame(Off),'ok');
plot(Tendv/sfreq,exame(Tendv),'o');plot(T/sfreq,exame(T),'go');

%%
keyboard;
Pattern_Matrix=[];
Pattern_Delay=[];
IntracardiacSignal_AfterSubtraction = [];

for n=1:length (Intracardiac_signals_raw(1,:)) 
    
    close all;
    Final_yw_intra= Intracardiac_signals_raw(:,n);
    Final_yw_intra = Final_yw_intra(1:length(ECG));

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



for ii=1:length(picos_R)
    
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
         ii
         n
         Final_yw_intra2(picos_R(ii)-round(0.30*RRmin)-delaypeak:picos_R(ii)+round(0.70*RRmin)-delaypeak) = Final_yw_intra2(picos_R(ii)-round(0.30*RRmin)-delaypeak:picos_R(ii)+round(0.70*RRmin)-delaypeak) - pattern;
    end
        
    
end
IntracardiacSignal_AfterSubtraction(:,n) = Final_yw_intra2;

n
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
for ii=1:length(picos_R)
    if Modelo_Subtractor == 1
        plot((On(ii):On(ii)+lp-1)/sfreq,Final_yw_intra2(On(ii):On(ii)+lp-1),'k');
    else plot((picos_R(ii)+VInd)/sfreq,Final_yw_intra2(picos_R(ii)+VInd),'k');
    end
end
hold on;plot(picos_R/sfreq,Final_yw_intra2(picos_R),'ro');



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


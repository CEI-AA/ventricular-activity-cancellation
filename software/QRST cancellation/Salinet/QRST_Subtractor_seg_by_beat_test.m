%function[picos_R,EscalaSel] = Det_QRS_H(exame)
global escala Fam sfreq anot ECG IBB limiar picos_R m0 On Off
global mediaIntervalo desvioIntervalo tipo
global Se Pr Der FN FP TP QRSver pproc
global nrotFP nrotFPc1 nrotFPc2 nrotFPc3
global MetricaTeste m1FPc1 m2FPc1 EscalaSel
global Nbatmin Nbatmax FPlimiarEscala exame
global  P POn POff T TOn TOff Final_yw_intra RRmin VInd Tendv Iter
%%% Aquisicao do sinal ECG %%%

tic
Modelo_Subtractor = 1;
%% 1 - QRS onset and T end
%% 2 - Fixed window average
%% 3 - Fixed window median


%% External ECG Data Input
exame = ECGs(:,1);
ECG = exame;

%% Intracardiac Noncontact Signals Input
%load Pat8_base1_v.mat;
Intracardiac_signals_raw= AEG;
time_window = 7; 
%%

index = anot;
% sfreq = 1200; %2034.5; %
sfreq = 2034.5; %
Fam = sfreq;
Nsoma = 0;
Nprod = 0;
Nsub = 0;
Ndiv = 0;
Nrq = 0;
Iter = 0;
vbol = 0;

%%% Inicializacao de Variaveis %%%
tic
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

T0 = 9*Fam;   %%% Tomando-se 10 segundos de intervalo inicial
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
       
    HFInt0 = imag(hilbert(FInt0)); %%% Transformada de Hilbert da derivada do sinal filtrado
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

HFInt0 = imag(hilbert(FInt0));
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

end

toc

%%
keyboard;


if length(Off) < length(picos_R)
    ECG = ECG(1:round(picos_R(end)-0.090*sfreq));
    exame = ECG;
    picos_R = picos_R(1:end-1);
    if length(Off) < length(On)
        On = On(1:end-1);
    end
    
end


% Det_Peak_T2;

keyboard;

RR = picos_R(2:end) - picos_R(1:end-1);

plot((1:length(ECG))/sfreq,ECG,'r');hold on;plot(picos_R/sfreq,ECG(picos_R),'ro');plot(On/sfreq,ECG(On),'ok');plot(Off/sfreq,ECG(Off),'ok');
figure;plot((1:length(ECG))/sfreq,ECG,'r');hold on;plot(picos_R/sfreq,ECG(picos_R),'ro');plot(On/sfreq,ECG(On),'ok');plot(Off/sfreq,ECG(Off),'ob');
plot(Tendv/sfreq,ECG(Tendv),'go');%plot(T/sfreq,exame(T),'go');
keyboard;
%%
sfreq = round (sfreq);
keyboard;
Pattern_Matrix=[];
Pattern_Delay=[];
IntracardiacSignal_AfterSubtraction = [];

% if Modelo_Subtractor == 1
%         Prot_Beat_Seg;
%     else if Modelo_Subtractor == 2
%             Prot_Beat2;
%         else Prot_Beat3;
%         end
% end
    
win_track_pattern=[sfreq*time_window]; 


for n=1:length (Intracardiac_signals_raw(1,:)) 
    Point=n
    
    close all;
    Final_yw_intra= Intracardiac_signals_raw(:,n);
    Final_yw_intra = Final_yw_intra(1:length(ECG));
    
    numberR = length(find(picos_R > win_track_pattern)); 
    refR = find(picos_R > win_track_pattern);
    
    nwins_pattern=floor(length(Final_yw_intra)/win_track_pattern);
    newbegin = 0; newend = 0;
    Final_yw_intra_sub_final = [];
    
     for n_patterns=1:numberR+1
       
       if n_patterns == 1
            signal_seg=Final_yw_intra(((n_patterns-1)*win_track_pattern+1):n_patterns*win_track_pattern);
            begin_seg (n_patterns) = ((n_patterns-1)*win_track_pattern+1);
            end_seg (n_patterns) = n_patterns*win_track_pattern;
       else signal_seg = Final_yw_intra((Tendv(refR(n_patterns-1))-win_track_pattern:Tendv(refR(n_patterns-1))));
            begin_seg (n_patterns) = (Tendv(refR(n_patterns-1))-win_track_pattern);
            end_seg (n_patterns) = Tendv(refR(n_patterns-1))+1;
       end
       refsignal = signal_seg;
       %begin_seg (n_patterns) = ((n_patterns-1)*win_track_pattern+1);
       %end_seg (n_patterns) = n_patterns*win_track_pattern;
       %if newend > 0 begin_seg (n_patterns) = newend + 1; end
       if Modelo_Subtractor == 1
        [pattern,newseg,Onr,Tendvr,newbegin,newend] = Prot_Beat_Seg_by_beat(begin_seg(n_patterns),end_seg(n_patterns),signal_seg,n_patterns);
       else [pattern,newseg,Onr,Tendvr,newbegin,newend] = Prot_Beat_Seg2_by_beat(begin_seg(n_patterns),end_seg(n_patterns),signal_seg);
       end
        signal_seg = newseg;
        patterns{n_patterns} = pattern;
        if n_patterns == 1 
            for n_QRST = 1:length(Onr)
                win2pro = signal_seg(Onr(n_QRST):Tendvr(n_QRST));
                delaypeak = aligncc2(pattern,win2pro);
                if delaypeak >= Onr(n_QRST)
                    delaypeak = Onr(n_QRST)-1;
                end
                lseg=length(win2pro);
                lp = length(pattern); [nl,nc] = size(pattern); if nc > 1 pattern = pattern'; end;
                %if lseg <= lp
                if n_QRST == length(Onr)
                   signal_seg(max(Onr(n_QRST)-delaypeak,1):min(Onr(n_QRST)-delaypeak+lp-1,length(signal_seg))) = signal_seg(max(Onr(n_QRST)-delaypeak,1):min(Onr(n_QRST)-delaypeak+lp-1,length(signal_seg))) - pattern(1:min(length(signal_seg),Onr(n_QRST)-delaypeak+lp-1)-Onr(n_QRST)+delaypeak+1);
                else signal_seg(max(Onr(n_QRST)-delaypeak,1):Onr(n_QRST)-delaypeak+lp-1) = signal_seg(max(Onr(n_QRST)-delaypeak,1):Onr(n_QRST)-delaypeak+lp-1) - pattern;
                     %else Final_yw_intra2(On(ii)-delaypeak:On(ii)-delaypeak+lp-1) = Final_yw_intra2(On(ii)-delaypeak:On(ii)-delaypeak+lp-1) - pattern;
                end
            end
            %figure;subplot(2,1,1);plot((newbegin:newend)/sfreq,newseg);
            %subplot(2,1,2);plot((newbegin:newend)/sfreq,signal_seg); 
            
        else 
               for n_QRST = 1:length(Onr)
                win2pro = signal_seg(Onr(n_QRST):Tendvr(n_QRST));
                delaypeak = aligncc2(pattern,win2pro);
                if delaypeak >= Onr(n_QRST)
                    delaypeak = Onr(n_QRST)-1;
                end
                lseg=length(win2pro);
                lp = length(pattern); [nl,nc] = size(pattern); if nc > 1 pattern = pattern'; end;
                %if lseg <= lp
                if n_QRST == length(Onr)
                   signal_seg(max(Onr(n_QRST)-delaypeak,1):min(Onr(n_QRST)-delaypeak+lp-1,length(signal_seg))) = signal_seg(max(Onr(n_QRST)-delaypeak,1):min(Onr(n_QRST)-delaypeak+lp-1,length(signal_seg))) - pattern(1:min(length(signal_seg),Onr(n_QRST)-delaypeak+lp-1)-Onr(n_QRST)+delaypeak+1);
                else signal_seg(max(Onr(n_QRST)-delaypeak,1):Onr(n_QRST)-delaypeak+lp-1) = signal_seg(max(Onr(n_QRST)-delaypeak,1):Onr(n_QRST)-delaypeak+lp-1) - pattern;
                     %else Final_yw_intra2(On(ii)-delaypeak:On(ii)-delaypeak+lp-1) = Final_yw_intra2(On(ii)-delaypeak:On(ii)-delaypeak+lp-1) - pattern;
                end
               end
               oldend = oldend - newbegin+1;
               %figure;subplot(2,1,1);plot((newbegin:newend)/sfreq,newseg);
               %subplot(2,1,2);plot((newbegin:newend)/sfreq,signal_seg); 
               signal_seg = signal_seg(oldend+1:end); 
               
               
         end
            
           
        Final_yw_intra2_cell{n_patterns}{n} =signal_seg;
        %figure;subplot(2,1,1);plot((1:length(refsignal))/sfreq,refsignal);
        
        oldend = newend;
        %keyboard;
     end
     
     


%% Concatenating segments

Final_yw_intra_sub_final = [];
 for icat = 1:n_patterns
     a = Final_yw_intra2_cell{1,icat}{n};
     Final_yw_intra_sub_final = cat(1,Final_yw_intra_sub_final,a);
 end

Final_yw_intra_sub_ready (:,n) = Final_yw_intra_sub_final;

end

toc;

save('Pat6_QRS_Tsub_1min.mat', 'Final_yw_intra_sub_ready')

aaa=1

save Pat6_QRS_Tsub_allvariabl_1min.mat


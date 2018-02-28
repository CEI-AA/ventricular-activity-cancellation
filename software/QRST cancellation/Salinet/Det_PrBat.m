function[vetorPFid,Desvio_ampl_pfid,fplim] = Det_PrBat(sinal,sinal_ref)

global limiar Fam
global Nbatmin Nbatmax Iter



MaxSinal = max(sinal);
sinal = sinal/MaxSinal;  %%% Normalizações
vbolParou = 0;
dlim = round(2*Fam); %%% distancia limite de 2s
pinic = round(0.1*Fam); %%% Ponto de partida a 100ms do inicio do sinal
pfin = length(sinal) - pinic;
fplim = 0.2;
fplimI = 0.2;
limiarAmp = fplim*max(sinal(pinic:pfin));
dmin = 0.120*Fam; %%% Distância mínima entre dois pontos fiduciais


while vbolParou == 0

limiarAmp = fplim*max(sinal(pinic:pfin));    
count = 1;
i = pinic;
mudlim = 0;
vetorPFid = [];
vetorPFid_Sfilt = [];

while i < pfin
    
    
    PosAcima = ultrapassa_limiar1(i,sinal,pfin,limiarAmp,vetorPFid_Sfilt,dmin);
    
    if PosAcima == 0
        break;
    end
    
    PosAbaixo = ultrapassa_limiar2(PosAcima,sinal,pfin,limiarAmp);
    
    if PosAbaixo == 0
        break;
    end
    
    PontoFiducial = Pico_Intervalo(PosAcima,PosAbaixo,sinal);
    
    if length(vetorPFid_Sfilt) > 0
        if abs(PontoFiducial - vetorPFid_Sfilt(end)) > dmin
            vetorPFid_Sfilt = [vetorPFid_Sfilt PontoFiducial];
            PicoMapeado = Mapeia_Pico(PontoFiducial,sinal_ref);
            PontoFiducial = PicoMapeado;
            vetorPFid = [vetorPFid PontoFiducial];
        else if sinal(PontoFiducial) > sinal(vetorPFid_Sfilt(end))
                vetorPFid_Sfilt(end) = PontoFiducial;
                PicoMapeado = Mapeia_Pico(PontoFiducial,sinal_ref);
                PontoFiducial = PicoMapeado;
                vetorPFid(end) = PontoFiducial;
            end
        end
    else vetorPFid_Sfilt = [vetorPFid_Sfilt PontoFiducial];
         PicoMapeado = Mapeia_Pico(PontoFiducial,sinal_ref);
         PontoFiducial = PicoMapeado;
         vetorPFid = [vetorPFid PontoFiducial];
    end
    
%     PicoMapeado = Mapeia_Pico(PontoFiducial,sinal_ref);
%     PontoFiducial = PicoMapeado;
%     vetorPFid = [vetorPFid PontoFiducial];    
        
     i = PosAbaixo+1;
     
end
 
Desvio_ampl_pfid = std(sinal(vetorPFid_Sfilt)); 
if (length(vetorPFid) >= Nbatmin & length(vetorPFid) < Nbatmax) 
%if length(vetorPFid) < Nbatmax
    vbolParou = 1;
else %fplim = fplim*1.2;
    if length(vetorPFid) < Nbatmin
        fplim = fplim*0.8;
    else fplim = fplim*1.2;
   end
end
end        


% limiarAmp = fplim*max(sinal(pinic:pfin));
%   figure;
%   subplot(2,1,1);
%   plot((1:length(sinal))/Fam,sinal);
%   hold on;
%   plot((1:length(sinal))/Fam,limiarAmp*ones(1,length(sinal)),'k');
%   hold on;
%   plot(vetorPFid_Sfilt/Fam,sinal(vetorPFid_Sfilt),'ro');
%   xlabel('tempo (segundos)');
%   ylabel('Unidades arbitrárias');
%   
%   subplot(2,1,2);
%   plot((1:length(sinal_ref))/Fam,sinal_ref);
%   hold on;
%   plot(vetorPFid/Fam,sinal_ref(vetorPFid),'ro');
%   xlabel('tempo (segundos)');
%   ylabel('Unidades arbitrárias');
% 
% pause; 
 


 
 
            
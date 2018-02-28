function[bproc] = RotinaFalsoPositivo_H(IndexPico)

global picos_R ECG IBB escala limiar
global Fam m0 nrotFP nrotFPc1 nrotFPc2 nrotFPc3
global mediaIntervalo desvioIntervalo MetricaTeste
global m1FPc1 m2FPc1
global EscalaSel FPlimiarEscala

nrotFP = nrotFP + 1;
Rot = 1;
EscTemp = escala;
VetEsc = [1];
t = 1;
IndexBat = picos_R(length(picos_R)); 
dmin = round(0.120*Fam);
bproc = 0;
nfilt = 0;
%keyboard;
if IndexPico - IndexBat < dmin
   m1 = ECG(IndexBat) - m0;
   m2 = ECG(IndexPico) - m0;
   if abs(m2) > abs(m1)
      picos_R(length(picos_R)) = IndexPico;
   end
   
else 
        
        offset = round((picos_R(end) - picos_R(end-1))/3);
        fim_intervalo = min(IndexPico + offset,length(ECG));
        inicio_intervalo = IndexBat - offset; 
        Intervalo = ECG(inicio_intervalo : fim_intervalo) - m0;
        SinalModif = Intervalo - mean(Intervalo);
        Wavelet_Filha = CHAPEU_MEXICANO(EscalaSel);
        SinalModif1 = conv(Wavelet_Filha,SinalModif);
        bproc = bproc + length(Intervalo);
        atraso = round(length(Wavelet_Filha)/2);
        SinalModif1 = SinalModif1(atraso:length(SinalModif1)-atraso);
        
        FInt0 = diff(SinalModif1);  %%% Derivada do sinal filtrado
        HFInt0 = imag(hilbert(FInt0));
        %EnvInt0 = sqrt(HFInt0.^2+FInt0.^2); %%% Envelope do sinal analitico correspondente 
        EnvInt0 = HFInt0.^2+FInt0.^2;
        SinalModif2 = EnvInt0;
        
                
        
        pinic = dmin;
        pfin = length(SinalModif2)-pinic;

	

    		
    		

%%% Calculo de um limiar temporario - a partir da 1ª escala %%%

    	MaxSinal = max(SinalModif2);
        SinalModif2 = SinalModif2/MaxSinal;  %%% Normalizações
	
        limiarTemp = FPlimiarEscala*max(SinalModif2(pinic:pfin));
        
    		

		i = 1;
    	maxi = length(SinalModif2);
    	bool = 0;
    	pico = 0;
    	dmin = ceil(0.120*Fam);
    	PicoIdent = 0;

        
       
        PicoIdent = det_pico_sinal(pinic,SinalModif2,pfin,limiarTemp,dmin);
        PontoFiducial=[];
        QminBat = floor(length(Intervalo)/Fam);
        QmaxBat = QminBat*2;
        
        %%% Aqui, mapear os dois picos a serem validados %%%
        
        %%% Construir função para, dentre os pontos fiduciais identificados, mapear os picos em teste %%%
        
        
        for iPM = 1:length(PicoIdent)
                PicoMapeado = Mapeia_Pico(PicoIdent(iPM),Intervalo);
                if length(PontoFiducial) > 0 
                    if PicoMapeado - PontoFiducial(end) > dmin
                        PontoFiducial = [PontoFiducial PicoMapeado];
                    end
                else
                    
                    PontoFiducial = [PontoFiducial PicoMapeado];
                    
                end
        end
         
        
        
        
    	if length(PontoFiducial) >= 2
            
            
               vbolPF1 = 0;
               PF1 = 0;
               vbolPF2 = 0;
               PF2 = 0;
               iPF = 1;
               while iPF <= length(PontoFiducial)
                   if abs(PontoFiducial(iPF) + inicio_intervalo - 1 - IndexBat) < dmin & PF1 == 0
                       vbolPF1 = 1;
                       PF1 = iPF;
                   end
                   if abs(PontoFiducial(iPF) + inicio_intervalo - 1 - IndexPico) < dmin & PF2 == 0
                       vbolPF2 = 1;
                       PF2 = iPF;
                   end
                   iPF = iPF + 1;
               end
               
               if vbolPF1 == 1 & vbolPF2 == 1
                   PontoFiducial = PontoFiducial([PF1 PF2]);
                   picos_R(length(picos_R)+1) = PontoFiducial(2) + inicio_intervalo - 1;
               end
                   
                   
%               figure(3);
%               subplot(2,1,1);
%               plot(SinalModif2);hold on;plot(PicoIdent,SinalModif2(PicoIdent),'ro');plot((1:length(SinalModif2)),limiarTemp*ones(1,length(SinalModif2)),'k');
%               subplot(2,1,2);
%               plot(Intervalo);hold on;plot(PontoFiducial,Intervalo(PontoFiducial),'ro');
%               keyboard;
%               close(figure(3));
            
             
   
%         	if (PontoFiducial(2) + inicio_intervalo - 1 - IndexBat) > dmin & (PontoFiducial(1) + inicio_intervalo - 1 - IndexBat) < dmin & (PontoFiducial(2) + inicio_intervalo - 1 - IndexPico) < dmin
%                                                
%                 picos_R(length(picos_R)+1) = PontoFiducial(2) + inicio_intervalo - 1;
%                 nrotFPc2 = nrotFPc2+1;
%                 
%         	else if abs(Intervalo(PontoFiducial(2))-m0) > abs(Intervalo(PontoFiducial(1))-m0)
%                 	picos_R(length(picos_R)) = PontoFiducial(2) + inicio_intervalo - 1;
%           
%             	 else picos_R(length(picos_R)) = PontoFiducial(1) + inicio_intervalo - 1;
%            
%                 end
% 
%             end
%      

           
        else

            if length(PontoFiducial) > 0
                
%              figure(3);
%              subplot(2,1,2);
%              plot(SinalModif2);hold on;plot(PicoIdent,SinalModif2(PicoIdent),'ro');plot((1:length(SinalModif2)),limiarTemp*ones(1,length(SinalModif2)),'k');
%              subplot(2,1,1);
%              plot(Intervalo);hold on;plot(PontoFiducial,Intervalo(PontoFiducial),'ro');
%              keyboard;
%              close(figure(3));

                
                    
                    if PontoFiducial ~= 0 & PontoFiducial + inicio_intervalo - 1 > IndexBat 
                        picos_R(length(picos_R)) = PontoFiducial + inicio_intervalo - 1;
      
      
                    end
                    
            end
   
   
      
   
       end
            
        



    		
    
    


end
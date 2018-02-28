function[fim,bproc] = RotinaFalsoNegativo_H(IndexPico)

global ECG escala picos_R IBB Fam m0
global mediaIntervalo desvioIntervalo limiar
global EscalaSel FPlimiarEscala

Rot = 2;
i = 1;
Tini = length(picos_R);
IndexBat = picos_R(length(picos_R));
OrdBat = length(picos_R);
Int1 = IndexBat - picos_R(OrdBat-1);
Int2 = IndexPico - IndexBat;
offset = round((mean(IBB))/2);
limiarTemp = 0;
bproc = 0;
jan_m0 = ECG(IndexBat + offset:IndexPico - offset);mo_a = m0;
m0 = mean(jan_m0);
mlim = 0.7*max(abs(ECG(IndexBat + offset:IndexPico - offset)-m0));

if Int2 > 3*mean(IBB) & mlim >=0.3*limiar & mlim <= 1.7*limiar
    fim = 0; limiarAnterior = limiar;
    limiar = 0.7*max(abs(ECG(IndexBat + offset:IndexPico - offset)-m0));
    
    %%% Alterando limiar e filtrando %%%
    
    offset1 = round((picos_R(end) - picos_R(end-1))/2);
    offset1 = round(0.120*Fam);
    offset2 = round((mean(IBB))/3);
    
    inicio = IndexBat + offset1;
    fim = IndexPico-offset1;
    Intervalo = ECG(inicio : fim) - m0;
    if length(Intervalo) > 0.090*Fam
        Wavelet_Filha = CHAPEU_MEXICANO(EscalaSel);
        SinalModif1 = conv(Wavelet_Filha,Intervalo);
        bproc = bproc + length(Intervalo);
        atraso = round(length(Wavelet_Filha)/2);
        SinalModif1 = SinalModif1(atraso:length(SinalModif1)-atraso);
        
        %FInt0 = DERIVADOR_SINAL(SinalModif1,2);  %%% Derivada do sinal filtrado
        FInt0 = diff(SinalModif1);  %%% Derivada do sinal filtrado
        
        %bproc = bproc + length(Intervalo);
        HFInt0 = imag(Hilbert(FInt0));
        %bproc = bproc + length(Intervalo);
        %EnvInt0 = sqrt(HFInt0.^2+FInt0.^2); %%% Envelope do sinal analitico correspondente 
        EnvInt0 = HFInt0.^2+FInt0.^2;
        
        SinalModif2 = EnvInt0;
        MaxSinal1 = max(SinalModif2);
        %SinalModif2 = SinalModif2/MaxSinal;  %%% Normalizações

        IndexBatAnt = picos_R(length(picos_R)-1);
        
        BatRecon = ECG(IndexBatAnt - offset2 : IndexBatAnt + offset2)-m0;
        BatReconMod1 = conv(Wavelet_Filha,BatRecon);
        bproc = bproc + length(BatRecon);
        BatReconMod1 = BatReconMod1(atraso:length(BatReconMod1)-atraso);
        %FInt0 = DERIVADOR_SINAL(BatReconMod1,2);  %%% Derivada do sinal filtrado
        FInt0 = diff(BatReconMod1);  %%% Derivada do sinal filtrado
        %bproc = bproc + length(BatRecon);
        HFInt0 = imag(Hilbert(FInt0));
        %EnvInt0 = sqrt(HFInt0.^2+FInt0.^2);
        EnvInt0 = HFInt0.^2+FInt0.^2;
        
        BatReconMod2 = EnvInt0;
        %bproc = bproc + length(BatReconMod2);

        pinic = round(0.120*Fam);
        pfin1 = length(BatReconMod2) - pinic;
        pfin2 = length(SinalModif2) - pinic;
        
        MaxSinal2 = max(BatReconMod2);
        MaxSinal = max(MaxSinal1,MaxSinal2);
        SinalModif2 = SinalModif2/MaxSinal;  %%% Normalizações
        BatReconMod2 = BatReconMod2/MaxSinal;  %%% Normalizações

    
        %limiarTemp = 0.5 * max(abs(BatReconMod2(pinic:pfin1)));
            
        limiarTemp = 0.3*max(BatReconMod2(pinic:pfin1));
        limiarTemp = FPlimiarEscala*max(BatReconMod2(pinic:pfin1));
        
    
                 
        %%% Busca de picos no intervalo filtrado %%%

   
         dmin = 0.120*Fam;
         PicoIdent = 0;
         PicoIdent = det_pico_sinal(pinic,SinalModif2,pfin2,limiarTemp,dmin);
         PontoFiducial=[];
         QminBat = floor(length(Intervalo)/Fam);
         QmaxBat = QminBat*2;
            
         if length(PicoIdent) > 0
             
            for iPM = 1:length(PicoIdent)
                PicoMapeado = Mapeia_Pico(PicoIdent(iPM),Intervalo);
                PontoFiducial = [PontoFiducial PicoMapeado];
            end
            
            if length(PontoFiducial) < QminBat |  length(PontoFiducial) > QmaxBat
                PontoFiducial = PontoFiducial(1);
            end
            
%             figure(3);
%             subplot(2,1,2);
%             %plot(SinalModif2);hold on;plot(PicoIdent,SinalModif2(PicoIdent),'ro');plot((1:length(SinalModif2)),limiarTemp*ones(1,length(SinalModif2)),'k');
%             plot((1:length(SinalModif2))/Fam + 2.3-offset1/Fam,SinalModif2);hold on;plot(PicoIdent/Fam + 2.3-offset1/Fam,SinalModif2(PicoIdent),'ro');plot((1:length(SinalModif2))/Fam + 2.3-offset1/Fam,limiarTemp*ones(1,length(SinalModif2)),'k');
%             subplot(2,1,1);
%             %plot(Intervalo);hold on;plot(PontoFiducial,Intervalo(PontoFiducial),'ro');
%             plot((1:length(Intervalo))/Fam + 2.3-offset1/Fam,Intervalo);hold on;plot(PontoFiducial/Fam + 2.3-offset1/Fam,Intervalo(PontoFiducial),'ro');
%             keyboard;
%             close(figure(3));            
            
            if PontoFiducial(1) ~= 0 & PontoFiducial(1) + inicio - 1 - IndexBat > dmin & abs(PontoFiducial(end) + inicio - 1 - IndexPico) > dmin
         
                for iPP = 1:length(PontoFiducial)
                    Pico = PontoFiducial(iPP) + inicio - 1;
                    Intt = Pico - IndexBat; 
                    MetricaTeste = (Intt - mediaIntervalo)/desvioIntervalo;
                    picos_R(length(picos_R)+1) = Pico;
                end
                         
            end
         end
            
        
    end
   
        if length(picos_R) == Tini
            fim = 0;
            %picos_R(length(picos_R)+1) = IndexPico;
        else fim = picos_R(length(picos_R))+10;
        end
        

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
       
    %%% Testar a possibilidade de se colocar offset1 como 0,120 s
	  
    offset1 = round((picos_R(end) - picos_R(end-1))/2);
    offset1 = round(0.120*Fam);
    offset2 = round((mean(IBB))/3);
    inicio = IndexBat + offset1;
    fim = IndexPico-offset1;
    Intervalo = ECG(inicio : fim) - m0;
    if length(Intervalo) > 0.090*Fam
        Wavelet_Filha = CHAPEU_MEXICANO(EscalaSel);
        SinalModif1 = conv(Wavelet_Filha,Intervalo);
        bproc = bproc + length(Intervalo);
        atraso = round(length(Wavelet_Filha)/2);
        SinalModif1 = SinalModif1(atraso:length(SinalModif1)-atraso);
        
        %FInt0 = DERIVADOR_SINAL(SinalModif1,2);  %%% Derivada do sinal filtrado
        FInt0 = diff(SinalModif1);  %%% Derivada do sinal filtrado
        
        %bproc = bproc + length(Intervalo);
        HFInt0 = imag(hilbert(FInt0));
        %bproc = bproc + length(Intervalo);
        %EnvInt0 = sqrt(HFInt0.^2+FInt0.^2); %%% Envelope do sinal analitico correspondente 
        EnvInt0 = HFInt0.^2+FInt0.^2;
        
        SinalModif2 = EnvInt0;
        MaxSinal1 = max(SinalModif2);
        %SinalModif2 = SinalModif2/MaxSinal;  %%% Normalizações

        IndexBatAnt = picos_R(length(picos_R)-1);
        BatRecon = ECG(IndexBatAnt - offset2 : IndexBatAnt + offset2)-m0;
        BatReconMod1 = conv(Wavelet_Filha,BatRecon);
        bproc = bproc + length(BatRecon);
        BatReconMod1 = BatReconMod1(atraso:length(BatReconMod1)-atraso);
        %FInt0 = DERIVADOR_SINAL(BatReconMod1,2);  %%% Derivada do sinal filtrado
        FInt0 = diff(BatReconMod1);  %%% Derivada do sinal filtrado
        %bproc = bproc + length(BatRecon);
        HFInt0 = imag(hilbert(FInt0));
        %EnvInt0 = sqrt(HFInt0.^2+FInt0.^2);
        EnvInt0 = HFInt0.^2+FInt0.^2;
        
        BatReconMod2 = EnvInt0;
        %bproc = bproc + length(BatReconMod2);
        
        MaxSinal2 = max(BatReconMod2);
        MaxSinal = max(MaxSinal1,MaxSinal2);
        SinalModif2 = SinalModif2/MaxSinal;  %%% Normalizações
        BatReconMod2 = BatReconMod2/MaxSinal;  %%% Normalizações


        pinic = round(0.120*Fam);
        pfin1 = length(BatReconMod2) - pinic;
        pfin2 = length(SinalModif2) - pinic;
        
    
        %limiarTemp = 0.5 * max(abs(BatReconMod2(pinic:pfin1)));
            
        limiarTemp = 0.3*max(BatReconMod2(pinic:pfin1));
        limiarTemp = FPlimiarEscala*max(BatReconMod2(pinic:pfin1));
        
    
                 
        %%% Busca de picos no intervalo filtrado %%%
        

   
         dmin = 0.120*Fam;
         PicoIdent = 0;
         PicoIdent = det_pico_sinal(pinic,SinalModif2,pfin2,limiarTemp,dmin);
         QminBat = floor(length(Intervalo)/Fam);
         QmaxBat = QminBat*2;
            
         PontoFiducial=[];   
            
         if length(PicoIdent) > 0
             
            for iPM = 1:length(PicoIdent)
                PicoMapeado = Mapeia_Pico(PicoIdent(iPM),Intervalo);
                PontoFiducial = [PontoFiducial PicoMapeado];
            end
            
            if length(PontoFiducial) < QminBat |  length(PontoFiducial) > QmaxBat
                PontoFiducial = PontoFiducial(1);
            end
            
%             figure(3);
%             subplot(2,1,2);
%             %plot(SinalModif2);hold on;plot(PicoIdent,SinalModif2(PicoIdent),'ro');plot((1:length(SinalModif2)),limiarTemp*ones(1,length(SinalModif2)),'k');
%             plot((1:length(SinalModif2))/Fam + 2.3-offset1/Fam,SinalModif2);hold on;plot(PicoIdent/Fam + 2.3-offset1/Fam,SinalModif2(PicoIdent),'ro');plot((1:length(SinalModif2))/Fam + 2.3-offset1/Fam,limiarTemp*ones(1,length(SinalModif2)),'k');
%             subplot(2,1,1);
%             %plot(Intervalo);hold on;plot(PontoFiducial,Intervalo(PontoFiducial),'ro');
%             plot((1:length(Intervalo))/Fam + 2.3-offset1/Fam,Intervalo);hold on;plot(PontoFiducial/Fam + 2.3-offset1/Fam,Intervalo(PontoFiducial),'ro');
%             keyboard;
%             close(figure(3));
            
            
            
            if PontoFiducial(1) ~= 0 & PontoFiducial(1) + inicio - 1 - IndexBat > dmin & abs(PontoFiducial(end) + inicio - 1 - IndexPico) > dmin
         
                for iPP = 1:length(PontoFiducial)
                    Pico = PontoFiducial(iPP) + inicio - 1;
                    Intt = Pico - IndexBat; 
                    MetricaTeste = (Intt - mediaIntervalo)/desvioIntervalo;
                    picos_R(length(picos_R)+1) = Pico;
                end
            end
         end
            
        
    end
   
        if length(picos_R) == Tini
            picos_R(length(picos_R)+1) = IndexPico;
        end
        fim = picos_R(length(picos_R))+10;
        

end
    

function[EscalaSeg]=Treino_OnOff_v5(NBatTreino)

global anot sfreq anot tipo sfreq ECG
global Se Pr Der FN FP TP QRSver pproc
global picos_R exame On Off EscalaSel exame

%%%% Rotina de Treinamento - resoluçao apropriada para o sinal Chapeu de
%%%% Mexicano


Fam = sfreq;
K = NBatTreino;
VetEsc = [8 16 32 64];
dminQRS = 0.085; % Limiar mínimo em milisegundos
dmaxQRS = 0.120; % Limiar máximo em milisegundos

offset0 = round(0.45*Fam);
offset1 = round(0.3*Fam); % para QRS onset
offset2 = round(0.19*Fam); % para QRS offset
Wmax = round(0.150*Fam); % Tamanho inicial da janela de busca
Wmax1 = round(0.140*Fam);
Wmax2 = round(0.140*Fam);
difoff1 = offset0-offset1;
difoff2 = offset0-offset2;

for t = 1 : length(VetEsc)

    Escala = VetEsc(t);
    Wavelet_Filha = CHAPEU_MEXICANO(Escala);
    mint = 0;
 
    for i = 1:K
        
        inicio(i,t) = max(picos_R(i) - offset0,1);
        fim = picos_R(i) + offset0;
        janela = ECG(inicio(i,t):fim);
                
        Fjan = conv(Wavelet_Filha,janela);
        atraso = round(length(Wavelet_Filha)/2);
        Fjan = Fjan(atraso:length(Fjan)-atraso);
        Fjan = diff(Fjan);  %%% Derivada do sinal filtrado
        HFjan = imag(hilbert(Fjan)); %%% Transformada de Hilbert da derivada do sinal filtrado
        Envjan = sqrt(HFjan.^2+Fjan.^2); %%% Envelope do sinal analitico correspondente 
        jfilt = Envjan;
        
        if i > 1
            jfilt = Envjan(difoff1:end-difoff2);
        else jfilt = Envjan(1:end-difoff2);
        end
  
        
        if inicio(i,t) == 1
            relR = picos_R(i);
        else relR = picos_R(i) - inicio(i,t)-difoff1 + 1; % Localização do ponto fiducial na janela
        end
        %pinic = round(0.100*sfreq); % Despreza-se amostras iniciais e finais da janela
        %pfin = length(Envjan) - pinic;
       
   
     %% Pesquisa de Pontos Criticos %%
        Wmax=Wmax1;
        picoR_map = find(jfilt == max(jfilt));
        picoR_map = picoR_map(1);
        
        %picoR_map = relR;
        
        index = picoR_map;
        On = 1; % Valores iniciais para as localizações ONSET e OFFSET
        Off = length(jfilt);
        W = Wmax; % Etapa 1: Tamanho da janela = Tamanho máximo
        idArea = 1; % Indexador do vetor Area
        gap = picoR_map; % Localização inicial da extremidade direita da janela
        idjan = max(1,gap-W+1); % Localização inicial da extremidade esquerda da janela
        gap = idjan+W-1;
        SomArea = 0;
        
        while gap(end) <= length(jfilt)  % Localização final da extremidade esquerda da janela
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan < gap(end) % Calculo da área correspondente ao intervalo [t-W,t], abaixo do sinal e acima da linha horizontal (t, Env(t))
                SomArea(idArea) = SomArea(idArea) + (jfilt(iddjan)-jfilt(gap(end)))*(1/Fam);
                iddjan = iddjan+1;
            end
            
            idArea = idArea + 1;
            idjan = idjan + 1; % Deslocamento da janela para a esquerda
            gap = [gap gap(end)+1];
            
        end
        
        MaxArea = find(SomArea == max(SomArea)); % Identificando a posição da extremidade direita da janela correspondente à área máxima
        OffTemp = gap(MaxArea(1));
    
        
        % Etapa 2: Tamanho da janela = Tamanho adequado ao QRS em exame
        
        idjan = picoR_map; 
        W = OffTemp - idjan;
        gap = idjan + W - 1;
        
        SomArea = 0;
        idArea = 1;
        
        while idjan >= 1
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan < gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt(iddjan)-jfilt(gap(end)))*(1/Fam);
                iddjan = iddjan+1;
            end
            
            idArea = idArea + 1;
            idjan = idjan - 1;
            gap = [gap gap(end)-1];
        end
        
        MaxArea = find(SomArea == max(SomArea));
        Off = gap(MaxArea(1));
        %keyboard;
        
        
        Wmax=Wmax2;
        idArea = 1;
        gap = picoR_map; % Posição inicial da extremidade esquerda da janela
        W = Wmax;
        idjan = min(gap + W - 1,length(jfilt)); % Posição inicial da extremidade direita da janela
        gap = idjan - W + 1;
        SomArea = 0;
        
       
        while gap(end) >= 1 && idjan >= picoR_map
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan > gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt(iddjan)-jfilt(gap(end)))*(1/Fam);
                iddjan = iddjan-1;
            end
            
            idArea = idArea + 1;
            idjan = idjan - 1; % Deslocamento da janela para a direita
            gap = [gap gap(end)-1];
        end
            
        MaxArea = find(SomArea == max(SomArea)); % Identificando a posição da extremidade esquerda da janela correspondente à área máxima
        OnTemp = gap(MaxArea(1));
        
        
        
        W = picoR_map - OnTemp; % Tamanho da Janela = Tamanho adequado ao QRS
        idjan = picoR_map;
        gap = idjan - W + 1;
        
        SomArea = 0;
        idArea = 1;
        
        while idjan <= length(jfilt) && gap(end)<= picoR_map
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan > gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt(iddjan)-jfilt(gap(end)))*(1/Fam);
                iddjan = iddjan-1;
            end
                                  
            idArea = idArea + 1;
            idjan = idjan + 1;
            gap = [gap gap(end)+1];
        end
            
        MaxArea = find(SomArea == max(SomArea));
        On = max(1,gap(MaxArea(1)));
        
            
        metrica1(i,t) = (Off - On)/sfreq;
        metrica2(i,t) = sum(janela(On:Off).^2);
        metrica3(i,t) = std(janela(On:Off));
        metrica4(i,t) = (Off + On -2*picoR_map)/sfreq;
        
        janela = janela(difoff1:end-difoff2);
        
%         if i == 1
%             
%             figure(4);
%             subplot(2,1,1);
%             plot((1:length(jfilt))/Fam,jfilt);
%             hold on; plot(On/Fam,jfilt(On),'ro');plot(Off/Fam,jfilt(Off),'ro');plot(picoR_map/Fam,jfilt(picoR_map),'ok');
%             subplot(2,1,2);
%             plot((1:length(janela))/Fam,janela);
%             hold on; plot(On/Fam,janela(On),'ro');plot(Off/Fam,janela(Off),'ro');plot(picoR_map/Fam,janela(picoR_map),'ok');
%             
%             close(figure(4));
%         end
        
     
    
    
     
     
    end

      merro(t) = mean(metrica1(:,t));
      merro4(t) = mean(metrica4(:,t));
      
end

% figure(5);
% plot(VetEsc,merro);
% figure(6);
% plot(VetEsc,merro4);
vbol1 = 0;




teste1 = find(merro > 0.15);
merro_mod = merro(2:end)-merro(1:end-1); 
teste2 = find(merro_mod > 0);

if length(teste1) < length(merro) & length(teste2) == length(merro)-1
    teste3 = find(merro >= 0.090);  %% Idéia: Se teste3(1) >= 3, então usar 1 escala para onset e 1 escala para offset
    if teste3(1) <= 2 & (abs(merro_mod(2)) < 10*abs(merro_mod(1))) & (abs(merro_mod(3)) < 10*abs(merro_mod(2)))
        EscalaSeg = VetEsc(teste3(1));
    else   if  abs(merro_mod(2)) > 2*abs(merro_mod(1)) | abs(merro_mod(3)) > 2*abs(merro_mod(2))
                EscalaSeg = [VetEsc(2) VetEsc(4)];
           else EscalaSeg = VetEsc(2);
           end
    end
else if length(teste2) < length(merro)-1
        
          if  abs(merro_mod(2)) > abs(10*merro_mod(1)) | abs(merro_mod(3)) > abs(10*merro_mod(2))
              
              if merro(4) > merro(3)
                EscalaSeg = [VetEsc(2) VetEsc(2)];
              else EscalaSeg = [VetEsc(2) VetEsc(3)];
              end
          
          else EscalaSeg = VetEsc(4);
              
          end

    else EscalaSeg = VetEsc(1);
    end
end




%%% Aplicando a escala de Wavelet adequada para segmentar todos os
%%% complexos QRS 


if length(EscalaSeg) == 1
    Wavelet_Filha1 = CHAPEU_MEXICANO(EscalaSeg);
else
    Wavelet_Filha1 = CHAPEU_MEXICANO(EscalaSeg(1));
    Wavelet_Filha2 = CHAPEU_MEXICANO(EscalaSeg(2));
end
mint = 0;
 

for ibat=1:NBatTreino
    
                   
       inicio(ibat) = max(picos_R(ibat) - offset0,1);
       fim = picos_R(ibat) + offset0;
       janela = ECG(inicio(ibat):fim);
       Fjan = conv(Wavelet_Filha1,janela);
       atraso = round(length(Wavelet_Filha1)/2);
       Fjan = Fjan(atraso:length(Fjan)-atraso);
       Fjan = diff(Fjan);  %%% Derivada do sinal filtrado
       HFjan = imag(hilbert(Fjan)); %%% Transformada de Hilbert da derivada do sinal filtrado
       Envjan = sqrt(HFjan.^2+Fjan.^2); %%% Envelope do sinal analitico correspondente 
       %Envjan = HFjan.^2+Fjan.^2;
       jfilt1 = Envjan;
       if ibat > 1
        jfilt1 = Envjan(difoff1:end-difoff2);
       else jfilt1 = Envjan(1:end-difoff2);
       end
        
       if length(EscalaSeg) == 1
           jfilt2 = jfilt1;
       else 
           Fjan = conv(Wavelet_Filha2,janela);
           atraso = round(length(Wavelet_Filha2)/2);
           Fjan = Fjan(atraso:length(Fjan)-atraso);
           Fjan = diff(Fjan);  %%% Derivada do sinal filtrado
           HFjan = imag(Hilbert(Fjan)); %%% Transformada de Hilbert da derivada do sinal filtrado
            Envjan = sqrt(HFjan.^2+Fjan.^2); %%% Envelope do sinal analitico correspondente 
            %Envjan = HFjan.^2+Fjan.^2;
            jfilt2 = Envjan;
            if ibat > 1
                jfilt2 = Envjan(difoff1:end-difoff2);
            else jfilt2 = Envjan(1:end-difoff2);
            end
           
           
       end
       
       
       
        if inicio(ibat)==1
            relR = picos_R(ibat);
        else relR = picos_R(ibat) - inicio(ibat)-difoff1 + 1;
        end


        %keyboard;    
     %% Pesquisa de Pontos Criticos %%
     
        %picoR_map = relR;
        picoR_map = find(jfilt2 == max(jfilt2));
        picoR_map = picoR_map(1);
        
        %picoR_map = relR;
        
        index = picoR_map;
        On(ibat) = 1;
        Off(ibat) = length(jfilt2);
        W = Wmax1;
        idArea = 1;
        idjan = picoR_map;
        gap = min(idjan + W - 1,length(jfilt2));
        SomArea = 0;
        
        %keyboard;
        
        while idjan >= 1
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan < gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt2(iddjan)-jfilt2(gap(end)))*(1/Fam);
                iddjan = iddjan+1;
                
            end
            
            idArea = idArea + 1;
            idjan = idjan - 1;
            gap = [gap gap(end)-1];
        end
        
        
        
        vtt = find(SomArea == max(SomArea));
        MaxArea = vtt;
        OffTemp = gap(MaxArea(1));
        
        
        %keyboard;
     
        idjan = picoR_map;
        
        %W = OffTemp - idjan + 4;
        W = OffTemp - idjan;
        
        %keyboard;
        
        gap = min(idjan + W - 1,length(jfilt2));
        SomArea = 0;
        idArea = 1;
        while idjan >= 1
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan < gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt2(iddjan)-jfilt2(gap(end)))*(1/Fam);
                iddjan = iddjan+1;
            end
            
            idArea = idArea + 1;
            idjan = idjan - 1;
            gap = [gap gap(end)-1];
        end
     
        vtt = find(SomArea == max(SomArea));
        MaxArea = vtt;
        Off(ibat) = gap(MaxArea(1));
      
        
        
        picoR_map = find(jfilt1 == max(jfilt1));
        picoR_map = picoR_map(1);
        
        %picoR_map = relR;
        
        idArea = 1;
        idjan = picoR_map;
        W = Wmax2;
        gap = max(idjan - W + 1,1);
        idjan = gap+W-1;
        SomArea = 0;
        
       
        while idjan <=length(jfilt1) & gap(end) <= picoR_map
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan > gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt1(iddjan)-jfilt1(gap(end)))*(1/Fam);
                iddjan = iddjan-1;
            end
            
            idArea = idArea + 1;
            idjan = idjan + 1;
            gap = [gap gap(end)+1];
        end
            
        MaxArea = find(SomArea == max(SomArea));
        OnTemp = gap(MaxArea(1));
        
        
        idjan = picoR_map;
        W = idjan - OnTemp;
        gap = idjan - W + 1;
        SomArea = 0;
        idArea = 1;
        
        while idjan <= length(jfilt1) && gap(end) <= picoR_map
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan > gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt1(iddjan)-jfilt1(gap(end)))*(1/Fam);
                iddjan = iddjan-1;
            end
            
            idArea = idArea + 1;
            idjan = idjan + 1;
            gap = [gap gap(end)+1];
        end
            
        MaxArea = find(SomArea == max(SomArea));
        On(ibat) = gap(MaxArea(1));        
        
        
        if ibat > 1
            On(ibat) = On(ibat) + inicio(ibat) +difoff1 - 1;
            Off(ibat) = Off(ibat) + inicio(ibat)+difoff1 -1;
        end
        
%         figure(2);
%         subplot(2,1,1);
%         plot((1:length(janela))/Fam,janela);
%         hold on;
%         plot((Off(ibat)-inicio(ibat)+1)/Fam,janela(Off(ibat)-inicio(ibat)+1),'ro');
%         plot((On(ibat)-inicio(ibat)+1)/Fam,janela(On(ibat)-inicio(ibat)+1),'ro');
%         subplot(2,1,2);
%         plot((1:length(jfilt2))/Fam,jfilt2);
%         hold on;
%         if ibat > 1
%             plot((Off(ibat)-inicio(ibat)-difoff1+1)/Fam,jfilt2(Off(ibat)-inicio(ibat)-difoff1+1),'ro');
%             plot((On(ibat)-inicio(ibat)-difoff1+1)/Fam,jfilt2(On(ibat)-inicio(ibat)-difoff1+1),'ro');
%         else plot((Off(ibat)-inicio(ibat)+1)/Fam,jfilt2(Off(ibat)-inicio(ibat)+1),'ro');
%             plot((On(ibat)-inicio(ibat)+1)/Fam,jfilt2(On(ibat)-inicio(ibat)+1),'ro');
%         end
%         keyboard;
%         close(figure(2));
        
end
       



% CM = CHAPEU_MEXICANO(Escala);
% DCM = Derivador(CM,1);
% 
% 
% 
% for i = 1:K
%    QRSi(i) = On(i,Iescala(1))+inicio(i,Iescala(1))-1;
%    QRSf(i) = Off(i,Iescala(1))+inicio(i,Iescala(1))-1;
%    eQS(i) = ErroEnergia(CM,exame(QRSi(i):QRSf(i)),DCM,i);
% end
 
%PesoHilbert = Pmax(Iescala);
 


 
    

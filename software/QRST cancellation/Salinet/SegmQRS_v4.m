function[Resultado] = SegmQRS_v3(ibat,escala)

global sfreq picos_R exame On Off PesoHilbert ECG

if length(escala) == 1
    Wavelet_Filha1 = CHAPEU_MEXICANO(escala);
else Wavelet_Filha1 = CHAPEU_MEXICANO(escala(1));
     Wavelet_Filha2 = CHAPEU_MEXICANO(escala(2));
end
dminQRS = 0.085;
dmaxQRS = 0.120;
offset0 = round(0.45*sfreq);
offset1 = round(0.3*sfreq);
offset2 = round(0.19*sfreq);
Resultado = 0;
W0 = round(0.02*sfreq);
Wmax = round(0.15*sfreq);
Wmax1 = round(0.140*sfreq);
Wmax2 = round(0.140*sfreq);
Fam=sfreq;
difoff1 = offset0-offset1;
difoff2 = offset0-offset2;

               
inicio(ibat) = picos_R(ibat) - offset0;
fim = min(picos_R(ibat) + offset0,length(ECG));
janela = ECG(inicio(ibat):fim);
                            
if length(janela) > dminQRS*sfreq
                
                    Fjan = conv(Wavelet_Filha1,janela);
                    atraso = round(length(Wavelet_Filha1)/2);
                    Fjan = Fjan(atraso:length(Fjan)-atraso);
                    Fjan = diff(Fjan);  %%% Derivada do sinal filtrado
                    HFjan = imag(hilbert(Fjan)); %%% Transformada de Hilbert da derivada do sinal filtrado
                    Envjan = sqrt(HFjan.^2+Fjan.^2); %%% Envelope do sinal analitico correspondente 
                    jfilt1 = Envjan;
                    jfilt1 = Envjan(difoff1:end-difoff2);
                    
                    
                    if length(escala) == 1
                        jfilt2 = jfilt1;
                    else 
                        
                        Fjan = conv(Wavelet_Filha2,janela);
                        atraso = round(length(Wavelet_Filha2)/2);
                        Fjan = Fjan(atraso:length(Fjan)-atraso);
                        Fjan = diff(Fjan);  %%% Derivada do sinal filtrado
                        HFjan = imag(Hilbert(Fjan)); %%% Transformada de Hilbert da derivada do sinal filtrado
                        Envjan = sqrt(HFjan.^2+Fjan.^2); %%% Envelope do sinal analitico correspondente 
                        jfilt2 = Envjan;
                        jfilt2 = Envjan(difoff1:end-difoff2);

                    end
                    
end



        relR = picos_R(ibat) - inicio(ibat)-difoff1 + 1;


    
     %% Pesquisa de Pontos Criticos %%
     
        picoR_map = find(jfilt2 == max(jfilt2));;
        picoR_map = picoR_map(1);
        
        %picoR_map = relR;
        
        index = picoR_map;
        On(ibat) = 1;
        Off(ibat) = length(jfilt2);
        W = Wmax1;
        idArea = 1;
        gap = picoR_map;
        idjan = max(1,gap-W+1);
        gap = idjan+W-1;
        SomArea = 0;
        
        while gap(end) <= length(jfilt2)
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan < gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt2(iddjan)-jfilt2(gap(end)))*(1/sfreq);
                iddjan = iddjan+1;
            end
            
            idArea = idArea + 1;
            idjan = idjan + 1;
            gap = [gap gap(end)+1];
        end
        
        vtt = find(SomArea == max(SomArea));
        MaxArea = vtt;
        OffTemp = gap(MaxArea(1));
        
               
        idjan = picoR_map;
        W = OffTemp - idjan;
        gap = min(idjan + W - 1,length(jfilt2));
        SomArea = 0;
        idArea = 1;
        
        while idjan >= 1
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan < gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt2(iddjan)-jfilt2(gap(end)))*(1/sfreq);
                iddjan = iddjan+1;
            end
            
            idArea = idArea + 1;
            idjan = idjan - 1;
            gap = [gap gap(end)-1];
        end
        
        MaxArea = find(SomArea == max(SomArea));
        Off(ibat) = gap(MaxArea(1));
        
        %keyboard;
        
         picoR_map = find(jfilt1 == max(jfilt1));
         picoR_map = picoR_map(1);
         
         %picoR_map = relR;
        
        gap=0;
        idArea = 1;
        gap = picoR_map;
        W = Wmax2;
        idjan = min(gap + W - 1,length(jfilt1));
        gap = idjan - W + 1;
        SomArea = 0;
        
        while gap(end) >= 1
            
            iddjan = idjan;
            SomArea(idArea) = 0;
            while iddjan > gap(end)
                SomArea(idArea) = SomArea(idArea) + (jfilt1(iddjan)-jfilt1(gap(end)))*(1/Fam);
                iddjan = iddjan-1;
            end
            
            idArea = idArea + 1;
            idjan = idjan - 1;
            gap = [gap gap(end)-1];
        end
            
        MaxArea = find(SomArea == max(SomArea));
        OnTemp = gap(MaxArea(1));
        
        %keyboard;
        
        idjan = picoR_map;
        W = idjan - OnTemp;
        gap = max(idjan - W + 1,1);
        SomArea = 0;
        idArea = 1;
        
        while idjan <= length(jfilt1)
            
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
%         if abs((Off(ibat)-On(ibat))/sfreq) < 0.1
%             keyboard;
%         end
%         if ibat == 142
%             keyboard;
%         end
%          figure(3);
%          plot(jfilt1);
%          pause;
%          close(figure(3));
               
        On(ibat) = On(ibat) + inicio(ibat)+difoff1 - 1;
        Off(ibat) = Off(ibat) + inicio(ibat)+difoff1 -1;
        Resultado = 1;
        
        
        
        
        
        
        
        
        
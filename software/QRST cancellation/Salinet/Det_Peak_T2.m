function[Tp]=Det_Peak_T2

global anot sfreq anot tipo sfreq ECG
global Se Pr Der FN FP TP QRSver pproc
global picos_R exame On Off EscalaSel exame picos_Rlst picos_Ri
global Onlb Offlb offset0 spikev Tp T Tendv


Tp=[];
Tendv=[];

SD1 = 0.1:0.05:2;
N0 = 10;

%alfa = pi/20-[0:pi/(N0*20):pi/20];
offset0 = round(0.100*sfreq);
%offset1 = round(0.090*sfreq);
%mSTD = zeros(length(SD1),length(alfa));
%vis = 0;

for iW = 1:length(Off)
    
%     mSTD = zeros(length(SD1),length(alfa));   
    if iW < length(Off)
        sW=ECG(Off(iW)-offset0:On(iW+1)+offset0);
    else sW=ECG(Off(iW)-offset0:min(Off(iW)+round(0.5*sfreq)+offset0,length(ECG)));
         
    end
%     
%     for iSD1 = 1:length(SD1)
%         
%         for ialfa = 1:length(alfa)
%             
%             g = Gaussia_Distortion_v4(SD1(iSD1),alfa(ialfa));
%             if length(g) == 1
%                 mSTD(iSD1,ialfa) = inf;
%             else 
%                 sft = conv(sW,g);
%                 atraso = round(length(g)/2);
%                 sft = sft(atraso:length(sft)-atraso);
%         
%                 
%                 sft = sft(offset0:end-offset0);
%                 
%             
%                 lim = 0.7*max(abs(sft));   
%                 j1=1;
%                 j2=0;
%                 j3=length(sft);
%                 
%                 while (sft(j1) > lim | sft(j3) > lim) & j1+5 < j3 
%                 %ROTINA PARA EVITAR A DETECCAO DE PICOS AO INICIO E/OU AO FINAL
%                     j1=j1+1;
%                     j3=j3-1;
%                     lim = 0.7*max(abs(sft(j1:j3)));
%                 end
%                 
%                 pck=0;
%                 cpck=0;
%             
%                 while j1 < length(sft)
%             
%                     if abs(sft(j1)) > lim 
%                 
%                         j2=j1;
%                 
%                         while (abs(sft(j2)) > lim & j2+1 <= length(sft))
%                             j2=j2+1;
%                         end
%                 
%                         jan = sft(j1:j2);
%                         pckt = find(abs(jan)==max(abs(jan)));
%                         cpck=cpck+1;
%                         pck(cpck)=pckt(1)+j1-1;
%                         break       
%                 
%                     end
%             
%                     j1=j1+1;
%             
%                 end
%                 if pck(1) == 0
%                     pck(1)=1;
%                 end
%             
                sW=sW-mean(sW);
                f1=0.01;f2=12;Wn=[2*f1/sfreq 2*f2/sfreq];order=min(32,floor(length(sW)/4));[b,a]=fir1(order,Wn);fsW = filtfilt(b,a,sW);sW=fsW;  
                signal = sW(offset0:end-offset0);
%                 [gmod,eg] = normalizag_3(g,signal,pck(1));  %% AO INVÉS DE PASSAR O PRIMEIRO PICO, PASSAR O DE MAIOR AMPLITUDE...
%                 irg = pck(1)-round(length(g)/2);
%             
%                 mSTD(iSD1,ialfa) = eg;
%            
%             
%                      
%                 if vis == 1 & length(Tp)==2
%             
%                     figure(3);
%              
%                     if signal(pck(1)) < 0
%                         gmod = gmod*(-1);
%                     end
%                     gpico = find(abs(gmod)==max(abs(gmod)));
%                     irg = pck(1)-gpico;
%                     pg = pck(1);
%                     subplot(3,1,1);
%                     plot((1:length(signal))/sfreq,signal)
%                     hold on;
%                     xini = max(irg,1);
%                     xfin = irg+length(g)-1;
%                     xinig = xini-irg + 1;
%                     yfing = length(gmod);
%                     plot((xini:xfin)/sfreq,gmod(xinig:yfing),'r');
%                     plot(pg/sfreq,signal(pg),'ro');
%                     title('Original Signal Window');
%                     axis([0 length(signal)/sfreq min(gmod) max(gmod)])
%                     subplot(3,1,2);plot((1:length(sft))/sfreq,sft);hold on;plot((1:length(sft))/sfreq,lim*ones(1,length(sft)),'k');plot((1:length(sft))/sfreq,-lim*ones(1,length(sft)),'k');
%                     title('Filtered Signal Window');
%                     subplot(3,1,3);
%                     plot((1:length(g))/sfreq,g);
%                     pause;
%                     close(figure(3));
%             
%                 end
%             
%             
%             end
%         
%         end
%         
%     end
%     
%     il=1;
%     ic=1;
%     
%     S=sort(mSTD);
%     MinV = S(1,:);
%     vmin = min(MinV);
%     Lsel = 0;
%     Csel = 0;
%     Vfind = 0;
%     
%     while il<=length(SD1)
%         ic = 1;
%         while ic<=length(alfa)
%             if mSTD(il,ic) == vmin
%                 Lsel = il;
%                 Csel = ic;
%                 Vfind = 1;
%                 break;
%             end
%             ic = ic+1;
%         end
%         if Vfind == 1
%             break;
%         end
%         il = il+1;
%     end
%     
%    
%     
%     
%    %% Applying selected standard deviation and deflection
%    if Vfind == 1
%    
%   
%       
%     g = Gaussia_Distortion_v4(SD1(il),alfa(ic));
%     sft = conv(sW,g);
%     atraso = round(length(g)/2);
%     sft = sft(atraso:length(sft)-atraso);
%         
%    
%     sft = sft(offset0:end-offset0);
%             
%     lim = 0.7*max(abs(sft));   
%     j1=1;
%     j2=0;
%     j3=length(sft);
%                 
%     while (sft(j1) > lim | sft(j3) > lim) & j1+5 < j3 
%         %ROTINA PARA EVITAR A DETECCAO DE PICOS AO INICIO E/OU AO FINAL
%         j1=j1+1;
%         j3=j3-1;
%         lim = 0.7*max(abs(sft(j1:j3)));
%     end
%     
%     pck=0;
%     cpck=0;
%             
%    while j1 < length(sft)
%             
%         if abs(sft(j1)) > lim 
%                 
%                     j2=j1;
%                 
%                     while (abs(sft(j2)) > lim & j2+1 <= length(sft))
%                         j2=j2+1;
%                     end
%                 
%                     jan = sft(j1:j2);
%                     pckt = find(abs(jan)==max(abs(jan)));
%                     cpck=cpck+1;
%                     pck(cpck)=pckt(1)+j1-1;
%                     break       
%                 
%         end
%             
%         j1=j1+1;
%             
%    end
%     
%    signal = sW(offset0:end-offset0);  
%    [gmod,eg] = normalizag_3(g,signal,pck(1));
%    
%    gmod2 = gmod+100; 
%    gpico = find(abs(gmod2)==max(abs(gmod2)));         
%     if signal(pck(1)) < 0
%         gmod = gmod*(-1);
%     end
%     
%    
%     gpico = gpico(1);
%     %gpico = round(length(gmod)/2);
%     
%     
%     %% Para correcao de defasagem na deteccao do pico %%
%     
%     lsrc1 = round(0.020*sfreq);
%     lsrc2 = round(0.020*sfreq);
%     s1=gmod(max(gpico-lsrc1,1):min(gpico+lsrc1,length(gmod)));
%     s2=signal(max(pck(1)-lsrc2,1):min(pck(1)+lsrc2,length(signal)));
%     L = length(s1);
%     cc=xcorr(s1,s2);
%     pmax = find(cc==max(cc));
%     pmaxc = pmax-lsrc1;
%     
%     if pck(1)-lsrc2 > 1
%         pck = pck(1)-lsrc2 + pmax-lsrc1; 
%         Tp(iW) = pck(1)+Off(iW)-1;
%         
%     else pck = pmax-lsrc1;
%          Tp(iW) = pck(1)+Off(iW)-1;
%     end
    
%   index = pck(1);


    p = 0.016*sfreq;
    w = 0.100*sfreq;
    k = p+1;
    lambda = 10;
    sW2 = signal;
    
    if length(sW2) >= w*1.5
    gap = w+1;
    idjan = gap-w;
    SomArea = 0; idArea = 1;
    while gap(end) <= length(sW2)-p-1
        
                    iddjan = idjan;
                    SomArea(idArea) = 0;
                    sk = mean(sW2(gap(end)-p:gap(end)+p));
                    while iddjan <= gap(end)
                        SomArea(idArea) = SomArea(idArea) + (sW2(iddjan)-sk)*(1/sfreq);
                        iddjan = iddjan+1;
                    end
            
                    idArea = idArea + 1;
                    idjan = idjan + 1;
                    gap = [gap gap(end)+1];
    end
        
    k1 = find(SomArea == max(SomArea));
    vbk1 = 0;
    
    if k1 == 1 | k1 == length(SomArea)
        k1=3;
        while k1 <= length(SomArea)-3
            if SomArea(k1) > SomArea(k1-2) & SomArea(k1) > SomArea(k1+2)
                vbk1 = 1;
                break;
            else k1 = k1+1;
            end
        end
        
    else vbk1 = 1;
    end
    
    k2 = find(SomArea == min(SomArea));
    vbk2 = 0;
    
    if k2 == 1 | k2 == length(SomArea)
        k2=3;
        while k2 <= length(SomArea)-3
            if SomArea(k2) < SomArea(k2-2) & SomArea(k2) < SomArea(k2+2)
                vbk2 = 1;
                break;
            else k2 = k2+1;
            end
        end
    else vbk2 = 1;
    end
    
    
    
    if (abs(SomArea(k1))/abs(SomArea(k2))) > (1/lambda) & (abs(SomArea(k1))/abs(SomArea(k2))) < lambda
        if vbk1 == 1 & vbk2 == 1
            vb_bif = 1;
            kend = max(k1,k2);
            kend = gap(kend);
        else
            vb_bif = 0;
            if abs(SomArea(k1)) >= abs(SomArea(k2))
                 kend = gap(k1);
            else kend = gap(k2);
            end
            
        end
    else vb_bif = 0;
         if abs(SomArea(k1)) >= abs(SomArea(k2))
             kend = gap(k1);
         else kend = gap(k2);
         end
             
    end


    Tend = kend;
        
    
    %Tp(iW) = pck(1)+Off(iW)-1;
   
    Tendv(iW) = min(Tend+Off(iW)-1,length(ECG));
    
   %else keyboard;
    
   end
    
    
  
end

nz = find(Tendv ~= 0);
Tendv = round(Tendv(nz));
%T=Tp;

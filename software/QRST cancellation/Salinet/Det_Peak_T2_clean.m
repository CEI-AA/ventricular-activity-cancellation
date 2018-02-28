function[Tp]=Det_Peak_T2

global anot sfreq anot tipo sfreq ECG
global Se Pr Der FN FP TP QRSver pproc
global picos_R exame On Off EscalaSel exame picos_Rlst picos_Ri
global Onlb Offlb offset0 spikev Tp T Tendv

Tp=[];
Tendv=[];
offset0 = round(0.100*sfreq);

for iW = 1:length(Off)
    
    if iW < length(Off)
        sW=ECG(Off(iW)-offset0:On(iW+1)+offset0);
    else sW=ECG(Off(iW)-offset0:min(Off(iW)+round(0.5*sfreq)+offset0,length(ECG)));
         
    end
                sW=sW-mean(sW);
                f1=0.01;f2=12;Wn=[2*f1/sfreq 2*f2/sfreq];order=min(32,floor(length(sW)/4));[b,a]=fir1(order,Wn);fsW = filtfilt(b,a,sW);sW=fsW;  
                signal = sW(offset0:end-offset0);

    p = 0.016*sfreq;
    w = 0.100*sfreq; % 
    k = p+1;
    lambda = 10;
    sW2 = signal;
    
    if length(sW2) >= w*1.5 % 180 ms (
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
    
    %% Testing if is Monofasic or Bifasic
    if k1 == 1 || k1 == length(SomArea)
        k1=3;
        while k1 <= length(SomArea)-3
            if SomArea(k1) > SomArea(k1-2) && SomArea(k1) > SomArea(k1+2)
                vbk1 = 1;
                break;
            else k1 = k1+1;
            end
        end
        
    else vbk1 = 1;
    end
    
    k2 = find(SomArea == min(SomArea));
    vbk2 = 0;
    
    if k2 == 1 || k2 == length(SomArea)
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
    Tendv(iW) = min(Tend+Off(iW)-1,length(ECG));
   
   end
end

nz = find(Tendv ~= 0);
Tendv = round(Tendv(nz));


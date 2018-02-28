function[Tp]=Det_Peak_T3

global anot sfreq anot tipo sfreq ECG
global Se Pr Der FN FP TP QRSver pproc
global picos_R exame On Off EscalaSel exame picos_Rlst picos_Ri
global Onlb Offlb offset0 spikev Tp T Tendv


Tp=[];
Tendv=[];

SD1 = 0.1:0.05:2;
N0 = 10;

alfa = pi/20-[0:pi/(N0*20):pi/20];
offset0 = round(0.100*sfreq);
offset1 = round(0.090*sfreq);
mSTD = zeros(length(SD1),length(alfa));
vis = 0;

for iW = 1:length(Off)
    
    mSTD = zeros(length(SD1),length(alfa));   
    if iW < length(Off)
        sW=ECG(Off(iW)-offset0:On(iW+1)-0.010*sfreq+offset0);
    else sW=ECG(Off(iW)-offset0:min(Off(iW)+round(0.5*sfreq)+offset0,length(ECG)));
    end
    
    
    sW2 = sW;
    sW2 = sW2-mean(sW2);
    
    p = 0.016*sfreq;
    w = 0.100*sfreq;
    k = p+1;
    lambda = 30;
    
    
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
    
    Tendv(iW) = Tend+Off(iW)-1;
    
       
   end
    
    
  
  
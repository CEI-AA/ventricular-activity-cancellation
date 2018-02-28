function[On,Off] = DetPcrt(janela,relR,sinal)

global exame picos_R sfreq


Fam = sfreq;
k1 = length(janela);
j = 1;
Pcrt = 0;
On = 0;
Off = 0;
dmax = Fam*0.150;
lim = 0.10*max(abs(janela));

%% Procura por pontos críticos %%

for i = 1:k1-2
   if ((janela(i+1)-janela(i))*(janela(i+2)-janela(i+1)) < 0) & abs(i-relR) < dmax & abs(janela(i+1)) > lim
      Pcrt(j) = i+1;
      j = j+1;
   end
end

%%% Pcrt - vetor de pontos críticos na ordem temporal

%% Classificação dos pontos críticos %%

if Pcrt ~= 0
	k2 = length(Pcrt);
	PcrtC = 0;
	index=1;
	vtemp = Pcrt(1);
	j=1;
    
    %%% PcrtC - Classificação dos pontos criticos do de menor amplitude para o de maior amplitude.

   janelatemp = janela;
	while j <= k2
		for i = 1:k2
   	        if abs(janela(Pcrt(i))) > abs(janela(vtemp))
         	    vtemp = Pcrt(i);
         	    index(j) = i; %%% index - posições no vetor Pcrt
      	    end
   	    end
   	    PcrtC(j) = vtemp;
        janela(vtemp)=0;
   	    j = j+1;
	end
   
    incQ = relR-5;
    incS = relR+5;
    for i=1:length(Pcrt)
        if Pcrt(i) < relR
            incQ = Pcrt(i);
        end
    end
    i=length(Pcrt);
    while i >=1
        if Pcrt(i) > relR
            incS = Pcrt(i);
        end
        i=i-1;
    end
    
    janela = janelatemp;
	j=1;
	while j <= length(PcrtC)
        if index(j)+1 <= k2
            if Pcrt(index(j)+1) < relR & janela(PcrtC(j))*janela(incQ) < 0 
                On = PcrtC(j);
                break;
            end
        end
        j=j+1;
	end

	j=1;
	while j <= length(PcrtC)
        if index(j)-1 >=1
            if Pcrt(index(j)-1) > relR & janela(PcrtC(j))*janela(incS) < 0 
                Off = PcrtC(j);
                break;
            end
        end
        j=j+1;
	end
    


	if On == 0
		j=1;
		while j <= length(PcrtC)
            if Pcrt(index(j)) < relR
                On = PcrtC(j);
                break;
            end
            j=j+1;
        end
   
	end

	if On == 0
    	On = 1;
	end



	if Off == 0
		j=1;
		while j <= length(PcrtC)
            if Pcrt(index(j)) > relR
                Off = PcrtC(j);
                break;
            end
            j=j+1;
        end
   
	end

	if Off == 0
    	Off = length(janela);
	end

else 
    
    On=1;
    Off = length(janela);
end
 

 


% subplot(2,1,1);
% plot((1:length(sinal))/Fam,sinal);
% hold on;
% plot(relR/Fam,sinal(relR),'ro');
% plot(On/Fam,sinal(On),'ok');
% plot(Off/Fam,sinal(Off),'ok');
% subplot(2,1,2);
% plot((1:length(janela))/Fam,janela);
% hold on;
% plot((1:length(janela))/Fam,lim*ones(1,length(janela)),'k');
% plot((1:length(janela))/Fam,-lim*ones(1,length(janela)),'k');
% plot(relR/Fam,janela(relR),'ro');
% plot(On/Fam,janela(On),'ok');
% plot(Off/Fam,janela(Off),'ok');
% pause;
% clf;
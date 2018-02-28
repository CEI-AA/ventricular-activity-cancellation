function[] = atualiza_limiar
global limiar ECG picos_R m0

Bat = length(picos_R);
FatAmp = 0.7;
PicoEst = limiar/FatAmp;
% VPico = ECG(picos_R(Bat)- 5:min(picos_R(Bat)+ 5,length(ECG)))-m0;
% PRPico = find(abs(VPico) == max(abs(VPico)));
% PRPico = PRPico-6;
% PicoEnc = abs(ECG(picos_R(Bat)+PRPico(1))-m0);

PicoEnc = abs(ECG(picos_R(Bat))-m0);

P_PEst = 1;
P_PEnc = 0.8;
limiar = (P_PEst*PicoEst + P_PEnc*PicoEnc)*FatAmp/(P_PEst + P_PEnc);




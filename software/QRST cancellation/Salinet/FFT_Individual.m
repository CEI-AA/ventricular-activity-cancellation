TABLE = QRST_3_10s;
fs=1200; % noncontact
factor=2.8571; 
sizefft=factor*8400;
fstep=fs/sizefft;
H=hamming (sizefft/factor);
freq_down=0.01;
freq_up=25;
position=1;



    Sfft=abs(fft(TABLE(:, position).*H, sizefft)); % Zero Padding
    Sfft=Sfft.*Sfft/length(Sfft);  % for power 
    Sffti (:,position)=Sfft(1:floor(freq_up*(1/fstep))); 
%     [P,F] = max(Sfft(floor(freq_down*(1/fstep)):floor(freq_up*(1/fstep))));  
%     Fsam=(F+floor(freq_down*(1/fstep))-1);  
%     Fhz=Fsam*fstep;
%         if  (Fhz>=freq_down&&Fhz<=freq_up)       
%             MFFTi(position, 1)=Fhz;
%         else
%             MFFTi(position, 1)=0;
%         end %if
%         
%    
% plot ((1:length(Sfft(1:freq_up/fstep)))*(fstep), Sfft(1:freq_up/fstep)); 
% 
%  

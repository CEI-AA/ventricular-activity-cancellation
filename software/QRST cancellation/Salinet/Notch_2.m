% clear;
% load CC_CS1_10_preatro_notched_50Hz.mat
f= 50; [50:50:1000] ; % frequency notch
steps = length(f);
fs=1000; %fs=2034.50;
alpha = 0.96; %0.95; %0.985; % 0.993
data = data_mV;

%for nn=1:steps
    
 omega = 2*pi*f/fs; % fs = sampling frequency
 num = [1 -2*cos(omega) 1];
 den = [1 -2*alpha*cos(omega) alpha*alpha];
    
    
    for n=1:5
        n   
        TABLE = data(: ,n);    
        Notched (:,n) = filtfilt(num,den,TABLE);    %Zero-phase forward and reverse digital filtering
    end
    
    %end

clear omega num n fs f den alpha TABLE data
data = Notched; 
clear Notched; 

plot (1:length(data(:,1)), data (:,5),'r',1:length(data_mV(:,1)), data_mV (:,5), 'k')
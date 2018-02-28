% clear;
% load CC_CS1_10_preatro_notched_50Hz.mat
f=11.25;                                % frequency notch
fs=1200;
%fs=2034.50;
alpha = 0.95; %0.985; % 0.993
omega = 2*pi*f/fs; % fs = sampling frequency
num = [1 -2*cos(omega) 1];
den = [1 -2*alpha*cos(omega) alpha*alpha];

for n=1:2048
 n   
TABLE = data(: ,n);    
Notched (:,n) = filtfilt(num,den,TABLE);    %Zero-phase forward and reverse digital filtering
   
end


clear omega num n fs f den alpha TABLE data

data = Notched; 

clear Notched; 

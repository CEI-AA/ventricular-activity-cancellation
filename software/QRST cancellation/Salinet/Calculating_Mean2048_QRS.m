     
%% Average Spectrum
%% Joao Loures Salinet Jr, 19/04/2012
%% University of Leicester - Engineering Department

%% Method 1
load Pat8_FFT_Raw_20s.mat;
% Area = S;
% 
% for n=1:length (Area)
%     
% Sfft11 (:,n) = Sfft1 (:, Area(n));
% Sfft21 (:,n) = Sfft2 (:, Area(n));
% Sfft31 (:,n) = Sfft3 (:, Area(n));
% Sfft41 (:,n) = Sfft4 (:, Area(n));
% Sfft51 (:,n) = Sfft5 (:, Area(n));
% Sfft61 (:,n) = Sfft6 (:, Area(n));
% Sfft71 (:,n) = Sfft7 (:, Area(n));
% Sfft81 (:,n) = Sfft8 (:, Area(n));
% Sfft91 (:,n) = Sfft9 (:, Area(n));
%     
% end
% 
% clear n Area;
% clear Sfft1 Sfft2 Sfft3 Sfft4 Sfft5 Sfft6 Sfft7 Sfft8 Sfft9;
% 
% 
% Sfft1 = Sfft11;
% Sfft2 = Sfft21;
% Sfft3 = Sfft31;
% Sfft4 = Sfft41;
% Sfft5 = Sfft51;
% Sfft6 = Sfft61;
% Sfft7 = Sfft71;
% Sfft8 = Sfft81;
% Sfft9 = Sfft91;
% 
% clear Sfft11 Sfft21 Sfft31 Sfft41 Sfft51 Sfft61 Sfft71 Sfft81 Sfft91;

Sfft1 = Sfft1';
Sfft2 = Sfft2';
Sfft3 = Sfft3';
Sfft4 = Sfft4';
Sfft5 = Sfft5';
Sfft6 = Sfft6';
Sfft7 = Sfft7';
Sfft8 = Sfft8';
Sfft9 = Sfft9';


Mean_SFFT1 = mean(Sfft1);
Mean_SFFT2 = mean(Sfft2);
Mean_SFFT3 = mean(Sfft3);
Mean_SFFT4 = mean(Sfft4);
Mean_SFFT5 = mean(Sfft5);
Mean_SFFT6 = mean(Sfft6);
Mean_SFFT7 = mean(Sfft7);
Mean_SFFT8 = mean(Sfft8);
Mean_SFFT9 = mean(Sfft9);



Mean_SFFTall=(Mean_SFFT1+Mean_SFFT2+Mean_SFFT3+Mean_SFFT4+Mean_SFFT5+Mean_SFFT6+Mean_SFFT7+Mean_SFFT8+Mean_SFFT9)/9;
        
Mean_SFFTall = Mean_SFFTall';
%Mean_SFFTall_2048_QRST = Mean_SFFTall';

%clear Mean_SFFTall;
clear Mean_SFFT1 Mean_SFFT2 Mean_SFFT3 Mean_SFFT4 Mean_SFFT5 Mean_SFFT6 Mean_SFFT7 Mean_SFFT8 Mean_SFFT9;
clear Sfft1 Sfft2 Sfft3 Sfft4 Sfft5 Sfft6 Sfft7 Sfft8 Sfft9;

%% Method 2


%SFFT_mean=(Sfft1+Sfft2+Sfft3+Sfft4+Sfft5+Sfft6+Sfft7+Sfft8+Sfft9)/9;

% %% Phase 3
% 
%  %694 1338 1454	442	1776 731 1180 1439
% 
% 
% SFFT_mean_Raw = SFFT_mean_Raw';
% SFFT_mean_Anita = SFFT_mean_Anita';
% SFFT_mean_averageABS = SFFT_mean_averageABS';
% SFFT_mean_medianABS = SFFT_mean_medianABS';
% SFFT_mean_averageQRST = SFFT_mean_averageQRST';
% SFFT_mean_medianQRST = SFFT_mean_medianQRST';
% 
% P1_SFFT_mean_Raw = SFFT_mean_Raw(:,694);
% P1_SFFT_mean_Anita = SFFT_mean_Anita (:,694);
% P1_SFFT_mean_averageABS = SFFT_mean_averageABS (:,694);
% P1_SFFT_mean_medianABS = SFFT_mean_medianABS (:,694);
% P1_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,694);
% P1_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,694);
% 
% P2_SFFT_mean_Raw = SFFT_mean_Raw(:,1338);
% P2_SFFT_mean_Anita = SFFT_mean_Anita (:,1338);
% P2_SFFT_mean_averageABS = SFFT_mean_averageABS (:,1338);
% P2_SFFT_mean_medianABS = SFFT_mean_medianABS (:,1338);
% P2_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,1338);
% P2_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,1338);
% 
% P3_SFFT_mean_Raw = SFFT_mean_Raw(:,1454);
% P3_SFFT_mean_Anita = SFFT_mean_Anita (:,1454);
% P3_SFFT_mean_averageABS = SFFT_mean_averageABS (:,1454);
% P3_SFFT_mean_medianABS = SFFT_mean_medianABS (:,1454);
% P3_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,1454);
% P3_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,1454);
% 
% P4_SFFT_mean_Raw = SFFT_mean_Raw(:,442);
% P4_SFFT_mean_Anita = SFFT_mean_Anita (:,442);
% P4_SFFT_mean_averageABS = SFFT_mean_averageABS (:,442);
% P4_SFFT_mean_medianABS = SFFT_mean_medianABS (:,442);
% P4_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,442);
% P4_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,442);
% 
% P5_SFFT_mean_Raw = SFFT_mean_Raw(:,1776);
% P5_SFFT_mean_Anita = SFFT_mean_Anita (:,1776);
% P5_SFFT_mean_averageABS = SFFT_mean_averageABS (:,1776);
% P5_SFFT_mean_medianABS = SFFT_mean_medianABS (:,1776);
% P5_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,1776);
% P5_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,1776);
% 
% P6_SFFT_mean_Raw = SFFT_mean_Raw(:,731);
% P6_SFFT_mean_Anita = SFFT_mean_Anita (:,731);
% P6_SFFT_mean_averageABS = SFFT_mean_averageABS (:,731);
% P6_SFFT_mean_medianABS = SFFT_mean_medianABS (:,731);
% P6_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,731);
% P6_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,731);
% 
% P7_SFFT_mean_Raw = SFFT_mean_Raw(:,1180);
% P7_SFFT_mean_Anita = SFFT_mean_Anita (:,1180);
% P7_SFFT_mean_averageABS = SFFT_mean_averageABS (:,1180);
% P7_SFFT_mean_medianABS = SFFT_mean_medianABS (:,1180);
% P7_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,1180);
% P7_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,1180);
% 
% P8_SFFT_mean_Raw = SFFT_mean_Raw(:,1439);
% P8_SFFT_mean_Anita = SFFT_mean_Anita (:,1439);
% P8_SFFT_mean_averageABS = SFFT_mean_averageABS (:,1439);
% P8_SFFT_mean_medianABS = SFFT_mean_medianABS (:,1439);
% P8_SFFT_mean_averageQRST = SFFT_mean_averageQRST (:,1439);
% P8_SFFT_mean_medianQRST = SFFT_mean_medianQRST (:,1439);

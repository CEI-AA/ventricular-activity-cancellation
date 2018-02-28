
fs= 2034.5; %1200; 
t = 4;
factor=5; 
sizefft=factor*t*fs;
fstep=fs/sizefft;
freq_up=20;

%% Anita_AW; ABS_AW; QRST_AW;
figure;
plot ((1:length(Raw_AW(1:freq_up/fstep)))*(fstep), Raw_AW(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_AW(1:freq_up/fstep)))*(fstep), Anita_AW(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_AW(1:freq_up/fstep)))*(fstep), ABS_AW(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_AW(1:freq_up/fstep)))*(fstep), QRST_AW(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Anterior Wall', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_IW; ABS_IW; QRST_IW
figure;
plot ((1:length(Raw_IW(1:freq_up/fstep)))*(fstep), Raw_IW(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_IW(1:freq_up/fstep)))*(fstep), Anita_IW(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_IW(1:freq_up/fstep)))*(fstep), ABS_IW(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_IW(1:freq_up/fstep)))*(fstep), QRST_IW(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Inferior Wall', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_LW; ABS_LW; QRST_LW
figure;
plot ((1:length(Raw_LW(1:freq_up/fstep)))*(fstep), Raw_LW(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_LW(1:freq_up/fstep)))*(fstep), Anita_LW(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_LW(1:freq_up/fstep)))*(fstep), ABS_LW(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_LW(1:freq_up/fstep)))*(fstep), QRST_LW(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Lateral Wall', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_MA; ABS_MA; QRST_MA;
figure;
plot ((1:length(Raw_MA(1:freq_up/fstep)))*(fstep), Raw_MA(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_MA(1:freq_up/fstep)))*(fstep), Anita_MA(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_MA(1:freq_up/fstep)))*(fstep), ABS_MA(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_MA(1:freq_up/fstep)))*(fstep), QRST_MA(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Mitral Annulus', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_PV; ABS_PV; QRST_PV;
figure;
plot ((1:length(Raw_PV(1:freq_up/fstep)))*(fstep), Raw_PV(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_PV(1:freq_up/fstep)))*(fstep), Anita_PV(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_PV(1:freq_up/fstep)))*(fstep), ABS_PV(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_PV(1:freq_up/fstep)))*(fstep), QRST_PV(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Close to the Pulmonary Veins', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_PW; ABS_PW; QRST_PW;
figure;
plot ((1:length(Raw_PW(1:freq_up/fstep)))*(fstep), Raw_PW(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_PW(1:freq_up/fstep)))*(fstep), Anita_PW(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_PW(1:freq_up/fstep)))*(fstep), ABS_PW(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_PW(1:freq_up/fstep)))*(fstep), QRST_PW(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Posterior Wall', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_R; ABS_R; QRST_R;
figure;
plot ((1:length(Raw_R(1:freq_up/fstep)))*(fstep), Raw_R(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_R(1:freq_up/fstep)))*(fstep), Anita_R(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_R(1:freq_up/fstep)))*(fstep), ABS_R(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_R(1:freq_up/fstep)))*(fstep), QRST_R(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Roof', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Anita_S; ABS_S; QRST_S;
figure;
plot ((1:length(Raw_S(1:freq_up/fstep)))*(fstep), Raw_S(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Anita_S(1:freq_up/fstep)))*(fstep), Anita_S(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(ABS_S(1:freq_up/fstep)))*(fstep), ABS_S(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(QRST_S(1:freq_up/fstep)))*(fstep), QRST_S(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 7 - Septum', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%% Mean_SFFTall_2048_Anita; Mean_SFFTall_2048_ABS; Mean_SFFTall_2048_QRST;
figure;
plot ((1:length(Mean_SFFTall_2048_Raw(1:freq_up/fstep)))*(fstep), Mean_SFFTall_2048_Raw(1:freq_up/fstep), 'g', 'LineWidth',2); 
hold on;
plot ((1:length(Mean_SFFTall_2048_Anita(1:freq_up/fstep)))*(fstep), Mean_SFFTall_2048_Anita(1:freq_up/fstep), 'k', 'LineWidth',2); 
hold on;
plot ((1:length(Mean_SFFTall_2048_ABS(1:freq_up/fstep)))*(fstep), Mean_SFFTall_2048_ABS(1:freq_up/fstep), 'b', 'LineWidth',2);
hold on;
plot ((1:length(Mean_SFFTall_2048_QRST(1:freq_up/fstep)))*(fstep), Mean_SFFTall_2048_QRST(1:freq_up/fstep), 'r', 'LineWidth',2);

xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
title ('Patient 8 - Whole Atrium', 'Color','k', 'FontSize',14)
legend('Raw', 'QRS','QRST-ABS','QRST-Segmentation')
set(legend, 'Box', 'off')

%%

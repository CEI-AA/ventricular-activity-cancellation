
fs= 2034.5; % 1200; 
t = 4;
factor=5; 
sizefft=factor*t*fs;
fstep=fs/sizefft;
freq_down_T=3;
freq_up_T=5.5;
freq_down=3;
freq_up=12;




%load Pat_all_2048_average.mat;

%Matrix = cat (2, Raw_AW, Raw_IW, Raw_LW, Raw_MA, Raw_PV, Raw_PW, Raw_R, Raw_S, Mean_SFFTall_2048_Raw); 
% Matrix = cat (2, Anita_AW, Anita_IW, Anita_LW, Anita_MA, Anita_PV, Anita_PW, Anita_R, Anita_S, Mean_SFFTall_2048_Anita); 
 %Matrix = cat (2, ABS_AW, ABS_IW, ABS_LW, ABS_MA, ABS_PV, ABS_PW, ABS_R, ABS_S, Mean_SFFTall_2048_ABS); 
 %Matrix = cat (2, QRST_AW, QRST_IW, QRST_LW, QRST_MA, QRST_PV, QRST_PW, QRST_R, QRST_S, Mean_SFFTall_2048_QRST); 

% Matrix = cat (2, Pat1_Raw, Pat2_Raw, Pat3_Raw, Pat4_Raw, Pat5_Raw, Pat6_Raw, Pat7_Raw); 
% Matrix = cat (2, Pat1_Anita, Pat2_Anita, Pat3_Anita, Pat4_Anita, Pat5_Anita, Pat6_Anita, Pat7_Anita); 
% Matrix = cat (2, Pat1_ABS, Pat2_ABS, Pat3_ABS, Pat4_ABS, Pat5_ABS, Pat6_ABS, Pat7_ABS); 
% Matrix = cat (2, Pat1_QRST, Pat2_QRST, Pat3_QRST, Pat4_QRST, Pat5_QRST, Pat6_QRST, Pat7_QRST); 
 
% Matrix = Pat8_Raw; 
% Matrix = Pat8_ABS; 
  Matrix = Mean_SFFTall_Anita; 


 
for var = 1%:7;
    
Patall_SNR(1,var)= sum(Matrix(floor(freq_down_T*(1/fstep)):floor(freq_up_T*(1/fstep)),var));  
Patall_SNR (2, var) = sum(Matrix(floor(freq_down*(1/fstep)):floor(freq_up*(1/fstep)), var)); 
Patall_SNR (3,var) = Patall_SNR(1, var)/Patall_SNR (2, var);

end

%clear R*; clear A*; clear Q*; clear M*;

%Pat7_SNR_All = cat (1, Pat7_SNR_Raw (3,:), Pat7_SNR_Anita (3,:), Pat7_SNR_ABS (3,:), Pat7_SNR_QRST (3,:));    


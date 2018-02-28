%% Spectral Analysis

MAT_filter = Anita;  %AEG; %  

FS=2034.5; % 1200; 
sizefft = FS*4; % noncontact
%sizefft = 4800; % noncontact
l = length(MAT_filter(:,1));
Int_part = floor(l/sizefft);
Rem_no_of_lines = l-(Int_part*sizefft);
Start_point = 0; End_point = 0;
j=0;

for i = 1: 40 %(2*Int_part-1);
    if (i==1)
    Start_point = 1 + End_point;
    End_point = End_point + sizefft;
    else
        Start_point = Start_point + sizefft/2;
        End_point = End_point + (sizefft/2);
    end
    A = MAT_filter(Start_point:End_point, :);
    j=j+1;
    eval(sprintf('[MFFT%d, Sfft%d] = Perform_FFT_4s_max_4_20Hz_indices(A);', [j j]));
    

end

%Copy the remaining part in T
Start_point = Start_point + sizefft/2;
End_point = End_point+Rem_no_of_lines;
T=MAT_filter(Start_point:End_point, :);

 clear End_point
 clear Int_part Rem_no_of_lines Start_point i l sizefft A Wn2 f3 f4 j End_point;
 clear MAT_filter

 %MFFT_all= cat (1, MFFT1, MFFT2, MFFT3, MFFT4, MFFT5, MFFT6, MFFT7, MFFT8, MFFT9);
 %Sfft_all= cat (1, Sfft1', Sfft2', Sfft3', Sfft4', Sfft5', Sfft6', Sfft7', Sfft8', Sfft9');

 

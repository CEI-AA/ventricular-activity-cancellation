% This code reads the raw data files (each files about 7 s) from the External ECGs and Ablations
% contact recordings and concatenating them in just one file.
% Please, have a look in the raw files ".txt" to identify what each column is.
% If you want more than 5 files, just copy and past and change 
   % 1) fid - with next name of the file
   % 2) variable yw_*, which usually correspond for each file (yw_1-file 1,
   % yw_2-file 2 and so on...)
   
   %PS: Freq Sampling = 1200 Hz
   
   % By Joao Salinet 28/02/12 - University of Leicester


% ********  First File - Wave ******************************************
tic
fid = fopen('Jaya_post_seg1_ecg_Wav_20140523_131220_0001.txt');
if fid==-1
  error('File not found or permission denied');
end	%if

beginwavData = 0;
while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
    mywavline = fgetl(fid);
    beginwavData = beginwavData + 1;
    % Exit the loop when the keyword is found.
     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%     if strfind(mywavline, 'Begin data') % compare Begin data with the line
       break;        
    end %if
      
end %while

mywavdata=[];
mywavdata = [mywavdata; fscanf(fid, '%f,')];
fclose(fid);

D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[2048,length(mywavdata)/2048]);
  D=reshape(mywavdata,[10,length(mywavdata)/10]);

    whos D;
    D=D';   % transpose it
    yw_1=D(1:length(D),:);
    
    clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans;  
    
    Reading_file_done=toc
    
 % ********  Second File - Wave ******************************************
 
 tic
fid = fopen('Jaya_post_seg2_ecg_Wav_20140523_132535_0001.txt');
if fid==-1
  error('File not found or permission denied');
end	%if

beginwavData = 0;
while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
    mywavline = fgetl(fid);
    beginwavData = beginwavData + 1;
    % Exit the loop when the keyword is found.
    if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
       break;        
    end %if
      
end %while

mywavdata=[];
mywavdata = [mywavdata; fscanf(fid, '%f,')];
fclose(fid);

D=mywavdata; % Create Matrix C = B       
    D=reshape(mywavdata,[10,length(mywavdata)/10]);
    whos D;
    D=D';   % transpose it
    yw_2=D(1:length(D),:);
   
    clear D; clear mywavdata; clear mywavline; clear beginwavData; 
   
  
    
Reading_file_done=toc     

    % ********  Third File - Wave ******************************************
    
     tic
fid = fopen('Jaya_post_seg3_ecg_Wav_20140523_134021_0001.txt');
if fid==-1
  error('File not found or permission denied');
end	%if

beginwavData = 0;
while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
    mywavline = fgetl(fid);
    beginwavData = beginwavData + 1;
    % Exit the loop when the keyword is found.
    if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
       break;        
    end %if
      
end %while

mywavdata=[];
mywavdata = [mywavdata; fscanf(fid, '%f,')];
fclose(fid);

D=mywavdata; % Create Matrix C = B       
    D=reshape(mywavdata,[10,length(mywavdata)/10]);
    whos D;
    D=D';   % transpose it
    yw_3=D(1:length(D),:);
     
    clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
  
   Reading_file_done=toc
   
     % ********  Fourth File - Wave ******************************************
%     
%      tic
% fid = fopen('base_wav_4.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_4=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%      % ********  Fifth File - Wave ******************************************
%     
% tic
% fid = fopen('base_wav_5.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_5=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
% 
%  % ********  Sixth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_6.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_6=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
%    
%    
% % ********  Seventh File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_7.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_7=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  eighth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_8.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_8=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  ninth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_9.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_9=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% 
% % ********  tenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_10.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_10=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
%    
%    
% % ********  eleventh File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_11.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_11=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  twelth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_12.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_12=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  thirteenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_13.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_13=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
%    
% % ********  fourteenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_14.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_14=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  fifteenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_15.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_15=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  sixteenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_16.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_16=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
% 
% % ********  seventeenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_17.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_17=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
%    
%    
% 
% % ********  eighteenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_18.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_18=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
%    
%    
% % ********  nineteenth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_19.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_19=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
% 
%    % ********  twentieth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_20.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_20=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
% % ********  twenty first File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_21.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_21=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% 
% % ********  twenty second File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_22.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_22=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
%    
% % ********  twenty third File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_23.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_23=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% % ********  twenty fourth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_24.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_24=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
%    
%    
% % ********  twenty fifth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_25.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_25=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
% 
% % ********  twenty sixth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_26.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_26=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
%    
%    
% 
% % ********  twenty seventh File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_27.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_27=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
%    
% % ********  twenty eighth File - Wave ******************************************
%    
%    
% tic
% fid = fopen('base_wav_28.txt');
% if fid==-1
%   error('File not found or permission denied');
% end	%if
% 
% beginwavData = 0;
% while ~feof(fid)% Test for end-of-file (ignore line blank and '*'
%     mywavline = fgetl(fid);
%     beginwavData = beginwavData + 1;
%     % Exit the loop when the keyword is found.
%     if strfind(mywavline, 'Wave  0, Wave  1, Wave  2, Wave  3, Wave  4, Wave  5, Wave  6, Wave  7, Wave  8, Wave  9,') % compare Begin data with the line
%        break;        
%     end %if
%       
% end %while
% 
% mywavdata=[];
% mywavdata = [mywavdata; fscanf(fid, '%f,')];
% fclose(fid);
% 
% D=mywavdata; % Create Matrix C = B       
%     D=reshape(mywavdata,[12,length(mywavdata)/12]);
%     whos D;
%     D=D';   % transpose it
%     yw_28=D(1:length(D),:);
%      
%     clear D; clear mywavdata; clear mywavline; clear beginwavData; clear ans; 
%   
%    Reading_file_done=toc
% 
% 
% 
%    
%  %% Concatenating all 28 files (which correspond for about 28*7 s - 196 sec)
%    
%    
%    Final_yw=cat(1,yw_1,yw_2,yw_3,yw_4,yw_5,yw_6,yw_7,yw_8,yw_9,yw_10,yw_11,yw_12,yw_13,yw_14,yw_15,yw_16,yw_17,yw_18,yw_19,yw_20,yw_21,yw_22,yw_23,yw_24,yw_25,yw_26,yw_27,yw_28); clear yw*; clear Reading_file_done; clear fid; whos;
%     
%     
  

    
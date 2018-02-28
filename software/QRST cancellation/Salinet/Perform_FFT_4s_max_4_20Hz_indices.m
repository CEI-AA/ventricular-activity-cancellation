function [MFFTi,Sffti] = Perform_FFT_4s_max_4_20Hz_indices(TABLE)

% fs=1000; %contact
fs= 2034.5; %1200; % noncontact  %
factor=5; 
sizefft=factor*fs*4;
fstep=fs/sizefft;
H=hamming (sizefft/factor);
freq_down=4;
freq_up=20;
freq_up_2=20;
down_index=4/fstep;



for position=1:2048
    Sfft=abs(fft(TABLE(:, position).*H, sizefft)); % Zero Padding
    Sfft=Sfft.*Sfft/length(Sfft);  % for power 
    Sffti (:,position)=Sfft(1:floor(freq_up_2*(1/fstep))); 
    [P,F] = max(Sfft(floor(freq_down*(1/fstep)):floor(freq_up_2*(1/fstep))));  
    Fsam=(F+floor(freq_down*(1/fstep))-1);  
    Fhz=Fsam*fstep;
        if  (Fhz>=freq_down&&Fhz<=freq_up)       
            MFFTi(position, 1)=Fhz;
        else
            MFFTi(position, 1)=0;
        end %if
        
         %% Regularity Index
  
     Area_Peak_RI08_FFT(position) = sum(Sfft(Fsam-8:Fsam+8)); 
     Area_Band_RI_FFT(position) = sum(Sfft(floor(freq_down*(1/fstep)):floor(freq_up*(1/fstep))));
     MFFTi (position, 2)= Area_Peak_RI08_FFT(position)/ Area_Band_RI_FFT(position);
      
     
     %% Organization Index
     
   
       if  (Fhz<=4.89) 
           Harmonics = [Fsam 2*Fsam 3*Fsam 4*Fsam];
              for Pos_Harm = 1: length (Harmonics)
                  for samples_FFT= Harmonics(Pos_Harm)+1:450 
                      if (Sfft(samples_FFT)<Sfft(samples_FFT-1)) && (Sfft(samples_FFT)<Sfft(samples_FFT+1)) 
                          Fsam_up (Pos_Harm) = samples_FFT; 
                      break
                      end
                  end     

                  clear samples_FFT;

                   for samples_FFT= Harmonics(Pos_Harm)-1:-1:60
                       if (Sfft(samples_FFT) < Sfft(samples_FFT+1)) && (Sfft(samples_FFT)<Sfft(samples_FFT-1))
                           Fsam_down (Pos_Harm) = samples_FFT; 
                       break
                       end
                   end 
              end 
             
         OI_FFT(position) = sum(Sfft(Fsam_down(1):Fsam_up(1))) + sum(Sfft(Fsam_down(2):Fsam_up(2))) + sum(Sfft(Fsam_down(3):Fsam_up(3))) + sum(Sfft(Fsam_down(4):Fsam_up(4))); 
         Band_OI_FFT(position) = sum(Sfft(down_index:400));
         MFFTi (position, 3)= OI_FFT(position)/ Band_OI_FFT(position);    
       end
       
       
        if   (Fhz>=4.9&&Fhz<=6.49)   
           Harmonics = [Fsam 2*Fsam 3*Fsam];
              for Pos_Harm = 1: length (Harmonics)
                  for samples_FFT= Harmonics(Pos_Harm)+1:450 
                      if (Sfft(samples_FFT)<Sfft(samples_FFT-1)) && (Sfft(samples_FFT)<Sfft(samples_FFT+1)) 
                          Fsam_up (Pos_Harm) = samples_FFT; 
                      break
                      end
                  end     

                  clear samples_FFT;

                   for samples_FFT= Harmonics(Pos_Harm)-1:-1:60
                       if (Sfft(samples_FFT) < Sfft(samples_FFT+1)) && (Sfft(samples_FFT)<Sfft(samples_FFT-1))
                           Fsam_down (Pos_Harm) = samples_FFT; 
                       break
                       end
                   end 
              end 
             
         OI_FFT(position) = sum(Sfft(Fsam_down(1):Fsam_up(1))) + sum(Sfft(Fsam_down(2):Fsam_up(2))) + sum(Sfft(Fsam_down(3):Fsam_up(3))); 
         Band_OI_FFT(position) = sum(Sfft(down_index:400));
         MFFTi (position, 3)= OI_FFT(position)/ Band_OI_FFT(position);    
        end
  
       
         if   (Fhz>=6.5&&Fhz<=9.79)   
           Harmonics = [Fsam 2*Fsam];
              for Pos_Harm = 1: length (Harmonics)
                  for samples_FFT= Harmonics(Pos_Harm)+1:450 
                      if (Sfft(samples_FFT)<Sfft(samples_FFT-1)) && (Sfft(samples_FFT)<Sfft(samples_FFT+1)) 
                          Fsam_up (Pos_Harm) = samples_FFT; 
                      break
                      end
                  end     

                  clear samples_FFT;

                   for samples_FFT= Harmonics(Pos_Harm)-1:-1:70
                       if (Sfft(samples_FFT) < Sfft(samples_FFT+1)) && (Sfft(samples_FFT)<Sfft(samples_FFT-1))
                           Fsam_down (Pos_Harm) = samples_FFT; 
                       break
                       end
                   end 
              end 
             
         OI_FFT(position) = sum(Sfft(Fsam_down(1):Fsam_up(1))) + sum(Sfft(Fsam_down(2):Fsam_up(2))); 
         Band_OI_FFT(position) = sum(Sfft(down_index:400));
         MFFTi (position, 3)= OI_FFT(position)/ Band_OI_FFT(position);    
         end  
       
 
        if   (Fhz>=9.8&&Fhz<=19.3)   
           Harmonics = [Fsam];
              for Pos_Harm = 1: length (Harmonics)
                  for samples_FFT= Harmonics(Pos_Harm)+1:450 
                      if (Sfft(samples_FFT)<Sfft(samples_FFT-1)) && (Sfft(samples_FFT)<Sfft(samples_FFT+1)) 
                          Fsam_up (Pos_Harm) = samples_FFT; 
                      break
                      end
                  end     

                  clear samples_FFT;

                   for samples_FFT= Harmonics(Pos_Harm)-1:-1:70
                       if (Sfft(samples_FFT) < Sfft(samples_FFT+1)) && (Sfft(samples_FFT)<Sfft(samples_FFT-1))
                           Fsam_down (Pos_Harm) = samples_FFT; 
                       break
                       end
                   end 
              end 
             
         OI_FFT(position) = sum(Sfft(Fsam_down(1):Fsam_up(1)));  
         Band_OI_FFT(position) = sum(Sfft(down_index:400));
         MFFTi (position, 3)= OI_FFT(position)/ Band_OI_FFT(position);    
       end 
 end
 

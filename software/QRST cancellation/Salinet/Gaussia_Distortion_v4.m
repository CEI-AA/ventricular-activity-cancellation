function[Gsyn] = Gaussia_Distortion_v4(SD,alfa)

global sfreq
N=0.500*sfreq;  %% Resolution

%function[Y,X]=GAUSSIAN_DISTORTION(N,SD,alfa)
%% N is the gaussian resolution (number of samples)
%% SD is the original standard deviation
%% alfa is the angle related to the first degree of distortion
Gsyn = 0;
a = tan(alfa); %% The angle declivity of the function



    %% Generating original gaussian function
    
limmax=3;
limmin=-3;
size=limmax-limmin;
passo=size/N;
X=limmin:passo:limmax;
Y=(1/(sqrt(2*pi)*SD)).*exp((-1.*X.^2)/(2*SD^2));
Y = 10*(Y/max(Y)); % Gaussian amplitude

a = tan(alfa);
i = 1;
Y2=0;
X2=[];

while i<=length(X) %% For every original sample of the gaussian function
    at = Y(i)*(a+eps);  %% Delay value for building distorted gaussian function
        %if at >= 1 %% The delay starts from a given sample
            %if i - round(at) >= 1 %% Protecting to not have negative indices
    Y2(i) = Y(i);
    X2(i) = X(i)+at;
            %end
        %else Y2(i) = Y(i); %% For low values of the original gaussian amplitudes
            %X2 = [X2 X(i)];
        %end
    i=i+1;
end
    
espX = X2(2:end)-X2(1:end-1);
tval = find(espX < 0);
if length(tval) == 0
    
    X2i = X;
    Y2i = INTERP1(X2,Y2,X2i,'cubic');
    X2 = X2i;
    Y2 = Y2i;
    Gsyn = Y2;
    
    

end
    


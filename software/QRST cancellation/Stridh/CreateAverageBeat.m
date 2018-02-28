function [ X ] = CreateAverageBeat( Y, Fs, R, del, N )
%CreateAverageBeat creates an average beat based on a signal
%Inputs:
%Y: Signal input, Nbeats samples by L leads
%Fs: Sampling Frequency
%R: Vector of R-peak indices
%del: correction factor
%N: Number of Samples in a beat
%Outputs:
%X: Average Beat 
[~, L] = size(Y);
[~, Nbeats] = size(R);
X = cell(1,L);
s = zeros(Nbeats,1);
e = zeros(Nbeats,1);
RRmin=intmax;
for i = 2:Nbeats
    Rdif = R(i)-R(i-1);
    RRmin = min(RRmin,Rdif); %minimum R-R difference
end

for i = 1:Nbeats
    s(i) = floor(R(i)-0.3*RRmin); %vector of window open indices
    e(i) = floor(R(i)+0.7*RRmin); %vector of window closing indices
end
% for i = 1:Nbeats
%     s(i) = floor(R(i)-0.04*Fs); %vector of window open indices
%     e(i) = floor(R(i)+RRmin-0.04*Fs-1); %vector of window closing indices
% end
if s(1)<=0
    R(1)=[];
    s(1)=[];
    e(1)=[];
    Nbeats=Nbeats-1;
end
if e(Nbeats)>size(Y,1)
    Nbeats=Nbeats-1;
end
%N = e(1)-s(1); %N = number of Samples
for i = 1:L
    x=zeros(N+2*del,1);
    for j = 1:Nbeats
        for k = 1:N
            Z(j,k) = Y(s(j)+k-1,i); %each row of Z is an individual beat
        end
    end
    for j = 1:N
        x(j+del) = mean(Z(:,j)); %averages the entire column
    end
    X{i} = x;
end
X=cell2mat(X);
end


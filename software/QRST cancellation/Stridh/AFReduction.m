function [ YAtilde ] = AFReduction( Y, Fs, R )
%AFReduction calculates the intermediate AF reduction step for the
%   StridhQRST method
%Inputs:
%Y: ECG input
%Fs: Sampling Frequency
%R: Vector of R-peak indices
%
%Outputs:
%YAtilde: AF reduced wave
[~,L] = size(Y); %L = number of leads
RRmin = intmax;
[~, Nbeats] = size(R); %Nbeats = number of beats
s = zeros(Nbeats,1); 
e = zeros(Nbeats,1);
YAtilde=zeros(size(Y,1),L);
for i = 2:Nbeats
    Rdif = R(i)-R(i-1);
    RRmin = min(RRmin,Rdif); %minimum R-R difference
end
for i = 1:Nbeats
    s(i) = floor(R(i)-0.3*RRmin); %vector of window open indices
    e(i) = floor(R(i)+0.7*RRmin); %vector of window closing indices
end
if s(1)<=0
    R(1)=[];
    s(1)=[];
    e(1)=[];
    Nbeats=Nbeats-1;
end
if(e(Nbeats)>size(Y,1))
    Nbeats=Nbeats-1;
end
Q=zeros(1,Nbeats);
T=zeros(1,Nbeats);
for i = 1:Nbeats
    Q(i)=floor(R(i)-0.1*Fs); %Q-onset Point for each beat
    T(i)=e(i); %T-endpoint for each beat
end
TQintervals=cell(L,Nbeats);
TQintervalsIndices=cell(1,Nbeats);
for i=1:L
    for j=1:Nbeats-1
        TQintervals{i,j}=Y((T(j):Q(j+1)),i);%collection of all TQintervals
    end
end
for i=1:Nbeats-1
    TQintervalsIndices{i}=T(i):Q(i+1);
end
for i=1:L
    for j=1:Nbeats
        n=T(j)-Q(j); %size of QRS complex
        w=linspace(0,1,n); %weighting vector
        if(j==1)
            b=size(TQintervals{i,j},1);
            if (b<n)
                xq=linspace(T(j),Q(j+1),n); %create n values in TQ-interval
                vq=interp1(TQintervalsIndices{j},TQintervals{i,j},xq); %values 
            else
                vq=Y(T(j):T(j)+n); %values to the right
            end
            wvq=zeros(1,n);
            for k=1:n
                wvq(k)=vq(k)*-w(k);
            end
        elseif (j==Nbeats)
            a=size(TQintervals{i,j-1},1);
            if (a<n)
                xq=linspace(T(j-1),Q(j),n); %create n values in TQ-interval
                vq=interp1(TQintervalsIndices{j-1},TQintervals{i,j-1},xq); %values 
            else
                vq=Y(Q(j)-n:Q(j)); %values to the left
            end
            wvq=zeros(1,n);
            for k=1:n
                wvq(k)=vq(k)*w(k);
            end
        else
            a=size(TQintervals{i,j-1});
            b=size(TQintervals{i,j},1);
            if (a<n)
                xql=linspace(T(j-1),Q(j),n); %create n values in TQ-interval
                vql=interp1(TQintervalsIndices{j-1},TQintervals{i,j-1},xql); %values 
            else
                vql=Y(Q(j)-n:Q(j)); %values to left
            end
            if (b<n)
                xqr=linspace(T(j),Q(j+1),n); %create n values in TQ-interval
                vqr=interp1(TQintervalsIndices{j},TQintervals{i,j},xqr); %values 
            else
                vqr=Y(T(j):T(j)+n-1); %values to right
            end
            wvq=zeros(1,n);
            for k=1:n
                wvq(k)=vql(k)*w(k)+vqr(k)*-w(k);
            end
        end
        for k=1:n
            YAtilde(Q(j)+k-1,i)=wvq(k);
        end
    end
end
end


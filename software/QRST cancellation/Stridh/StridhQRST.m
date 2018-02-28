function [ Ya ] = StridhQRST( Y, R, del, MaxIt, Fs, eps )
%StridhQRST is a function to calculate and cancel the QRS complex of an ECG
%   to leave the Atrial Fibrillation data behind
%inputs:
%Y: NxL matrix of values corresponding to the ecg, N is number of samples,
%   L is number of leads
%R: vector of R-peak locations;
%del: maximum synchronisation in time
%MaxIt: maximum number of iterations
%Fs: Sampling Frequency
%eps: maximum error value
Ya=Y;
[~, L] = size(Y); %L = number of leads
if(L>12)
    L = 12;
end
del = del;
[~,Nbeats] = size(R); %Nbeats = number of beats
s = zeros(Nbeats,1);
e = zeros(Nbeats,1);
RRmin = intmax;
% RRmin = min(diff(R))
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
if (e(Nbeats)>size(Y,1))
    Nbeats=Nbeats-1;
end
N = e(1)-s(1); %N = number of samples in each beat

JTaus = cell(del*2+1,1);
X = CreateAverageBeat(Y,Fs,R,del,N);
YAtilde = AFReduction(Y,Fs,R);
Z = Y-YAtilde;
%Z = Y;
taus = -del:del;
Params = cell(Nbeats,numel(taus));

for i=1:Nbeats
    Zbeat = zeros(N,L);
    for k = 1:N
        for l = 1:L
            Zbeat(k,l) = Z((s(i)+k-1),l);
        end
    end
    for tau=1:numel(taus)
        count = 1;
        J1 = zeros(N, del+taus(tau));
        size(J1);
        J2 = eye(N);
        size(J2);
        J3 = zeros(N, del-taus(tau));
        size(J3);
        JTau = cat(2,J1,J2,J3);
        D = eye(L);
        JTaus{tau} = JTau;
        err = intmax;
        params = cell(MaxIt);
        while err>eps && count<=MaxIt
            T = D'*X'*JTau'*Zbeat;
            [U,~,V] = svd(T);
            Q = U*V';
            K = JTau*X;
            ZPrime = Zbeat*inv(Q);
            for j=1:L                
                D(j,j) = inv(K(:,j)'*K(:,j))*(K(:,j)'*ZPrime(:,j));
            end
            params{count} = Zbeat-(JTau*X*D*Q);
%             if(count>1)
%                 err = norm(params{count-1},'fro')-norm(params{count},'fro');%
%             else
%             end
            err=norm(params{count},'fro')^2;
            count = count+1;
        end
        Params{i,tau} = params{count-1};
    end
end
bestParams = cell(Nbeats,1);
norms = zeros(size(Params,1), size(Params,2));
for i=1:size(Params,1)
    for j=1:size(Params,2)
        norms(i,j) = norm(Params{i,j},'fro');
    end    
end
[~, I] = min(norms,[],2);%indices of the min norm, therefore the min JTau
for i=1:Nbeats
    bestParams(i) = Params(i,I(i));
end

for j=1:Nbeats
    Ya(s(j):s(j)+N-1,:)=bestParams{j,1}(:,:);
end
end
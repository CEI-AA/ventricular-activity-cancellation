function[s1mod,erro] = normalizag(s1,s2,P)

testegauss = find(s1 < 0);

if length(testegauss) == 0
ps1=find(s1==max(s1));
ofs1=ps1-1;
ss2=s2(max(P-ps1,1):min(P+ps1,length(s2)));

if s2(P) < 0
    ss2 = -1*ss2;
end



AmpMax=max(ss2);
AmpMin=min(ss2);
AmpmaxS1=max(s1);
AmpminS1=min(s1);
if AmpmaxS1-AmpminS1 ~=0
	s1mod=((AmpmaxS1-s1)/(AmpmaxS1-AmpminS1))*(AmpMax-AmpMin);
	s1mod=AmpMax-s1mod;
end

K = min(length(s1mod),length(ss2));
soma1 = 0;
soma2 = 0;
for i = 1:K
   soma1 = soma1 + (s1mod(i) - ss2(i))^2;
   soma2 = soma2 + (ss2(i))^2; 
end
erro = soma1/soma2;


else 
    
    ps1=find(s1==min(s1));
    ofs1=ps1-1;
    ss2=s2(max(P-ps1,1):min(P+ps1,length(s2)));





    AmpMax=max(ss2);
    AmpMin=min(ss2);
    AmpmaxS1=max(s1);
    AmpminS1=min(s1);
    if AmpmaxS1-AmpminS1 ~=0
        s1mod=((AmpmaxS1-s1)/(AmpmaxS1-AmpminS1))*(AmpMax-AmpMin);
        s1mod=AmpMax-s1mod;
    end

    K = min(length(s1mod),length(ss2));
    soma1 = 0;
    soma2 = 0;
    for i = 1:K
        soma1 = soma1 + (s1mod(i) - ss2(i))^2;
        soma2 = soma2 + (s1mod(i))^2; 
    end
    erro = soma1/soma2;
    
end

    
    
    
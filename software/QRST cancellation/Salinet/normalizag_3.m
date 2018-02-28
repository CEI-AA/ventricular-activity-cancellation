function[s1mod,erro] = normalizag(s1,s2,P)

ps1=find(s1==max(s1));
ofs1=ps1-1;
ss2=s2(max(P-ps1,1):min(P+(length(s1)-ps1),length(s2)));

if P == 0
    keyboard;
end

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


gpico = find(abs(s1mod)==max(abs(s1mod)));
irg = P-gpico;
pg = P;
xini = max(irg,1);
xfin = irg+length(s1)-1;
xinig = xini-irg + 1;
yfing = length(s1mod);
ws1mod = s1mod(xinig:yfing);

K = min(length(ws1mod),length(ss2));
soma1 = 0;
soma2 = 0;

for i = 1:K
   soma1 = soma1 + (ws1mod(i) - ss2(i))^2;
   soma2 = soma2 + (ss2(i))^2; 
end

erro = soma1/soma2;



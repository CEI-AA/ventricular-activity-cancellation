function[sout] = aligncc(s1,s2)

C = xcorr(s1,s2);
L = length(s1);
inicio = L;
vbol=0;
maxf = 0;
maxback = 0;
ii=inicio;

while vbol == 0
    if C(ii) > C(ii+1) && C(ii) > C(ii-1)
        vbol=1;
    else ii = ii+1;
    end
end

ii2 = inicio;
vbol=0;

while vbol == 0
    if C(ii2) > C(ii2+1) && C(ii2) > C(ii2-1)
        vbol=1;
    else ii2 = ii2-1;
    end
end

if abs(ii-inicio) < abs(ii2-inicio)
    delay = ii-inicio;
    s2m = zeros(1,length(s2));
    s2m(delay+1:end)=s2(1:end-delay);
else delay = inicio-ii2;
     s2m = zeros(1,length(s2));
     s2m(1:end-delay)=s2(delay+1:end);
end

%figure;plot(s1);hold on;plot(s2,'k');plot(s2m,'r');keyboard;close all; 
s2m=s2m';
sout = s1+s2m;

    
    

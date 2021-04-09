clc;
snrdb =-5:1:85;
snr = 10.^(snrdb/10);
lambda = 100;
r = 45e-9;
d = 500e-9;
D = 4.265e-10;
delta_T = 9e-6;
T = 30*delta_T;
L = 5;
global sum_Cj;
sum=0;
sum1=0;
x =0:1:65;
sum_Cj=0;
Co = zeros(1,length(snrdb));
s=0;
temp=0;
data1 = rand(1,L)>0.5;
display(data1);
P_0 = (r/d)*(erfc((d-r)/sqrt(4*D*T)));

Ntx = 2.*lambda.*T.*(10.^(snr./10))./P_0;
P_i1(1)=P_0;
final_term = 100000;
ri = zeros(90,5);
for ii=1:1:length(snrdb)
    for i = 2:5
        P_i1(i) = (r/d)*(erfc((d-r)/sqrt(4*D*i*T))-erfc((d-r)/sqrt(4*D*(i-1)*T)));
    end
    for j = 1:L
        cj(j) = Ntx(ii)*P_i1(j);
    end
    Co(ii) = snr(ii)*2*lambda*T;
    average = (lambda*T + sum1);
    average1 = average + data1.*Co(ii);
    ri(ii,:) = poissrnd(average1);
end

ri=reshape(ri', [5*90,1]);
data2=repmat(data1',90,1);
data2=reshape(data2,[5*90,1]);
writematrix(data2,'bit_detect.csv');
writematrix(ri,'ri.csv');
clc
close all

ri_final = [];
sj_final = [];
tau_final = [];
%snr = 40;
for snr = 0:50
  lambda0 = 100;
  T = 30*9*(10^-6);
  d = 500*(10^-9);
  r = 45*(10^-9);
  L=5;
  time_slot = [1:5000];
  capD = 4.265*(10^-10); % diffusion coefficient
  i=1;
  ri = zeros(1,5000);
  sj = randn(1, length(time_slot))>0.5;

  P_0 = (r/d)*(erfc((d-r)/sqrt(4*capD*i*T))-erfc((d-r)/sqrt(4*capD*(i-1)*T)));
  Ntx = 2.*lambda0.*T.*(10.^(snr./10))./P_0;
  c0 = Ntx*P_0;

  for i = 1:length(time_slot)
  	  P_i1(i) = (r/d)*(erfc((d-r)/sqrt(4*capD*i*T))-erfc((d-r)/sqrt(4*capD*(i-1)*T)));
  end

  for j = 1:length(time_slot)
	  cj(j) = Ntx.*P_i1(j);
  end

%   c0 = 54;
  sum1 = sum(sj.*cj);
  avg = (lambda0*T + sum1);
  avg1 = avg + sj.*c0;
  ri(1,:) = poissrnd(avg1);

  tau = c0./log(1 + (c0/(sum(cj/2)+(lambda0*T))));

  ri_final = [ri_final, ri];
  sj_final = [sj_final, sj];
  tau_final = [tau_final, tau];
end
%ri = reshape(ri, [1, 56*5])

%ri = reshape
% sj = repmat(sj, 1, 56)

A = [ri_final'; sj_final'];
writematrix(A, 'Data_MC_Ntx.csv');

%Ntx
%ri
%sj

index = 1;
bit_detct = [];

for j = 1:51
%disp(tau_final(j));
%disp(j)
    for i=1:100
        %ri(index) = tmp(i);
        if(ri_final(index)>tau_final(j))
             bit_detect(index) = 1;
        else
             bit_detect(index) =0;
    end
  index = index +1;
    end
end
bit_detect = bit_detect';
writematrix(bit_detect, 'Test_Bit_detect.csv');
%writematrix(sj, 'MC_Ntx.csv');

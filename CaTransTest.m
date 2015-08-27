%% Ca Transient Test
% Author: David Wyrick
% Created: 8/21/15
% Last Updated: 8/21/15
%%
clc;
Ca_0 = 0.1399; % 0.1399 Resting Diastolic Calcium concentration in micromolar
Ca_max = 0.483; %0.3431 Ca Transient Max concentration at Time to Peak (48.2ms)
TPT = 48.2; %ms
RT50 = 175.9; %Time to 50 percent relaxation
RT90 = 343.1;
tao = 48.2;
% 
% %% Piecewise function fitting
% % for time less than TPT
% fca1 = Ca_0 + (Ca_max - Ca_0).*(TPT/t).*exp(1-(TPT/t));
% % for time gt TPT but lt RT50
% fca2 = Ca_0 + (Ca_max - Ca_0).*(Tao1/(t-t1)).*exp(1-(Tao1/(t-t1)));
% % for time gt RT50 but lt 2*RT90
% fca2 = Ca_0 + (Ca_max - Ca_0).*(Tao2/(t-t2)).*exp(1-(Tao2/(t-t2)));
% % for time g
%%
tstep = 0:10:1000;
timebtwPulse=700;
Ca_t = Ca_0 + (Ca_max - Ca_0).*(tstep/tao).*exp(1-(tstep/tao)+timebtwPulse); %Ca Transient function

p = plot(tstep,Ca_t,'-.k')
xlabel('Time(ms)')
ylabel('[Ca2+]')
xlim([-50,1000])
ylim('auto')
title('Ca Transient Profile Test')
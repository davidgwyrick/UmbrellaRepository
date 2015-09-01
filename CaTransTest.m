%% Ca Transient Test
% Author: David Wyrick
% Created: 8/21/15
% Last Updated: 8/21/15
%%
clc;
Ca_0 = 0.1399*10^(-6); % 0.1399 Resting Diastolic Calcium concentration in micromolar
Ca_dif = 0.3431*10^(-6); %0.3431 Ca Transient Max concentration at Time to Peak (48.2ms)
TPT = 48.2; %ms
RT50 = 175.9; %Time to 50 percent relaxation
RT90 = 343.1;
tao = 48.2;
% 
% %% Piecewise function fitting
% % for time less than TPT
fca1 = Ca_0 + (Ca_dif).*(t/TPT).*exp(1-(t/TPT));
% % for time gt TPT but lt RT50
fca2 = Ca_0 + (Ca_dif).*((t-t1)/tao1).*exp(1-(t-t1)/tao1);
% % for time gt RT50 but lt 2*RT90
fca3 = Ca_0 + (Ca_dif).*((t-t2)/tao2).*exp(1-(t-t2)/tao2);
% % for time gt RT90 but lt 1000
fca4 = a-b*t
%%
tstep = 0:10:1000;
timebtwPulse=700;
Ca_t = Ca_0 + (Ca_dif).*(tstep/tao).*exp(1-(tstep/tao)+timebtwPulse); %Ca Transient function

p = plot(tstep,Ca_t,'-.k')
xlabel('Time(ms)')
ylabel('[Ca2+]')
xlim([-50,1000])
ylim('auto')
title('Ca Transient Profile Test')
% 
%         CaChange_index = (0.25*10^-7)*ones(NSTEPS,1);
%         Ca_0 = 0.14*10^(-6); % 0.1399*10^(-6)molar Resting Diastolic Calcium concentration 
%         Ca_max =10^(-New_pCa); %0.3431*10^(-6)molar is pCa of 6.46 Ca Transient Max concentration at Time to Peak (48.2ms)
%         tao=.5*pulse_width;
%         for indy=1:NSTEPS  
%             CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*(indy/tao)*exp(1-(indy/tao));
%         end
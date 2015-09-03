%% Ca Transient Test
% Author: David Wyrick
% Created: 8/21/15
% Last Updated: 8/21/15
%%
clc;
% Ca_0 = 0.1399*10^(-6); % 0.1399 Resting Diastolic Calcium concentration in micromolar
% Ca_dif = 0.3431*10^(-6); %0.3431 Ca Transient Max concentration at Time to Peak (48.2ms)
% TPT = 48.2; %ms
% RT50 = 175.9; %Time to 50 percent relaxation
% RT90 = 343.1;
% tao = 48.2;
% % 
% % %% Piecewise function fitting
% % % for time less than TPT
% fca1 = Ca_0 + (Ca_dif).*(t/TPT).*exp(1-(t/TPT));
% % % for time gt TPT but lt RT50
% fca2 = Ca_0 + (Ca_dif).*((t-t1)/tao1).*exp(1-(t-t1)/tao1);
% % % for time gt RT50 but lt 2*RT90
% fca3 = Ca_0 + (Ca_dif).*((t-t2)/tao2).*exp(1-(t-t2)/tao2);
% % % for time gt RT90 but lt 1000
% fca4 = a-b*t
%%
% 
% Tao=[5 10 15 20 25 50 75 90 100 150 200 250 400]
% tstep = 0:10:2000;
% timebtwPulse=700;
% 
% Ca_t = Ca_0 + (Ca_dif).*(tstep/tao).*exp(1-(tstep/tao)); %Ca Transient function
% 
% p = plot(tstep,Ca_t,'-.k')
% xlabel('Time(ms)')
% ylabel('[Ca2+]')
% xlim([-50,1000])
% ylim('auto')
% title('Ca Transient Profile Test')
% 
%         CaChange_index = (0.25*10^-7)*ones(NSTEPS,1);
%         Ca_0 = 0.14*10^(-6); % 0.1399*10^(-6)molar Resting Diastolic Calcium concentration 
%         Ca_max =10^(-New_pCa); %0.3431*10^(-6)molar is pCa of 6.46 Ca Transient Max concentration at Time to Peak (48.2ms)
%         tao=.5*pulse_width;
%         for indy=1:NSTEPS  
%             CaChange_index(indy) = Ca_0 + (Ca_max - Ca_0)*(indy/tao)*exp(1-(indy/tao));
%         end

%% 
pCa_max=4.0;
tao_cell = zeros(100,2);
Ca_0=0.14*10^(-6); %*10^(-6)
 %Ca_dif=0.3431*10^(-6);
Ca_max=10^(-pCa_max)
for ii=1:100,
    
    tao_cell(ii,1)=ii*5.0;
    tao=ii*5.0*10^(-3); %Tao parameter
    tstep=0:0.001:2.0; %1ms timestep 
    %tstep = 0:1:2000;
    Ca_t = Ca_0 + (Ca_max - Ca_0).*(tstep/tao).*exp(1-(tstep/tao)); %Ca Transient function
%     plot(tstep,Ca_t,'-.k')
%     plot(tstep,-log10(Ca_t),'-.k')
%     axis([-.1 2.0 5 7],'on')
%     set(gca, 'ydir', 'reverse')
%     xlabel('Time(ms)')
%     ylabel('[Ca2+]')
% 
%     title('Ca Transient Profile Test')
    PW_found=0;
    [mm mm_index] =max(Ca_t);
    for jj=1:length(Ca_t),
        if jj > mm_index && PW_found == 0 &&  or(Ca_t(jj) < .01*mm,Ca_t(jj) <1.01*Ca_0)
            tao_cell(ii,2)=jj; %Time for Ca profile curve to reach 97.5% relaxation
            PW_found=1;
        end
    end
end
del_ind=0;
for ii=1:100,
    if tao_cell(ii,2) == 0 && del_ind==0, del_ind=ii;end
end
tao_param=tao_cell(1:(del_ind-1),1);
PW_val=tao_cell(1:(del_ind-1),2);
% plot(PW_val,tao_param,'k.')
% ylabel('Tao parameter')
% xlabel('Actual Pulse Width of Calcium profile')
% Title(['Pulse Width - Tao relationship for pCa of ' num2str(pCa_max)])

f=fittype('a*x+b');
[c,gof,outp]=fit(PW_val,tao_param,f);
slope=c.a;
yint=c.b;

%% Making a pulse 
pulse_width=90;

tao=(0.12805*pulse_width-0.20183)*10^(-3);
tstep=0:0.001:2.0; %1ms timestep 
Ca_t = Ca_0 + (Ca_max - Ca_0).*(tstep/tao).*exp(1-(tstep/tao)); 
plot(tstep,-log10(Ca_t),'-.k')
axis([-.1 2.0 5 7],'on')
set(gca, 'ydir', 'reverse')
xlabel('Time(ms)')
ylabel('[Ca2+]')

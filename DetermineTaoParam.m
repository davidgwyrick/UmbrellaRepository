function [Tao_final] = DetermineTaoParam(Pulse_Width,pCa_max)
%To make a Ca transient profile of a given pulse_width, we need to find the
%Tao parameter of the equation: Ca_t = Ca_0 + (Ca_max - Ca_0).*(tstep/tao).*exp(1-(tstep/tao));
%This function determines this tao parameter based on the pulse_width
%desired. Pulse_Width input must be a multiple of 5
%Created by: David Wyrick
%9/3/15

tao_cell = zeros(100,2);
Ca_0=0.14*10^(-6); %*10^(-6)
Ca_max=10^(-pCa_max);

for ii=1:100,
    
    tao_cell(ii,1)=ii*5.0; 
    tao=ii*5.0*10^(-3); %Tao parameter
    tstep=0:0.001:2.0; %1ms timestep 

    Ca_t = Ca_0 + (Ca_max - Ca_0).*(tstep/tao).*exp(1-(tstep/tao)); %Ca Transient function
%     plot(tstep,Ca_t,'-.k')
%     plot(tstep,-log10(Ca_t),'-.k')
%     axis([-.1 2.0 5 7],'on')
%     set(gca, 'ydir', 'reverse')
%     xlabel('Time(ms)')
%     ylabel('[Ca2+]')
%     title('Ca Transient Profile Test')
    PW_found=0;
    [mm mm_index] =max(Ca_t);
    for jj=1:length(Ca_t),
        if jj > mm_index && PW_found == 0 &&  or(Ca_t(jj) < .01*mm,Ca_t(jj) <1.01*Ca_0)
            tao_cell(ii,2)=jj; %Time for Ca profile curve to reach 99% relaxation
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
Tao_final=(c.a)*Pulse_Width+c.b;





end


%% Ca Transient Profile
% 9/9/15
%
%% Load Data

filename=['HumanCai.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
humanCai=importdata(filepath);
filename=['mouseCai.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
mouseCai=importdata(filepath);
filename=['ratCai.txt'];
filepath=['J:\TannerBertGroup\Sims\David_Wyrick\UmbrellaRepository\DataFiles\CaTransient\',filename];
ratCai=importdata(filepath);

%% 
NSTEPS=1001
CaWyHuman = (0.25*10^-7)*ones(NSTEPS,1);
Ca_0 = 0.14; % 0.1399*10^(-6)molar Resting Diastolic Calcium concentration 
Ca_max =.5; %0.3431*10^(-6)molar is pCa of 6.46 Ca Transient Max concentration at Time to Peak (48.2ms)
 tao=DetermineTaoParam(600,6.3);
 disp('Human')
 disp(tao)
  tao=80;
for indy=1:NSTEPS 
    CaWyHuman(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
end
CaWyHuman=CaWyHuman';
%%
NSTEPS=168;
CaWyMouse=ones(NSTEPS,1);
Ca_0 = 0.199;
Ca_max=.4965;
 tao=DetermineTaoParam(55,6.3);
  disp('Mouse')
 disp(tao)
tao=25;
for indy=1:NSTEPS 
    CaWyMouse(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
end
CaWyMouse=CaWyMouse';
%%
NSTEPS=168;
CaWyRat=ones(NSTEPS,1);
Ca_0 = 0.237;
Ca_max=1.423;
 tao=DetermineTaoParam(65,5.85);
  disp('Rat')
 disp(tao)
tao=25;
for indy=1:NSTEPS 
    CaWyRat(indy) = Ca_0 + (Ca_max - Ca_0)*((indy-1)/tao)*exp(1-(indy/tao));
end
CaWyRat=CaWyRat';
%%

figure(1); clf; hold on; 
plot(0:1001,[mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1) mouseCai(1:end-1)],'-','Color',[0.8 0.2 0.2],'LineWidth',2);
plot(0:1001,[ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1) ratCai(1:end-1)],'-','Color',[0.2 0.2 0.8],'LineWidth',2);
plot(0:1000,humanCai,'-','Color',[0.2 0.8 0.2],'LineWidth',2);


plot(0:1001,[CaWyMouse(1:end-1) CaWyMouse(1:end-1) CaWyMouse(1:end-1) CaWyMouse(1:end-1) CaWyMouse(1:end-1) CaWyMouse(1:end-1)],'-c','LineWidth',2);
plot(0:1001,[CaWyRat(1:end-1) CaWyRat(1:end-1) CaWyRat(1:end-1) CaWyRat(1:end-1) CaWyRat(1:end-1) CaWyRat(1:end-1)],'-m','LineWidth',2);
plot(0:1000,CaWyHuman,'-k','LineWidth',2);

xlim([0 1000]); 
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i \mu{}M');
legend('Mouse', 'Rat', 'Human','WyMouse','WyRat','WyHuman');
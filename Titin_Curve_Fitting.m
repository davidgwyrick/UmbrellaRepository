%% Examples of Fitting Titin Data

close
clear all
clc

soleus.SL=[1.9, 2.5, 3, 3.2, 3.5, 3.9]; %% in um
soleus.HS_Force=[0, 12, 30, 48, 78, 240]; %% in pN
soleus.Titin_Length=[90, 390, 640, 740, 890, 1090]; %% in nm
soleus.Titin_Force=[0, 0.5, 1.25, 2, 3.25, 10]; %% in pN
soleus.TF_Force=[0, 3, 7.5, 12, 19.5, 60]; %% in pN

psoas.SL=[1.9, 2.3, 2.6, 2.8, 3.0, 3.1]; %% in um
psoas.HS_Force=[0, 14.4, 24, 50.4, 122.4, 240]; %% in pN
psoas.Titin_Length=[90, 290, 440, 540, 640, 690]; %% in nm
psoas.Titin_Force=[0, 0.6, 1.0, 2.1, 5.1, 10]; %% in pN
psoas.TF_Force=[0, 3.6, 6, 12.6, 30.6, 60]; %% in pN

N2B.SL=[1.9, 2.2, 2.3, 2.4, 2.5, 2.55]; %% in um
N2B.HS_Force=[0, 12, 36, 72, 168, 240]; %% in pN
N2B.Titin_Length=[90, 240, 290, 340, 390, 415]; %% in nm
N2B.Titin_Force=[0, 0.5, 1.5, 3, 7, 10]; %% in pN
N2B.TF_Force=[0, 3, 9, 18, 42, 60]; %% in pN

Rest_Length_SL = 1.9;
Rest_Lenth_Titin = 90;

%% Plot our curves
clf(figure(1))
%%Soleus
subplot(3,3,3)
plot(soleus.SL, soleus.HS_Force, 'ko')
ylabel('HS Passive Force (pN)'), hold on
xlabel('Sarc. Length (um)')
set(gca, 'box', 'off')
title('Soleus: Force per Half-Sarcomere')

subplot(3,3,2)
plot(soleus.Titin_Length, soleus.TF_Force, 'ko')
ylabel('TF Passive Force (pN)'), hold on
xlabel('Titin Length (um)')
set(gca, 'box', 'off')
title('Soleus: Force per Thick Filament')

subplot(3,3,1)
plot(soleus.Titin_Length, soleus.Titin_Force, 'k*'), hold on
ylabel('Titin Force (pN)')
xlabel('Titin Length (nm)')
set(gca, 'box', 'off')
title('Soleus: Force per Titin')

%%Psoas
subplot(3,3,6)
plot(psoas.SL, psoas.HS_Force, 'ko')
ylabel('HS Passive Force (pN)'), hold on
xlabel('Sarc. Length (um)')
set(gca, 'box', 'off')
title('Psoas: Force per Half-Sarcomere')

subplot(3,3,5)
plot(psoas.Titin_Length, psoas.TF_Force, 'ko')
ylabel('TF Passive Force (pN)'), hold on
xlabel('Titin Length (um)')
set(gca, 'box', 'off')
title('Psoas: Force per Thick Filament')

subplot(3,3,4)
plot(psoas.Titin_Length, psoas.Titin_Force, 'k*'), hold on
ylabel('Titin Force (pN)')
xlabel('Titin Length (nm)')
set(gca, 'box', 'off')
title('Psoas: Force per Titin')

%%N2B
subplot(3,3,9)
plot(N2B.SL, N2B.HS_Force, 'ko')
ylabel('HS Passive Force (pN)'), hold on
xlabel('Sarc. Length (um)')
set(gca, 'box', 'off')
title('N2B: Force per Half-Sarcomere')

subplot(3,3,8)
plot(N2B.Titin_Length, N2B.TF_Force, 'ko')
ylabel('TF Passive Force (pN)'), hold on
xlabel('Titin Length (um)')
set(gca, 'box', 'off')
title('N2B: Force per Thick Filament')

subplot(3,3,7)
plot(N2B.Titin_Length, N2B.Titin_Force, 'k*'), hold on
ylabel('Titin Force (pN)')
xlabel('Titin Length (nm)')
set(gca, 'box', 'off')
title('N2B: Force per Titin')

%% Soleus
SL = soleus.SL;
HS_Force = soleus.HS_Force;
Titin_Length = soleus.Titin_Length;
Titin_Force = soleus.Titin_Force;
TF_Force = soleus.TF_Force;

%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess=[0.0673; ... % Sigma
    90; ... % Loffset
    200.2493]; % HalfLength

%% Now use this to fit the Titin as a 3 free parameter fit
%% using lsqcurvefit.  First define the inline function
Fun3=inline('x(1,1)*exp((xdata-x(2,1))/x(3,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
%x_Fit3 = lsqcurvefit(@myfun3, x_Guess, xdata, ydata)
x_Fit3 = lsqcurvefit(Fun3, x_Guess, Titin_Length, Titin_Force)


%% Now use this to fit the Titin as a 2 free parameter fit
%% using lsqcurvefit.  This would be using the 2 parameter fit, saying that
%% we are going to fix tithe L_offset=90.  This could be passed in, using
%% additional function inputs--or it can be fixed at = 90 as I've writtin it
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-90)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2 = lsqcurvefit(Fun2, [x_Guess(1,1); x_Guess(3,1)], Titin_Length, Titin_Force)

%% Now calculate our Titin values
x_Titin=[0:10:1200]; %% Just to get a better resolution
Titin_F_Guess=x_Guess(1,1)*exp((x_Titin-x_Guess(2,1))/x_Guess(3,1));
Titin_F3_Param=x_Fit3(1,1)*exp((x_Titin-x_Fit3(2,1))/x_Fit3(3,1));
Titin_F2_Param=x_Fit2(1,1)*exp((x_Titin-90)/x_Fit2(2,1));

subplot(3,3,1)
plot(x_Titin, Titin_F_Guess, 'b-')
plot(x_Titin, Titin_F3_Param, 'g-')
plot(x_Titin, Titin_F2_Param, 'k-')
%legend('Data', 'Guess', ['3-param fit=' num2str(x_Fit3(1,1) 'exp[(x-' num2str(x_Fit3(2,1) ')/' num2str(x_Fit3(3,1) ']'], '2-param fit', 'location', 'Northwest')
legend('Data', 'Guess', '3-param fit', '2-param fit', 'location', 'Northwest')


%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess_TF=[0.4037; ... % Sigma
    90; ... % Loffset
    200.2493]; % HalfLength
% 
% %% Now use this to fit the Titin as a 3 free parameter fit
% %% using lsqcurvefit.  First define the inline function
x_Fit3_TF = lsqcurvefit(Fun3, x_Guess_TF, Titin_Length, TF_Force)
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-90)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2_TF = lsqcurvefit(Fun2, [x_Guess_TF(1,1); x_Guess_TF(3,1)], Titin_Length, TF_Force)

%% Now calculate our Titin values
x_Titin=[0:10:1200]; %% Just to get a better resolution
TF_F_Guess=x_Guess_TF(1,1)*exp((x_Titin-x_Guess_TF(2,1))/x_Guess_TF(3,1));
TF_Fit3=x_Fit3_TF(1,1)*exp((x_Titin-x_Fit3_TF(2,1))/x_Fit3_TF(3,1));
TF_Fit2=x_Fit2_TF(1,1)*exp((x_Titin-90)/x_Fit2_TF(2,1));

subplot(3,3,2)
plot(x_Titin, TF_F_Guess, 'b-')
plot(x_Titin, TF_Fit3, 'g-')
plot(x_Titin, TF_Fit2, 'k-')


%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess_HS=[1.6149; ... % Sigma
    1.9; ... % Loffset
    0.4005]; % HalfLength
% 
% %% Now use this to fit the Titin as a 3 free parameter fit
% %% using lsqcurvefit.  First define the inline function
x_Fit3_HS = lsqcurvefit(Fun3, x_Guess_HS, SL, HS_Force)
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-1.9)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2_HS = lsqcurvefit(Fun2, [x_Guess_HS(1,1); x_Guess_HS(3,1)], SL, HS_Force)

%% Now calculate our Titin values
x_SL=[1.5:0.1:4]; %% Just to get a better resolution
HS_F_Guess=x_Guess_HS(1,1)*exp((x_SL-x_Guess_HS(2,1))/x_Guess_HS(3,1));
HS_Fit3=x_Fit3_HS(1,1)*exp((x_SL-x_Fit3_HS(2,1))/x_Fit3_HS(3,1));
HS_Fit2=x_Fit2_HS(1,1)*exp((x_SL-1.9)/x_Fit2_HS(2,1));

subplot(3,3,3)
plot(x_SL, HS_F_Guess, 'b-')
plot(x_SL, HS_Fit3, 'g-')
plot(x_SL, HS_Fit2, 'k-')

%% Psoas
SL = psoas.SL;
HS_Force = psoas.HS_Force;
Titin_Length = psoas.Titin_Length;
Titin_Force = psoas.Titin_Force;
TF_Force = psoas.TF_Force;

%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess=[0.0129; ... % Sigma
    90; ... % Loffset
    90.5240]; % HalfLength

%% Now use this to fit the Titin as a 3 free parameter fit
%% using lsqcurvefit.  First define the inline function
Fun3=inline('x(1,1)*exp((xdata-x(2,1))/x(3,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
%x_Fit3 = lsqcurvefit(@myfun3, x_Guess, xdata, ydata)
x_Fit3 = lsqcurvefit(Fun3, x_Guess, Titin_Length, Titin_Force)


%% Now use this to fit the Titin as a 2 free parameter fit
%% using lsqcurvefit.  This would be using the 2 parameter fit, saying that
%% we are going to fix tithe L_offset=90.  This could be passed in, using
%% additional function inputs--or it can be fixed at = 90 as I've writtin it
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-90)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2 = lsqcurvefit(Fun2, [x_Guess(1,1); x_Guess(3,1)], Titin_Length, Titin_Force)

%% Now calculate our Titin values
x_Titin=[0:10:800]; %% Just to get a better resolution
Titin_F_Guess=x_Guess(1,1)*exp((x_Titin-x_Guess(2,1))/x_Guess(3,1));
Titin_F3_Param=x_Fit3(1,1)*exp((x_Titin-x_Fit3(2,1))/x_Fit3(3,1));
Titin_F2_Param=x_Fit2(1,1)*exp((x_Titin-90)/x_Fit2(2,1));

subplot(3,3,4)
plot(x_Titin, Titin_F_Guess, 'b-')
plot(x_Titin, Titin_F3_Param, 'g-')
plot(x_Titin, Titin_F2_Param, 'k-')
%legend('Data', 'Guess', ['3-param fit=' num2str(x_Fit3(1,1) 'exp[(x-' num2str(x_Fit3(2,1) ')/' num2str(x_Fit3(3,1) ']'], '2-param fit', 'location', 'Northwest')
legend('Data', 'Guess', '3-param fit', '2-param fit', 'location', 'Northwest')

%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess_TF=[0.0776; ... % Sigma
    90; ... % Loffset
    90.5243]; % HalfLength
% 
% %% Now use this to fit the Titin as a 3 free parameter fit
% %% using lsqcurvefit.  First define the inline function
x_Fit3_TF = lsqcurvefit(Fun3, x_Guess_TF, Titin_Length, TF_Force)
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-90)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2_TF = lsqcurvefit(Fun2, [x_Guess_TF(1,1); x_Guess_TF(3,1)], Titin_Length, TF_Force)

%% Now calculate our Titin values
x_Titin=[0:10:800]; %% Just to get a better resolution
TF_F_Guess=x_Guess_TF(1,1)*exp((x_Titin-x_Guess_TF(2,1))/x_Guess_TF(3,1));
TF_Fit3=x_Fit3_TF(1,1)*exp((x_Titin-x_Fit3_TF(2,1))/x_Fit3_TF(3,1));
TF_Fit2=x_Fit2_TF(1,1)*exp((x_Titin-90)/x_Fit2_TF(2,1));

subplot(3,3,5)
plot(x_Titin, TF_F_Guess, 'b-')
plot(x_Titin, TF_Fit3, 'g-')
plot(x_Titin, TF_Fit2, 'k-')


%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess_HS=[0.3102; ... % Sigma
    1.9; ... % Loffset
    0.1810]; % HalfLength
% 
% %% Now use this to fit the Titin as a 3 free parameter fit
% %% using lsqcurvefit.  First define the inline function
x_Fit3_HS = lsqcurvefit(Fun3, x_Guess_HS, SL, HS_Force)
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-1.9)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2_HS = lsqcurvefit(Fun2, [x_Guess_HS(1,1); x_Guess_HS(3,1)], SL, HS_Force)

%% Now calculate our Titin values
x_SL=[1.5:0.1:3.3]; %% Just to get a better resolution
HS_F_Guess=x_Guess_HS(1,1)*exp((x_SL-x_Guess_HS(2,1))/x_Guess_HS(3,1));
HS_Fit3=x_Fit3_HS(1,1)*exp((x_SL-x_Fit3_HS(2,1))/x_Fit3_HS(3,1));
HS_Fit2=x_Fit2_HS(1,1)*exp((x_SL-1.9)/x_Fit2_HS(2,1));

subplot(3,3,6)
plot(x_SL, HS_F_Guess, 'b-')
plot(x_SL, HS_Fit3, 'g-')
plot(x_SL, HS_Fit2, 'k-')



%% N2B
SL = N2B.SL;
HS_Force = N2B.HS_Force;
Titin_Length = N2B.Titin_Length;
Titin_Force = N2B.Titin_Force;
TF_Force = N2B.TF_Force;

%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess=[0.0619; ... % Sigma
    90; ... % Loffset
    63.8005]; % HalfLength

%% Now use this to fit the Titin as a 3 free parameter fit
%% using lsqcurvefit.  First define the inline function
Fun3=inline('x(1,1)*exp((xdata-x(2,1))/x(3,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
%x_Fit3 = lsqcurvefit(@myfun3, x_Guess, xdata, ydata)
x_Fit3 = lsqcurvefit(Fun3, x_Guess, Titin_Length, Titin_Force)


%% Now use this to fit the Titin as a 2 free parameter fit
%% using lsqcurvefit.  This would be using the 2 parameter fit, saying that
%% we are going to fix tithe L_offset=90.  This could be passed in, using
%% additional function inputs--or it can be fixed at = 90 as I've writtin it
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-90)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2 = lsqcurvefit(Fun2, [x_Guess(1,1); x_Guess(3,1)], Titin_Length, Titin_Force)

%% Now calculate our Titin values
x_Titin=[0:10:500]; %% Just to get a better resolution
Titin_F_Guess=x_Guess(1,1)*exp((x_Titin-x_Guess(2,1))/x_Guess(3,1));
Titin_F3_Param=x_Fit3(1,1)*exp((x_Titin-x_Fit3(2,1))/x_Fit3(3,1));
Titin_F2_Param=x_Fit2(1,1)*exp((x_Titin-90)/x_Fit2(2,1));

subplot(3,3,7)
plot(x_Titin, Titin_F_Guess, 'b-')
plot(x_Titin, Titin_F3_Param, 'g-')
plot(x_Titin, Titin_F2_Param, 'k-')
%legend('Data', 'Guess', ['3-param fit=' num2str(x_Fit3(1,1) 'exp[(x-' num2str(x_Fit3(2,1) ')/' num2str(x_Fit3(3,1) ']'], '2-param fit', 'location', 'Northwest')
legend('Data', 'Guess', '3-param fit', '2-param fit', 'location', 'Northwest')

%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess_TF=[0.3713; ... % Sigma
    90; ... % Loffset
    63.8005]; % HalfLength
% 
% %% Now use this to fit the Titin as a 3 free parameter fit
% %% using lsqcurvefit.  First define the inline function
x_Fit3_TF = lsqcurvefit(Fun3, x_Guess_TF, Titin_Length, TF_Force)
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-90)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2_TF = lsqcurvefit(Fun2, [x_Guess_TF(1,1); x_Guess_TF(3,1)], Titin_Length, TF_Force)



%% Now calculate our Titin values
x_Titin=[0:10:500]; %% Just to get a better resolution
TF_F_Guess=x_Guess_TF(1,1)*exp((x_Titin-x_Guess_TF(2,1))/x_Guess_TF(3,1));
TF_Fit3=x_Fit3_TF(1,1)*exp((x_Titin-x_Fit3_TF(2,1))/x_Fit3_TF(3,1));
TF_Fit2=x_Fit2_TF(1,1)*exp((x_Titin-90)/x_Fit2_TF(2,1));

subplot(3,3,8)
plot(x_Titin, TF_F_Guess, 'b-')
plot(x_Titin, TF_Fit3, 'g-')
plot(x_Titin, TF_Fit2, 'k-')


%% Now fit paramters with equation
%% Force = Sigma*exp((x-Loffset)/HalfLength)
x_Guess_HS=[1.4851; ... % Sigma
    1.9; ... % Loffset
    0.1276]; % HalfLength
% 
% %% Now use this to fit the Titin as a 3 free parameter fit
% %% using lsqcurvefit.  First define the inline function
x_Fit3_HS = lsqcurvefit(Fun3, x_Guess_HS, SL, HS_Force)
%% here.  Note there are now only 2 parameters to optimize for X.
Fun2=inline('x(1,1)*exp((xdata-1.9)/x(2,1))', 'x', 'xdata');
%% Next Call the curfit on the usage data for 
x_Fit2_HS = lsqcurvefit(Fun2, [x_Guess_HS(1,1); x_Guess_HS(3,1)], SL, HS_Force)

%% Now calculate our Titin values
x_SL=[1.5:0.1:2.6]; %% Just to get a better resolution
HS_F_Guess=x_Guess_HS(1,1)*exp((x_SL-x_Guess_HS(2,1))/x_Guess_HS(3,1));
HS_Fit3=x_Fit3_HS(1,1)*exp((x_SL-x_Fit3_HS(2,1))/x_Fit3_HS(3,1));
HS_Fit2=x_Fit2_HS(1,1)*exp((x_SL-1.9)/x_Fit2_HS(2,1));

subplot(3,3,9)
plot(x_SL, HS_F_Guess, 'b-')
plot(x_SL, HS_Fit3, 'g-')
plot(x_SL, HS_Fit2, 'k-')
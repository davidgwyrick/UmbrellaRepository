%% BCWT and Mark McNabb
%% July 31, 2010
%%
%%
%% Import Data from a Simulation Set
%% Plot the force-pCa Data via function ton_pCa_Data_2013_7_31
%% Plot the ATP-pCa Data
%% Run the Time-On Analysis for a user defined pCa value
%% Write the Outfile for the data that will be used to plot/analyize
%% 
%%
clear
close all
clc

cd('J:\TannerBertGroup\Sims\Invited_MS_2013\WSU_Analysis')

%%-------------------------
SetNumber=140; 
pCa_In=4.0;%% User Enter Number for pCa to import
%pCa_In=5.8;%% User Enter Number
%%-------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User input the Set and pCa level that will be analyzed
% Assuming data at current directory-->ton_Data:
DataDir=['J:\TannerBertGroup\Sims\Vert_3Simple\DataFiles\'];
% Assuming output directory is the current directory-->\ton_Data:
Out_ton_Dir=([pwd filesep 'ton_Data' filesep]);

% Import Matlab Data Structure from the Simulation
InSet=[DataDir 'Set' num2str(SetNumber) '\']
InMatFile=sprintf('%sSimulationData.mat', InSet);
load(InMatFile)

%% Plot the Force, ATP, and Force/ATPase vs pCa at the defined figure:
Fig_Num=1;
View_3SimpleHillCurves_2013_7_31(pCaRange, OutTotal, Fig_Num, SetNumber, Tn_den, XB_den)

% %% Now Gather Force, ATP, Force/ATPase, ton data at a particular pCa
% %% level
% %% Also, output/plot a figure to look at the ton details
Fig_Num=3;
ton_pCa_Data_2013_7_31(OutTotal, Fig_Num, InSet, pCa_In, Out_ton_Dir, XB_den, SetNumber)


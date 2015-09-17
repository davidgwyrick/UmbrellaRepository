function [Time2Relax] = TwitchAnalysis(OutDir,pCa)
%Simple function to analyze twitch parameters like Time to Peak, time to
%50% and 90% relaxation, and max force. I use the Averaged max force and
%index instead of raw max
%David Wyrick
%9/16/15

%% Loads the data
filepath=[OutDir,'TimeSeriesAvg_pCa_',num2str(pCa,'%1.2f'),'.txt'];
Data=importdata(filepath);
TS_Time=Data.data(:,1);
TS_Force=Data.data(:,3);
TS_CaLevel=Data.data(:,12);
TS_TF3=Data.data(:,10);

%% Finds Max values, peak indices, T2P, RT90, and RT50
f7= ones(1,7)/7; %Average over 7 data points or 7ms
Time2Relax=zeros(2,8); %Matrix to store relaxation times for force and Ca2 fraction (Time to 50% activation, Time to 70% activation, T2P, RT50 along decay, RT70 along decay, RT90 along decay, Averaged max force, raw max force)(Force-row 1, Ca2-row 2)

% For the force transient
filtforce=filtfilt(f7,1,TS_Force);
[AvgMax avgMax_index]=max(filtforce);
Time2Relax(1,7)=AvgMax;
[Time2Relax(1,8) rawMax_index]=max(TS_Force);
Time2Relax(1,1)=TS_Time(find(TS_Force > 0.5*AvgMax,1,'first')); %Time to 50% activation
Time2Relax(1,2)=TS_Time(find(TS_Force > 0.7*AvgMax,1,'first')); %Time to 70% activation
Time2Relax(1,3)=TS_Time(avgMax_index); %Time to peak
Time2Relax(1,4)=TS_Time(find(TS_Force > 0.5*AvgMax,1,'last')); %Finds the last index of TS_Force that is 50% of the max and returns that index #, which is inputted into TS_Time
Time2Relax(1,5)=TS_Time(find(TS_Force > 0.3*AvgMax,1,'last')); %Time to 70% relaxation
Time2Relax(1,6)=TS_Time(find(TS_Force > 0.1*AvgMax,1,'last')); %Time to 90% relaxation

% For the Ca2 Fraction curve
filtCa2=filtfilt(f7,1,TS_TF3);
[AvgMax avgMax_index]=max(filtCa2);
Time2Relax(2,7)=AvgMax;
[Time2Relax(2,8) rawMax_index]=max(TS_TF3);
Time2Relax(2,1)=TS_Time(find(TS_TF3 > 0.5*AvgMax,1,'first')); %Time to 50% activation
Time2Relax(2,2)=TS_Time(find(TS_TF3 > 0.9*AvgMax,1,'first')); %Time to 90% activation
Time2Relax(2,3)=TS_Time(avgMax_index); %Time to peak
Time2Relax(2,4)=TS_Time(find(TS_TF3 > 0.5*AvgMax,1,'last')); %Finds the last index of TS_Force that is 50% of the max and returns that index #, which is inputted into TS_Time
Time2Relax(2,5)=TS_Time(find(TS_TF3 > 0.3*AvgMax,1,'last')); %Time to 70% relaxation
Time2Relax(2,6)=TS_Time(find(TS_TF3 > 0.1*AvgMax,1,'last')); %Time to 90% relaxation


%% Write to file in the OutDir given
OutFile=sprintf('%sTimeToRelax_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
fid=fopen(OutFile, 'wt' );
fprintf(fid, 'Tt50 Act\tTt70 Act\tTtPeak\tTt50 Relax\tTt70 Relax\tTt90 Relax\tAvgMax\tRawMax\n');	%header for outfile
FormatString=[];
[~, ColOut]=size(Time2Relax);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
fprintf(fid, FormatString, Time2Relax(1,:));
fprintf(fid, FormatString, Time2Relax(2,:));
fclose(fid);


end


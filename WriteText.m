function [] = WriteText(OutDir, pCa, dt, Binder, Steps, Stats, IndexThalf)

Stat_Names =  fieldnames(Stats);

%%  Just copying this here for reference, Idea being I can grab the state
%   names using this as an index -Axel
% Stat_Names =   
%     'MFMean'
%     'AFMean'
%     'FractXB1Mean'
%     'FractXB2Mean'
%     'FractCa0Mean'
%     'FractCa1Mean'
%     'FractCa2Mean'
%     'ATPuseMean'
%     'FractBoundMean'
%     'MFVar'
%     'AFVar'
%     'FractXB1Var'
%     'FractXB2Var'
%     'FractCa0Var'
%     'FractCa1Var'
%     'FractCa2Var'
%     'ATPuseVar'
%     'FractBoundVar'

TotalRuns = size(Stats.MFMean, 2);
TotalROIs = size(Stats.MFMean, 3);
ROIsAZ = 'A':'Z';

NSTEPS = size(Steps,2);

t = (dt*(1:NSTEPS))-dt;

% % Package Up the Averaged (Index) information into a vector.
% % Also, flip to column vectors, so we can write it to a outfile
% OutIndex=zeros(TotalRuns, 20);
% %% 24 rows for the first half results
% OutIndex(:, 1)=(1:TotalRuns)';
% OutIndex(:, 2)=pCa*ones(TotalRuns, 1);
% %% Concat the SS mean data, followed by Variance
% %for SS data before or without length change
% OutIndex(:, 3)=Means(1,:)';
% OutIndex(:, 4)=Vars(1,:)';
% OutIndex(:, 5)=Means(2,:)';
% OutIndex(:, 6)=Vars(2,:)';
% OutIndex(:, 7)=Means(3,:)';
% OutIndex(:, 8)=Vars(3,:)';
% OutIndex(:, 9)=Means(4,:)';
% OutIndex(:, 10)=Vars(4,:)';
% OutIndex(:, 11)=Means(5,:)';
% OutIndex(:, 12)=Vars(5,:)';
% OutIndex(:, 13)=Means(6,:)';
% OutIndex(:, 14)=Vars(6,:)';
% OutIndex(:, 15)=Means(7,:)';
% OutIndex(:, 16)=Vars(7,:)';
% OutIndex(:, 17)=Means(8,:)';
% OutIndex(:, 18)=Vars(8,:)';
% OutIndex(:, 19)=Means(9,:)';
% OutIndex(:, 20)=Vars(9,:)';


%The "Binder" that arrives is [1 x #runs] array, which then has sub
%arrays of size [1 x #ROIs], and each of these then holds the binding data.
%While this is a very nice organization, we want to concatenate the ROIS
%between runs for file output. ie, 1 file = ROIA from each each run.
Binder_Out = cell(1, TotalROIs);
for i = 1:TotalROIs
    for j = 1:TotalRuns
        Binder_Out{i} = [Binder_Out{i}; Binder{j}{i}];
    end
end




OutIndex = cell(1,TotalROIs); %Preallocating for funsies.
for i = 1:TotalROIs
    OutIndex{i} = zeros(TotalRuns, 20);
    OutIndex{i}(1:TotalRuns,1) = (1:TotalRuns);
    OutIndex{i}(1:TotalRuns,2) = pCa;
    OutIndex{i}(1:TotalRuns,3) = Stats.MFMean(1,:,i);
    OutIndex{i}(1:TotalRuns,4) = Stats.MFVar(1,:,i);
    OutIndex{i}(1:TotalRuns,5) = Stats.AFMean(1,:,i);
    OutIndex{i}(1:TotalRuns,6) = Stats.AFVar(1,:,i);
    OutIndex{i}(1:TotalRuns,7) = Stats.FractXB1Mean(1,:,i);
    OutIndex{i}(1:TotalRuns,8) = Stats.FractXB1Var(1,:,i);
    OutIndex{i}(1:TotalRuns,9) = Stats.FractXB2Mean(1,:,i);
    OutIndex{i}(1:TotalRuns,10) = Stats.FractXB2Var(1,:,i);
    OutIndex{i}(1:TotalRuns,11) = Stats.FractCa0Mean(1,:,i);
    OutIndex{i}(1:TotalRuns,12) = Stats.FractCa0Var(1,:,i);
    OutIndex{i}(1:TotalRuns,13) = Stats.FractCa1Mean(1,:,i);
    OutIndex{i}(1:TotalRuns,14) = Stats.FractCa1Var(1,:,i);
    OutIndex{i}(1:TotalRuns,15) = Stats.FractCa2Mean(1,:,i);
    OutIndex{i}(1:TotalRuns,16) = Stats.FractCa2Var(1,:,i);
    OutIndex{i}(1:TotalRuns,17) = Stats.ATPuseMean(1,:,i);
    OutIndex{i}(1:TotalRuns,18) = Stats.ATPuseVar(1,:,i);
    OutIndex{i}(1:TotalRuns,19) = Stats.FractBoundMean(1,:,i);
    OutIndex{i}(1:TotalRuns,20) = Stats.FractBoundVar(1,:,i);       
end



%% NOW Concat the half times
OutHalfTimes=zeros(TotalRuns, 11);
OutHalfTimes(:, 1)=(1:TotalRuns)';
OutHalfTimes(:, 2)=pCa*ones(TotalRuns, 1);
OutHalfTimes(:, 3)=IndexThalf(1,:)';
OutHalfTimes(:, 4)=IndexThalf(2,:)';
OutHalfTimes(:, 5)=IndexThalf(3,:)';
OutHalfTimes(:, 6)=IndexThalf(4,:)';
OutHalfTimes(:, 7)=IndexThalf(5,:)';
OutHalfTimes(:, 8)=IndexThalf(6,:)';
OutHalfTimes(:, 9)=IndexThalf(7,:)';
OutHalfTimes(:, 10)=IndexThalf(8,:)';
OutHalfTimes(:, 11)=IndexThalf(9,:)';
%%%%%%%%%%%%%%%%%%%% END INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%
% Where every row here is a data input.

TimeSeries_DataOut=[t; Steps(1,:); Steps(2,:); Steps(8,:); Steps(3,:)+Steps(4,:); Steps(3,:); Steps(4,:); Steps(5,:); Steps(6,:); Steps(7,:); Steps(9,:); -log10(Steps(10,:))]';
% Build output filename
OutFile=sprintf('%sTimeSeriesAvg_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
% Open File
fid=fopen(OutFile, 'wt' );	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Time (sec)\tThick F pN\tThin F (pN\tATP per dt\tFrac Bound\tFract. XB1\tFract. XB2\tActins Ca0\tActins Ca1\tActins Ca2\tHandle Position\tpCa\n');	%header for outfile
%Create the output file format string, based on the size and data type.
%Force the default precision is 10 places wide, with 6 decimal places,
% a floating point number, with a tab delimiter.
% This comes out as %10.6f\t, in the Format String.
FormatString=[];
[~, ColOut]=size(TimeSeries_DataOut);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
% Write Data File and close file
fprintf(fid, FormatString, TimeSeries_DataOut');
fclose(fid);     %close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN INDEXED SS OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
for i = 1:TotalROIs
    OutFile=sprintf('%sTwitchData_pCa_%s_ROI_%s.txt', OutDir, num2str(pCa, '%3.2f'),ROIsAZ(i));
    fid=fopen(OutFile,'w');	%open outfile--tab delimited text
    %Create the output file column header:
    fprintf(fid, 'Run Index\tpCa Value\tThickF(pN)\tVARThick F\tThinF (pN)\tVAR Thin F\tFract. XB1\tVAR F XB1\tFract. XB2\tVAR F XB2\tActins Ca0\tVAR ActCa0\tActins Ca1\tVAR ActCa1\tActins Ca2\tVAR ActCa2\tATP per dt\tVAR ATP dt\tFrct Bound\tVAR XB Bnd\n');	%header for outfile
    FormatString=[];
    [~, ColOut]=size(OutIndex{i});
    for j=1:ColOut-1 %for all but last
        FormatString=[FormatString, '%10.6f\t'];
    end
    FormatString=[FormatString, '%10.6f\n'];
    % Write Data File and close file
    fprintf(fid, FormatString, OutIndex{i}');
    fclose(fid);     %close the file
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXED SS OUTPUT TO FILE %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN INDEXED Time Half OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
OutFile=sprintf('%sHalfTimeData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Run Index\tpCa Value\tThickF(pN)\tThinF (pN)\tFract. XB1\tFract. XB2\tActins Ca0\tActins Ca1\tActins Ca2\tATP per dt\tFrct Bound\n');	%header for outfile
FormatString=[];
[~, ColOut]=size(OutHalfTimes);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
% Write Data File and close file
fprintf(fid, FormatString, OutHalfTimes');
fclose(fid);     %close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXED Time Half OUTPUT TO FILE %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
% Only do so if we had binding events
for i = 1:TotalROIs
    %if not(isempty(Binder_Out{i}))
        OutFile=sprintf('%sXBBindingData_pCa_%s_ROI_%s.txt', OutDir, num2str(pCa, '%3.2f'),ROIsAZ(i));
        fid=fopen(OutFile,'w');	%open outfile--tab delimited text
        %Create the output file column header:
        fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
        FormatString=[];
        ColOut=size(Binder_Out{i}, 2); % not sure if valid
        for j=1:ColOut-1 %for all but last
            FormatString=[FormatString, '%10.6f\t'];
        end
        FormatString=[FormatString, '%10.6f\n'];
        % Write Data File and close file
        fprintf(fid, FormatString, Binder_Out{i}');
        fclose(fid);     %close the file
    %end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%BEGIN Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% % Build output filename, open file, write header, write data, close file
% % Just as above:
% % Only do so if we had binding events
% if not(isempty(Binder_All.A))
%     OutFile=sprintf('%sXBBindingData_pCa_%s_ROI_A.txt', OutDir, num2str(pCa, '%3.2f'));
%     fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%     %Create the output file column header:
%     fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
%     FormatString=[];
%     ColOut=size(Binder_All.A, 2); % not sure if valid
%     for i=1:ColOut-1 %for all but last
%         FormatString=[FormatString, '%10.6f\t'];
%     end
%     FormatString=[FormatString, '%10.6f\n'];
%     % Write Data File and close file
%     fprintf(fid, FormatString, Binder_All.A');
%     fclose(fid);     %close the file
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%BEGIN Binding Events OUTPUT TO FILE (for binding events during length change) %%%%%%%%%%%%%%%%%%
% % Build output filename, open file, write header, write data, close file
% % Just as above:
% % Only do so if we had binding events
% if length(fieldnames(Stats)) > 1 && not(isempty(Binder_All.B))
%     OutFile=sprintf('%sXBBindingData_pCa_%s_ROI_B.txt', OutDir, num2str(pCa, '%3.2f'));
%     fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%     %Create the output file column header:
%     fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
%     FormatString=[];
%     ColOut=size(Binder_All.B, 2); % not sure if valid
%     for i=1:ColOut-1 %for all but last
%         FormatString=[FormatString, '%10.6f\t'];
%     end
%     FormatString=[FormatString, '%10.6f\n'];
%     % Write Data File and close file
%     fprintf(fid, FormatString, Binder_All.B');
%     fclose(fid);     %close the file
% end
% if length(fieldnames(Stats)) > 2 && not(isempty(Binder_All.C))
%     OutFile=sprintf('%sXBBindingData_pCa_%s_ROI_C.txt', OutDir, num2str(pCa, '%3.2f'));
%     fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%     %Create the output file column header:
%     fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
%     FormatString=[];
%     ColOut=size(Binder_All.C, 2); % not sure if valid
%     for i=1:ColOut-1 %for all but last
%         FormatString=[FormatString, '%10.6f\t'];
%     end
%     FormatString=[FormatString, '%10.6f\n'];
%     % Write Data File and close file
%     fprintf(fid, FormatString, Binder_All.C');
%     fclose(fid);     %close the file
% end
% if length(fieldnames(Stats)) > 3 && not(isempty(Binder_All.D))
%     OutFile=sprintf('%sXBBindingData_pCa_%s_ROI_D.txt', OutDir, num2str(pCa, '%3.2f'));
%     fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%     %Create the output file column header:
%     fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
%     FormatString=[];
%     ColOut=size(Binder_All.D, 2); % not sure if valid
%     for i=1:ColOut-1 %for all but last
%         FormatString=[FormatString, '%10.6f\t'];
%     end
%     FormatString=[FormatString, '%10.6f\n'];
%     % Write Data File and close file
%     fprintf(fid, FormatString, Binder_All.D');
%     fclose(fid);     %close the file
% end
% if length(fieldnames(Stats)) > 4 && not(isempty(Binder_All.E))
%     OutFile=sprintf('%sXBBindingData_pCa_%s_ROI_E.txt', OutDir, num2str(pCa, '%3.2f'));
%     fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%     %Create the output file column header:
%     fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
%     FormatString=[];
%     ColOut=size(Binder_All.E, 2); % not sure if valid
%     for i=1:ColOut-1 %for all but last
%         FormatString=[FormatString, '%10.6f\t'];
%     end
%     FormatString=[FormatString, '%10.6f\n'];
%     % Write Data File and close file
%     fprintf(fid, FormatString, Binder_All.E');
%     fclose(fid);     %close the file
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
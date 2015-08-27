function [MFMean, AFMean, FractXB1Mean, FractXB2Mean, FractCa0Mean, FractCa1Mean, FractCa2Mean, ATPuseMean, FractBoundMean, MFVar, AFVar, FractXB1Var, FractXB2Var, FractCa0Var, FractCa1Var, FractCa2Var, ATPuseVar, FractBoundVar] = Run_Stats(iTotRun, ROI_index, TempFractCa0, TempFractCa1, TempFractCa2, TempFractXB1, TempFractXB2, TempMFvec, TempAFvec, TempATPuse)


for i = 1:length(ROI_index)
    MFMean(1,iTotRun,i) = mean(TempMFvec(ROI_index{i}(1)));
    AFMean(1,iTotRun,i) = mean(TempAFvec(ROI_index{i}(1)));
    FractXB1Mean(1,iTotRun,i) = mean(TempFractXB1(ROI_index{i}(1)));
    FractXB2Mean(1,iTotRun,i) = mean(TempFractXB2(ROI_index{i}(1)));
    FractCa0Mean(1,iTotRun,i) = mean(TempFractCa0(ROI_index{i}(1)));
    FractCa1Mean(1,iTotRun,i) = mean(TempFractCa1(ROI_index{i}(1)));
    FractCa2Mean(1,iTotRun,i) = mean(TempFractCa2(ROI_index{i}(1)));
    ATPuseMean(1,iTotRun,i) = mean(TempATPuse(ROI_index{i}(1)));
    FractBoundMean(1,iTotRun,i) = FractXB1Mean(1,iTotRun,i) + FractXB2Mean(1,iTotRun,i); % sums: convenient
    
    MFVar(1,iTotRun,i)=var(TempMFvec(ROI_index{i}));
    AFVar(1,iTotRun,i)=var(TempAFvec(ROI_index{i}));
    FractXB1Var(1,iTotRun,i)=var(TempFractXB1(ROI_index{i}));
    FractXB2Var(1,iTotRun,i)=var(TempFractXB2(ROI_index{i}));
    FractCa0Var(1,iTotRun,i)=var(TempFractCa0(ROI_index{i}));
    FractCa1Var(1,iTotRun,i)=var(TempFractCa1(ROI_index{i}));
    FractCa2Var(1,iTotRun,i)=var(TempFractCa2(ROI_index{i}));
    ATPuseVar(1,iTotRun,i)=var(TempATPuse(ROI_index{i}));
    FractBoundVar(1,iTotRun,i)=var(TempFractXB1(ROI_index{i}) + TempFractXB2(ROI_index{i})); % these are correlated probably, so we can't just sum variances
end




% MFMean(iTotRun) = mean(TempMFvec(ROI_index));
% AFMean(iTotRun) = mean(TempAFvec(ROI_index));
% FractXB1Mean(iTotRun) = mean(TempFractXB1(ROI_index));
% FractXB2Mean(iTotRun) = mean(TempFractXB2(ROI_index));
% FractCa0Mean(iTotRun) = mean(TempFractCa0(ROI_index));
% FractCa1Mean(iTotRun) = mean(TempFractCa1(ROI_index));
% FractCa2Mean(iTotRun) = mean(TempFractCa2(ROI_index));
% ATPuseMean(iTotRun) = mean(TempATPuse(ROI_index));
% FractBoundMean(iTotRun) = FractXB1Mean(iTotRun) + FractXB2Mean(iTotRun); % sums: convenient
% 
% MFVar(iTotRun)=var(TempMFvec(ROI_index), 1);
% AFVar(iTotRun)=var(TempAFvec(ROI_index), 1);
% FractXB1Var(iTotRun)=var(TempFractXB1(ROI_index), 1);
% FractXB2Var(iTotRun)=var(TempFractXB2(ROI_index), 1);
% FractCa0Var(iTotRun)=var(TempFractCa0(ROI_index), 1);
% FractCa1Var(iTotRun)=var(TempFractCa1(ROI_index), 1);
% FractCa2Var(iTotRun)=var(TempFractCa2(ROI_index), 1);
% ATPuseVar(iTotRun)=var(TempATPuse(ROI_index), 1);
% FractBoundVar(iTotRun)=var(TempFractXB1(ROI_index) + TempFractXB2(ROI_index), 1); % these are correlated probably, so we can't just sum variances

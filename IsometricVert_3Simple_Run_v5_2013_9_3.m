%% BCWT, August 17, 2010-- Hack the Isometric Skeltal code:
%% Include new Geomtry--should be in the CaLookup File (4-start)
%% Prepare Istometric code for David to run some isometric force-pCa
%% simulations.  Make this more of a function use more initialization as a
%% script too.
%% Changes to this code:
%% 1) Force isometric, and only one time-half
%% 2) Functionalize
%% 3) Comment Out/Delete 2nd half averages
%% 4) Add in variances as output
%% 5) For ktr rate; remove the data fitting to an Exp and simply do the
%%    t1/2 from the SS isometric rate
%% 6) Bring into simulation runs a toggle for XBKO fraction.

%% Function Inputs:

%% Function Outputs:
% OutIndex%%(24, TotalRuns)
% OutExpFitInfo%%(28,1)
% OutSimpAve%%(20, 1)

%% Dependent scripts:


function [TimeSeries_DataOut, OutIndex, OutHalfTimes] = IsometricVert_3Simple_Run_v5_2013_9_3(OutDir, Psi, TnFraction, XB_Fraction, Tm_Type, pCa, Run_bool, SXB1, SXB2, SXB3, Coop_Type, ScaleK1, ScaleK2, ScaleK3)

% Build our lattice, stiffness values
IsometricVert_3Simple_Lattice_v1
% Load all the constants, coefficients for XB cycling
% paramters and Free energies and such.  Via script:
IsometricVert_3Simple_XB_pCa_Initialize_v1
IsometricVert_3Simple_Output_Variables_v1
IsometricVert_3Simple_CaLookup


% Sparse up or big matrix
Mcore=sparse(Mcore);

%% Ca should have a single value
iCa=1;

%% We will need to zero out these values to keep the
%% temparary, so that we can assign them at the end
%% of a Ca Loop.  We only need to keep track of the Averages
%% because the Standard Deviations can be taken later.
kCa_Ave_TEMP=zeros(TotalRuns, 1);
SS_FCa_Ave_TEMP=zeros(TotalRuns, 1);
kTr_Ave_TEMP=zeros(TotalRuns, 1);
SS_FTr_Ave_TEMP=zeros(TotalRuns, 1);
XB1_kCa_Ave_TEMP=zeros(TotalRuns, 1);
SS_XB1Ca_Ave_TEMP=zeros(TotalRuns, 1);
XB1_kTr_Ave_TEMP=zeros(TotalRuns, 1);
SS_XB1Tr_Ave_TEMP=zeros(TotalRuns, 1);
XB2_kCa_Ave_TEMP=zeros(TotalRuns, 1);
SS_XB2Ca_Ave_TEMP=zeros(TotalRuns, 1);
XB2_kTr_Ave_TEMP=zeros(TotalRuns, 1);
SS_XB2Tr_Ave_TEMP=zeros(TotalRuns, 1);
FA_kCa_Ave_TEMP=zeros(TotalRuns, 1);
SS_FACa_Ave_TEMP=zeros(TotalRuns, 1);

%% Big Binder
Big_Binder_3=[];

for iTotRun=1:TotalRuns
    
%     %% MOD HERE FOR TnKO and XBKO
%         %% Here is for a random TnKO
%         %% Initialize the matrix
%         TnKO=rand(NACT, NTn);
%         [iTest, jTest]=find(TnKO<=TnFraction);
%         %% zero out and update the matix
%         TnKO=zeros(NACT, NTn);
%         for iii=1:length(iTest)
%             %% Now the matrix is zeros, update with the Tn that will actually
%             %% exist. As a 1;
%             TnKO(iTest(iii,1), jTest(iii,1) )=1;
%         end
%     
    
%% MOD September 3, 2013
    TnKO=ones(NACT, NTn); %% Fill as though no KOs
    [iTest, jTest]=find(TnKO==1); %% Gather row and col indicies to
    %% filter on the KO
    TnKO=zeros(NACT, NTn);
    for iii=1:length(jTest) %% Loop down each row, looking at first half
        %% then the second half to fill with zero if NTn<(1-TnKO)
        %% Now the matrix is zeros, update with the Tn that will actually
        %% exist.
        if(iii<=(NTn/2)) %% First Half
            %if( (iii/(NTn/2))<(1-TnFraction) ) %% ERROR--9/3/2013
            if( (iii/(NTn/2))<=(TnFraction) )
                TnKO(:, iii)=1;
            end
        end
        if(iii>(NTn/2)) %% 2nd Half
            %if( ((iii-NTn/2)/(NTn/2))<(1-TnFraction) ) %% ERROR--9/3/2013
            if( ((iii-NTn/2)/(NTn/2))<=(TnFraction) )
                TnKO(:, iii)=1;
            end
        end
    end
    
    
    %% Now repeat for XBKO--for a random KO
    %size(LAUREL)
    XBKO=rand(size(LAUREL)); %% Get here the XBKO array--similar to a boolean with same size as LAUREL
    [iTest, jTest]=find(XBKO<=XB_Fraction);
    XBKO=zeros(size(LAUREL)); %% Reset to fill only the 1=Head Available
    for iii=1:length(iTest)
        XBKO(iTest(iii,1), jTest(iii,1) )=1;
    end
    
    
    %     %% Now repeat for XBKO for KO from the Free End
    %     XBKO=ones(size(LAUREL)); %% Fill as though no KOs
    %     [iTest, jTest]=find(XBKO==1); %% Gather row and col indicies to filter the KO on
    %     for iii=1:length(iTest)
    %         %% If we are XB_Fraction near the free end of the filament
    %         %% Set XBKO to 0.  This will also set XBK) to zeros
    %         %% for the 61st XB Node (handle) at the M-line
    %         if(mod(iii,NumBridges)<=NumBridges*(1-XB_Fraction))
    %             XBKO(iTest(iii,1), jTest(iii,1) )=0;
    %         end
    %     end
    
    
    CountOff=0; % To count the total number of reverse transitions from
    % unbound to strongly bound in the XB cycle.
    
    TFstate=zeros(ACTNodes, 3);  % Refresh initial value
    TFstate(1:NumSites:ACTNodes, :)=-1; % To set the handle nodes to -1, so as
    % to not interfere with thinking they are counted as an unboud state.
    Prior_TFstate=[];
    Binder_3=[];
    M = Mcore; % Refresh
    EndVec = EndVeccore; % Freshen up?
    X=zeros(TOTNodes, 1); %Because we need to know position.
    
    
    % Set up our temperary variables/zero them out.
    TempFractCa0=zeros(1, NSTEPS);
    TempFractCa1=zeros(1, NSTEPS);
    TempFractCa2=zeros(1, NSTEPS);
    TempFractXB1=zeros(1, NSTEPS);
    TempFractXB2=zeros(1, NSTEPS);
    TempMFvec=zeros(1, NSTEPS);
    TempAFvec=zeros(1, NSTEPS);
    TempATPuse=zeros(1, NSTEPS);
    ATPcycle=0;
    
    for iStep = 1:NSTEPS
        % Take in to account about our state of where things are and
        % Load in our statistics for what we need to save for the final data.
        %mprime = sparse(M);
        mprime = M;
        X = mprime\EndVec;
        % Must comment out the BigX, we don't have the memory without
        % more page file space.
        %BigX(:,iStep, iCa)=X;  %Save our position vector at each timestep.
        EndVec = EndVeccore;
        M = Mcore;
        
        % Begin here to load up our Temp vecters, that will be dumped at the end of
        % each run.  Load the temp each timestep.
        TempFractXB1(1, iStep)=length(find(TFstate(:,2)==1))/(N_Thick_Start*(MYONodes-NMYO)); % Fract XB Low force bearing
        TempFractXB2(1, iStep)=length(find(TFstate(:,2)==2))/(N_Thick_Start*(MYONodes-NMYO)); % Fract XB High force bearing
        TempFractCa0(1, iStep)=length(find(TFstate(:,1)==0))/(ACTNodes-NACT);
        TempFractCa1(1, iStep)=length(find(TFstate(:,1)==1))/(ACTNodes-NACT);
        TempFractCa2(1, iStep)=length(find(TFstate(:,1)==2))/(ACTNodes-NACT); % Fract Avail BS.
        
        
        %         % Force Calculation on Myosin and Actin
        % actin Forces (near handel end)  ( Handle+1 - Handle )
        TempAFvec(1,iStep)= ((X(2,1)-X(1,1)) - LACT)*KACT ...
            + ((X(1*NumSites+2,1)-X(1*NumSites+1,1)) - LACT)*KACT ...
            + ((X(2*NumSites+2,1)-X(2*NumSites+1,1)) - LACT)*KACT ...
            + ((X(3*NumSites+2,1)-X(3*NumSites+1,1)) - LACT)*KACT ...
            + ((X(4*NumSites+2,1)-X(4*NumSites+1,1)) - LACT)*KACT ...
            + ((X(5*NumSites+2,1)-X(5*NumSites+1,1)) - LACT)*KACT ...
            + ((X(6*NumSites+2,1)-X(6*NumSites+1,1)) - LACT)*KACT ...
            + ((X(7*NumSites+2,1)-X(7*NumSites+1,1)) - LACT)*KACT;
        
        % myosin Forces
        TempMFvec(1,iStep)=( (X(1*NumBridges + ACTNodes, 1) - X(1*NumBridges + (ACTNodes-1), 1) )-LMYO)*KMYO ...
            + ( (X(2*NumBridges + ACTNodes, 1) - X(2*NumBridges + (ACTNodes-1), 1) )-LMYO)*KMYO ...
            + ( (X(3*NumBridges + ACTNodes, 1) - X(3*NumBridges + (ACTNodes-1), 1) )-LMYO)*KMYO ...
            + ( (X(4*NumBridges + ACTNodes, 1) - X(4*NumBridges + (ACTNodes-1), 1) )-LMYO)*KMYO;
        
        % Load the average ATPuse;
        ATPuse(1, iStep)=ATPuse(1, iStep) + ATPcycle/TotalRuns;
        TempATPuse(1, iStep)=ATPcycle;
        ATPcycle=0;
        
        %%%********************************************************
        %% If we want to keep track of time-on; binding
        %% Only keep track of binding over the final portion of the
        %% simulation; where we also will average the steady-state
        %% information
        if iStep>=(1-DataSpan)*NSTEPS
            [Prior_TFstate, Binder_3] = Time_On_Calculator(X, TFstate, Prior_TFstate, Binder_3, (1-DataSpan)*NSTEPS, iStep, reach, kxscaler);
        end
        %%%********************************************************
        
        
        %%%*****************************************************************
        %This is the loop over the actin nodes to check and see their
        %state, and assing them Ca binding or not.  Regulation loop.
        
        switch Tm_Type
            case 0 %% 18 nm
                TFstate=CaRegCoop_0a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
            case 1
                TFstate=CaRegCoop_1a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
            case 2
                TFstate=CaRegCoop_2a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
            case 3
                TFstate=CaRegCoop_3a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
            case 4 %% 74 nm
                TFstate=CaRegCoop_3b_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
        end
        %%%*****************************************************************
        
        
        
        %%%*****************************************************************
        %%% Now that we are RU state, check about the XB Binding
        
        occupied=0; %this will set occupied to 0, and a 1 by 1 matrix.
        
        for iXB = 1:MYONodes
            for iStart=1:N_Thick_Start
                % loop over every myosin node--looking across the N_Thick_Start
                % Because LAUREL is no-longer a vector, we will have to shuffel
                % this to work with an Array--which will work out just fine now
                if LAUREL(iXB, iStart) % 0 implies handle node
                    % or no XB pointing to a thick filament so skip those.
                    
                    % LAUREL(iXB, iStart) is a number between 1-24--pointing to
                    % the row of HARDY that could be co-linearly aligned
                    % The XB which is is facing has been stored in the
                    % TFstate col. 3 as a number iXB+ACTNodes.
                    
                    
                    %% Check whether this XB is knocked out or now
                    %% If working = 1; if Knocked Out = 0
                    if XBKO(iXB, iStart)==1
                        
                        % Where the variable 'node' for the rest of this
                        % loop refers to the Actin Node for a given iXB
                        % might be looking at.
                        
                        % Check here to see whether our XB head is indeed bound
                        % Pull back the node (thin filament node we are looking for)
                        node=find(TFstate( (HARDY(LAUREL(iXB,iStart),:) ) , 3)==(iXB+ACTNodes));
                        % This will be an index to use in HARDY, returning
                        % the act node to look at.
                        if isempty(node)
                            % If it has returned an empty set, then we will need to
                            % search over the binding sites and find the nearest
                            % node for it to bind.
                            
                            % So here--we take the difference between the XB position and all the
                            % potential co-linear sites for our thin filament target zones
                            % where dist is the value of the minimum difference
                            % where node is the index used just below--a bit unclean but oh well
                            [dist, node] = min(abs(X(iXB+ACTNodes) - X(HARDY(LAUREL(iXB,iStart),:))));
                            % Set the node we are looking as row in TFstate
                            % here node is the thin filament node
                            node=HARDY(LAUREL(iXB,iStart), node);
                            % Set the distance as a vector, rather than abs
                            % value as above.
                            dist=X(iXB+ACTNodes,1) - X(node,1);
                        else
                            % it will already be bound somewhere.  Now we
                            % just need to set the node we are looking at
                            % in the row of TFstate
                            node=HARDY(LAUREL(iXB,iStart), node);
                            dist=X(iXB+ACTNodes,1) - X(node,1);
                        end
                        
                        
                        % Now we are all set with our node we are looking at and we
                        % have calulated the distance, then we
                        % can proceed with looking at the binding for a given iXB.
                        % The XB state it kept in the TFstate col. 2, and the
                        % XB it is bound with is kept in the TFstate col. 3.
                        % Case 0 --Unbound XB
                        % Case 1 --Stongly Bound XB, with low force into
                        % the matrix, such that the energy from ATP
                        % hydrolysis remains in the XB head.
                        % Case 2 --Stronly Bound XB, with high force into
                        % the latice, such that the E from ATP hydrolys is
                        % transfered from the myosin molecule into the
                        % lattice and there is a displacement.
                        %
                        % The XB dynamics are set up and defined in the
                        % mexDooItRates.m program.
                        % All Cases--Check for the Regulation TFstate col. 1 to
                        % be in state 2.
                        
                        % If the node has already been attached, then another XB cannot look at
                        % it as well.  So, check this by looking at whether occupied is
                        % empty?
                        if ( (isempty(find(occupied==node))) ) %& ( (dist > -12)&(dist<12) ) );
                            % The node has not been occupied if we get this
                            % far.  However, we should keep track that is
                            % the nod gets occupied, we will fill the free
                            % end of the occupied vector using: occupied(end+1,1)=node;
                            
                            % clear p01 p02 p12 p10 p20 p21 %PERFMOD
                            if (TFstate(node, 1)==2 & TFstate(node, 2)==0)  %%%%%%%% AVAIL. and UNBOUND
                                % computer forward and reverse transition likelihoods
                                %  [k01,k10,k12,k21,k20,k02] = mexDooItRates(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD);
                                [k01,k10,k12,k21,k20,k02] = IsometricVert_3Simple_XBRates3(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD, SXB1, SXB2, SXB3);
                                p01 = dt*k01;
                                p02 = dt*k02;
                                gx = rand;
                                
                                if (gx <= p01) % It goes forward to weak state
                                    % Update the TFState in the XB state col. 2
                                    % and in the head/site match, col. 3.
                                    TFstate(node, 2)=1;
                                    TFstate(node, 3)=(iXB+ACTNodes); %% Load thick filament node bound at target zone
                                    occupied(end+1,1)=node;  % The node has not been occupied
                                elseif (gx >= 1.0 - p02) % Reverse rarely happens.
                                    TFstate(node, 2)=2;
                                    TFstate(node, 3)=(iXB+ACTNodes);
                                    CountOff=CountOff+1;
                                    occupied(end+1,1)=node;  % The node has not been occupied
                                    
                                end
                            elseif (TFstate(node, 1)==2 & TFstate(node, 2)==1)  %%%%%%%% WEAK BOUND and AVAIL.
                                % computer forward and reverse transition likelihoods
                                %  [k01,k10,k12,k21,k20,k02] = mexDooItRates(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD);
                                [k01,k10,k12,k21,k20,k02] = IsometricVert_3Simple_XBRates3(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD, SXB1, SXB2, SXB3);
                                p12 = dt*k12;
                                p10 = dt*k10;
                                gx = rand;
                                
                                if (gx <= p12) % It goes forward to strong state
                                    % Update the TFState in the XB state col. 2
                                    % and in the head/site match, col. 3.
                                    TFstate(node, 2)=2;
                                    TFstate(node, 3)=(iXB+ACTNodes);
                                    occupied(end+1,1)=node;  % The node has not been occupied
                                elseif (gx >= 1.0 - p10) % Reverse rarely happens.
                                    TFstate(node, 2)=0;
                                    TFstate(node, 3)=0;
                                end
                            elseif (TFstate(node, 1)==2 & TFstate(node, 2)==2)  %%%%%%%% STRONG BOUND and AVAIL.
                                % computer forward and reverse transition likelihoods
                                %  [k01,k10,k12,k21,k20,k02] = mexDooItRates(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD);
                                [k01,k10,k12,k21,k20,k02] = IsometricVert_3Simple_XBRates3(dist, reach, kRT, f1, eta, dGtot, AA, BB, CC, DD, SXB1, SXB2, SXB3);
                                p20 = dt*k20;
                                p21 = dt*k21;
                                gx = rand;
                                
                                if (gx <= p20) % It goes forward to unbound state
                                    % Update the TFState in the XB state col. 2
                                    % and in the head/site match, col. 3.
                                    %CountCycle=CountCycle+1;
                                    ATPcycle=ATPcycle+1;
                                    TFstate(node, 2)=0;
                                    TFstate(node, 3)=0;
                                elseif (gx >= 1.0 - p21) % Reverse rarely happens.
                                    TFstate(node, 2)=1;
                                    TFstate(node, 3)=(iXB+ACTNodes);
                                    occupied(end+1,1)=node;  % The node has not been occupied
                                end
                            end
                            
                        end
                    end % End Loop on Wheather XB is KNocked Out
                    
                end % if over the 'whether' handel? (If over Hardy(iXB==0) to over myosin handles.
            end    % end loop over iStart
        end    % end loop over Myo. Heads
        
        %%%*****************************************************************
        
        % this block modifies the matrix which encodes the connections
        % between myosin and actin nodes along XB. (M)
        % For any force bearing XB/ACTIN pair (state 2==strong bound high force, and state 1
        % == strongly bound, with low force).
        % The matrix will be loaded.
        % The reason this is done outside of the above loop is to make
        % sure there is not superload of the matrix if I made any counting
        % mistake
        temp=find(TFstate(:,2)>0);
        
        for i=1:length(temp)
            %         if i==1
            %             disp([temp, TFstate(temp, 2), TFstate(temp, 3), BigBoundState(temp, (iStep)), BigBoundPair(temp, (iStep))])
            %         end
            % NOTE on what we have.
            %temp_site=temp(i);
            %temp_brdg=(TFstate(temp(i),3));
            % Fill the actin diagonal
            M(temp(i),temp(i)) = M(temp(i),temp(i)) + kxscaler;
            % Fill the myosin diagonal
            M(TFstate(temp(i),3),TFstate(temp(i),3)) = M(TFstate(temp(i),3),TFstate(temp(i),3)) + kxscaler;
            % Fill the cross term on the actin row
            M(temp(i),TFstate(temp(i),3)) = M(temp(i),TFstate(temp(i),3)) - kxscaler;
            % Fill the cross term on the myosin row
            M(TFstate(temp(i),3),temp(i)) = M(TFstate(temp(i),3),temp(i)) - kxscaler;
        end
        
        
    end % end of one run (iStep)
    
    % Output from a given run.
    
    % Update the fitting in the matrix.
    % Collapse the Temp. vectors into the big vector that will be used
    % for final averages:
    FractXB1(1, 1:NSTEPS)=FractXB1(1, 1:NSTEPS)+TempFractXB1/TotalRuns;
    FractXB2(1, 1:NSTEPS)=FractXB2(1, 1:NSTEPS)+TempFractXB2/TotalRuns;
    FractCa0(1, 1:NSTEPS)=FractCa0(1, 1:NSTEPS)+TempFractCa0/TotalRuns;
    FractCa1(1, 1:NSTEPS)=FractCa1(1, 1:NSTEPS)+TempFractCa1/TotalRuns;
    FractCa2(1, 1:NSTEPS)=FractCa2(1, 1:NSTEPS)+TempFractCa2/TotalRuns;
    MFvec(1, 1:NSTEPS)=MFvec(1, 1:NSTEPS)+TempMFvec/TotalRuns;
    AFvec(1, 1:NSTEPS)=AFvec(1, 1:NSTEPS)+TempAFvec/TotalRuns;
    
    %% Load All Our Index Variables
    %%%%%%%% Begin Index for the First Half of the Simulation
    %% First, Intermediate Indexing for the Averages:
    IndexSSMFvec(1, iTotRun)=mean(TempMFvec(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSAFvec(1, iTotRun)=mean(TempAFvec(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSFractXB1(1, iTotRun)=mean(TempFractXB1(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSFractXB2(1, iTotRun)=mean(TempFractXB2(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSFractCa0(1, iTotRun)=mean(TempFractCa0(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSFractCa1(1, iTotRun)=mean(TempFractCa1(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSFractCa2(1, iTotRun)=mean(TempFractCa2(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSATPuse(1, iTotRun)=mean(TempATPuse(1, (1-DataSpan)*NSTEPS:NSTEPS));
    IndexSSFractBound(1, iTotRun)=mean( TempFractXB1(1, (1-DataSpan)*NSTEPS:NSTEPS) + TempFractXB2(1, (1-DataSpan)*NSTEPS:NSTEPS) );
    % Continues with Intermediate Indexing for the
    % variance about the mean for the index forming the Averages:
    IndexVARMFvec(1, iTotRun)=var(TempMFvec(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARAFvec(1, iTotRun)=var(TempAFvec(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARFractXB1(1, iTotRun)=var(TempFractXB1(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARFractXB2(1, iTotRun)=var(TempFractXB2(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARFractCa0(1, iTotRun)=var(TempFractCa0(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARFractCa1(1, iTotRun)=var(TempFractCa1(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARFractCa2(1, iTotRun)=var(TempFractCa2(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARATPuse(1, iTotRun)=var(TempATPuse(1, (1-DataSpan)*NSTEPS:NSTEPS), 1);
    IndexVARFractBound(1, iTotRun)=var( TempFractXB1(1, (1-DataSpan)*NSTEPS:NSTEPS) + TempFractXB2(1, (1-DataSpan)*NSTEPS:NSTEPS) , 1);
    
    %% Also calculate the time to reach 1/2 of the SS value for each of the
    %% variables.
    %% Load them with time value, closest to half of the max avg value
    [val,Index] = min( abs( TempMFvec-(0.5*IndexSSMFvec(1, iTotRun)) ) );
    IndexThalfMFvec(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( TempAFvec-(0.5*IndexSSAFvec(1, iTotRun)) ) );
    IndexThalfAFvec(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( TempFractXB1-(0.5*IndexSSFractXB1(1, iTotRun)) ) );
    IndexThalfFractXB1(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( TempFractXB2-(0.5*IndexSSFractXB2(1, iTotRun)) ) );
    IndexThalfFractXB2(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( TempFractCa0-(0.5*IndexSSFractCa0(1, iTotRun)) ) );
    IndexThalfFractCa0(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( TempFractCa1-(0.5*IndexSSFractCa1(1, iTotRun)) ) );
    IndexThalfFractCa1(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( TempFractCa2-(0.5*IndexSSFractCa2(1, iTotRun)) ) );
    IndexThalfFractCa2(1, iTotRun)=dt*Index;
    [val,Index] = min( abs(TempATPuse-(0.5*IndexSSATPuse(1, iTotRun)) ) );
    IndexThalfATPuse(1, iTotRun)=dt*Index;
    [val,Index] = min( abs( (TempFractXB1(1, (1-DataSpan)*NSTEPS:NSTEPS) + TempFractXB2(1, (1-DataSpan)*NSTEPS:NSTEPS))- (0.5*IndexSSFractBound(1, iTotRun)) ) );
    IndexThalfFractBound(1, iTotRun)=dt*Index;
    
    
    %%%%%%%% Begin Index for the Second Half of the Simulation Outputs--
    Big_Binder_3=[Big_Binder_3; Binder_3];
    
end % End loop over the runs for a Ca level.

%%%%%%%%%%%%%%%%%%%% BEGIN INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%
% Package Up the Averaged (Index) information into a vector.
% Also, flip to column vectors, so we can write it to a outfile
OutIndex=zeros(TotalRuns, 20);
%% 24 rows for the first half results
OutIndex(:, 1)=[1:TotalRuns]';
OutIndex(:, 2)=pCa*ones(TotalRuns, 1);
%% Concat the SS mean data, followed by Variance
OutIndex(:, 3)=IndexSSMFvec';
OutIndex(:, 4)=IndexVARMFvec';
OutIndex(:, 5)=IndexSSAFvec';
OutIndex(:, 6)=IndexVARAFvec';
OutIndex(:, 7)=IndexSSFractXB1';
OutIndex(:, 8)=IndexVARFractXB1';
OutIndex(:, 9)=IndexSSFractXB2';
OutIndex(:, 10)=IndexVARFractXB2';
OutIndex(:, 11)=IndexSSFractCa0';
OutIndex(:, 12)=IndexVARFractCa0';
OutIndex(:, 13)=IndexSSFractCa1';
OutIndex(:, 14)=IndexVARFractCa1';
OutIndex(:, 15)=IndexSSFractCa2';
OutIndex(:, 16)=IndexVARFractCa2';
OutIndex(:, 17)=IndexSSATPuse';
OutIndex(:, 18)=IndexVARATPuse';
OutIndex(:, 19)=IndexSSFractBound';
OutIndex(:, 20)=IndexVARFractBound';


%% NOW Concat the half times
OutHalfTimes=zeros(TotalRuns, 11);
OutHalfTimes(:, 1)=[1:TotalRuns]';
OutHalfTimes(:, 2)=pCa*ones(TotalRuns, 1);
OutHalfTimes(:, 3)=IndexThalfMFvec';
OutHalfTimes(:, 4)=IndexThalfAFvec';
OutHalfTimes(:, 5)=IndexThalfFractXB1';
OutHalfTimes(:, 6)=IndexThalfFractXB2';
OutHalfTimes(:, 7)=IndexThalfFractCa0';
OutHalfTimes(:, 8)=IndexThalfFractCa1';
OutHalfTimes(:, 9)=IndexThalfFractCa2';
OutHalfTimes(:, 10)=IndexThalfATPuse';
OutHalfTimes(:, 11)=IndexThalfFractBound';
%%%%%%%%%%%%%%%%%%%% END INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%
% Where every row here is a data input.

TimeSeries_DataOut=[...
    t(1, 1:NSTEPS); MFvec(1, 1:NSTEPS); AFvec(1, 1:NSTEPS);...
    ATPuse(1, 1:NSTEPS); FractXB1(1, 1:NSTEPS)+FractXB2(1, 1:NSTEPS); ...
    FractXB1(1, 1:NSTEPS); FractXB2(1, 1:NSTEPS); ...
    FractCa0(1, 1:NSTEPS); FractCa1(1, 1:NSTEPS); FractCa2(1, 1:NSTEPS)]';
% Build output filename
OutFile=sprintf('%sTimeSeriesAvg_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
% Open File
fid=fopen(OutFile, 'wt' );	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Time (sec)\tThick F pN\tThin F (pN\tATP per dt\tFrac Bound\tFract. XB1\tFract. XB2\tActins Ca0\tActins Ca1\tActins Ca2\n');	%header for outfile
%Create the output file format string, based on the size and data type.
%Force the default precision is 10 places wide, with 6 decimal places,
% a floating point number, with a tab delimiter.
% This comes out as %10.6f\t, in the Format String.
FormatString=[];
[RowOut, ColOut]=size(TimeSeries_DataOut);
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
OutFile=sprintf('%sSSData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Run Index\tpCa Value\tThickF(pN)\tVARThick F\tThinF (pN)\tVAR Thin F\tFract. XB1\tVAR F XB1\tFract. XB2\tVAR F XB2\tActins Ca0\tVAR ActCa0\tActins Ca1\tVAR ActCa1\tActins Ca2\tVAR ActCa2\tATP per dt\tVAR ATP dt\tFrct Bound\tVAR XB Bnd\n');	%header for outfile
FormatString=[];
[RowOut, ColOut]=size(OutIndex);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
% Write Data File and close file
fprintf(fid, FormatString, OutIndex');
fclose(fid);     %close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXED SS OUTPUT TO FILE %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%BEGIN INDEXED Time Half OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
OutFile=sprintf('%sHalfTimeData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Run Index\tpCa Value\tThickF(pN)\tThinF (pN)\tFract. XB1\tFract. XB2\tActins Ca0\tActins Ca1\tActins Ca2\tATP per dt\tFrct Bound\n');	%header for outfile
FormatString=[];
[RowOut, ColOut]=size(OutHalfTimes);
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
if not(isempty(Big_Binder_3))
    OutFile=sprintf('%sXBBindingData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
    fid=fopen(OutFile,'w');	%open outfile--tab delimited text
    %Create the output file column header:
    fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
    FormatString=[];
    [RowOut, ColOut]=size(Binder_3);
    for i=1:ColOut-1 %for all but last
        FormatString=[FormatString, '%10.6f\t'];
    end
    FormatString=[FormatString, '%10.6f\n'];
    % Write Data File and close file
    fprintf(fid, FormatString, Big_Binder_3');
    fclose(fid);     %close the file
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXED Time Half OUTPUT TO FILE %%%%%%%%%%%%%%%%%%

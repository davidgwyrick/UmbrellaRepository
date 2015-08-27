%% BCWT, March 11, 2008
%% This function will take care of cooperativity for the TmSpan=74 nm
%% which is 2*37 nm.  This function, however, only deals with span and does
%% not sort the cooperative activation any differently than the TmSpan=56
%% nm for the different types of cooperativity.  This will not matter any
%% difference when CoopType==0; where the span seems to show up as
%% physiologically relevant.

function [TFstate]=CaRegCoop_3a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO)

%%% Unpack coopeartativity parameters
TF3TF12=CoopPass(1, 1);
XB2TF12=CoopPass(2, 1);
XB3TF12=CoopPass(3, 1);
XB2TF23=CoopPass(4, 1);
XB3TF23=CoopPass(5, 1);
TF3TF23=CoopPass(6, 1);
XB2TF21=CoopPass(7, 1);
XB3TF21=CoopPass(8, 1);
TF3TF21=CoopPass(9, 1);
XB2TF32=CoopPass(10, 1);
XB3TF32=CoopPass(11, 1);
TF3TF32=CoopPass(12, 1);
CooperativeType=CoopPass(13, 1);

%%% Unpack Thin Filament Rates
koff=TFRatePass(1,1);
kon=TFRatePass(2,1);
RuOff=TFRatePass(3,1);
RuOn=TFRatePass(4,1);
CaOff=TFRatePass(5,1);
CaOn=TFRatePass(6,1);



% Look over each BS
% Look over each BS
% We will house the Ca binding in the first col. of the
% TFstate matrix.  The binding is
% 0==CaFree + TnC + TnI.A
% 1==CaTnC + TnI.A
% 2==CaTnC.TnI + A --> This is the state needed to make
% actin nodes available to bind.
for iRow=1:NACT  % Loop over each filament
  for iTn=1:NTn % Loop down the Tn index of each filament
    %% Now test whether this Tn can be modified.  TnKO will be zero if
    %% there are not Tn avail, and 1 if Tn is avail
    if TnKO(iRow, iTn) == 1
      ii=Tn(iRow, iTn); % This will dump a given Actin node into ii
      if (mod(iTn,(NTn/2))==1) % if iTn==1 or 24
        %Then we are at the beginning near the Z-line and
        %we can only check the next RU down the filament
        %towards the M-line.  Within an RU at the Z-line
        %end we will have Tm
        % and that will be treated a little differently.
        % Each Tn will controll itself (Super Actin) and
        % 1 actin further down the filament from it.
        switch TFstate(ii, 1)
          case 0 %ACTNode off
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF12 rates--Ca2+binding (p01 in prob.)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Not at the free end of the actin filament
                % we have to search at a single
                % neighboring regulatory unit to the left and
                % right of the current RU. The previous
                % Tn site and the next Tn site will be found.
                % Then querried
                % to see whether any binding occurs.
                inextTn=Tn(iRow, iTn+1);
                %% June 25, 2007
                %% Here we are at the z-line end, so we only look
                %% to the right.
                % Now actin inextTn:2:inextTn+2 should
                % represent all the nodes within the
                % two neighboring RUs.
                %% June 25, 2007
                %% Using these indices on the TFState variable will
                %% look at the actin nodes themselves---rather than
                %% the Tn nodes, therefore allowing a query of
                %% what goes on within regulatory units, rather
                %% than the Tn sites.
                %And we can look at either
                % column 1 (thin fil state)  or column 2
                % (XB state) of TFstate to querry
                % cooperativity.--
                if (max(TFstate(inextTn:2:inextTn+2, 1)==2)>0)
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 2  %% 2: XB3RU Only
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0)
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 3    %% 3: Both RURU or XB3RU
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 4 %% 4: XB2RU Only
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0)
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 5 %% 5: Both RURU or XB2RU
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 6 %% 6: Both XB3RU or XB2RU
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                inextTn=Tn(iRow, iTn+1);
                if ( (max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p01 = dt*kon*Ca(iCa);
            end %% End switch CoopType influence pij

            %% p02 is reverse transition, un affected
            %% by cooperativity.
            p02 = dt*CaOn*Ca(iCa);
            gx = rand;  % seed this puppy ? -- seeded at top of program
            % If the rand is within the transition prob.
            % Then move to the second state
            % computer forward and reverse transition likelihoods
            if (gx <= p01)
              TFstate(ii,1)=1;  % The transition to the next state
              % This also will transition the next actin on the particular
              % helix
              %% June 25, 2007
              %% MOD: TFstate(ii+2,1)=1;

              %% MOD: March 11, 2008
              %% For the added span, about what to do at the next TnC
              %% down the pike.  So, if it is one of the TnKO, then it will
              %% be controlled only by the current Tn.  However, if it is
              %% it's individual (not a KO), then the same rule applies as
              %% the TFstate(ii+1, 1) just above.
              if TnKO(iRow, iTn+1) == 0
                %% Then it is a knockout and should be set at 1
                TFstate(ii+4,1)=1;
                TFstate(ii+2,1)=1;
              else
                %% Otherwise, it controls itself and we have to querry
                TFstate(ii+4,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                %% Unless the next Tn over is activated.  And therefore that
                %% Tmspan of that guy will be greater than just iteslf, and
                %% one to the right.  It will cover this guy (one to the left
                %% of that Tn guy) as well.
                %% So querry our State at the next Tn
                %% TFstate(Tn(iRow, iTn+1), 1);
                %% And modify the shared node accordingly
                TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
              end

            elseif (gx >= 1.0 - p02)
              TFstate(ii, 1)=2;  % The reverse transition, which should rarely happen.
              % This also will transition the next actin on the particular
              % helix
              TFstate(ii+2,1)=2;
              %% And the next Tn Over
              TFstate(ii+4,1)=2;
            end

          case 1 %ACTNode is in the intermediate position
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF23 rates--Ru Activation (p12, prob) &
            % TF21 rates--Ca dissociation (p10, prob)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Again, both neighboring RU
                inextTn=Tn(iRow, iTn+1);
                % Again actin inextTn:2:inextTn+2
                % represent all the nodes within the
                % previous RU.--

                if (max(TFstate(inextTn:2:inextTn+2, 1)==2)>0)
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*TF3TF23;
                  p10 = dt*koff*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 2  %% 2: XB3RU Only
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0)
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB3TF23;
                  p10 = dt*koff*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 3    %% 3: Both RURU or XB3RU
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) )
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 4 %% 4: XB2RU Only
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0)
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  p12 = dt*RuOn*XB2TF23;
                  p10 = dt*koff*XB2TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 5 %% 5: Both RURU or XB2RU
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 6 %% 6: Both XB3RU or XB2RU
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB2TF23*XB3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                inextTn=Tn(iRow, iTn+1);
                if ( (max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p12 = dt*RuOn;
                p10 = dt*koff;
            end %% End switch CoopType influence pij

            % Port the normal cases to otherwise of switch above
            %p12 = dt*RuOn;
            %p10 = dt*koff;
            gx = rand;  % seed this puppy ? -- seeded at top of program
            if (gx <= p12)
              TFstate(ii,1)=2;  % The transition to the next state, which makes
              % it availble for binding.
              % This also will transition the next actin on the particular
              % helix
              TFstate(ii+2,1)=2;
              %% Again for the next Tn over (MOD March 11, 2008


            elseif (gx >= 1.0 - p10)
              TFstate(ii, 1)=0;  % The reverse transition.
              % This also will transition the next actin on the particular
              % helix
              % TFstate(ii+2,1)=0;
              %% June 25, 2007
              %% Modify as above with
              %% TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
              %% But use the Coop_2aTo0 function
              %% TFstate(ii+2,1) = Coop_2aTo0( TFstate(Tn(iRow, iTn+1), 1) );

              %% MOD March 11, 2008
              %% Here we need to extend our check on whether the next Tn C
              %% over can control itself, or is a knocked out Tn
              if TnKO(iRow, iTn+1) == 0
                %% Then it is a knockout and should be set at 0, as well as
                %% the TFstate(ii+2,1) in the middle
                TFstate(ii+2,1)=0;
                TFstate(ii+4,1)=0;
              else
                %% Otherwise, it controls itself and we have to querry
                %% using of information from the earlier 56 nm span code:
                %% MOD on June 26, 2007
                %% The Coop_2aTo0 function simply works a switch statement on
                %% the state of the neighboring Tn. Therefore, we can collapse
                %% that into a single line by querying the neighboring Tn
                %% directly.  This only works for To0 (to zero) where the 0 is
                %% as low as we can go, and anything else helps out activation.
                TFstate(ii+2,1) = TFstate(Tn(iRow, iTn+1), 1);
                %% And therefore the following line is really referencing
                %% itself which seems really inefficient.
                TFstate(ii+4,1) = TFstate(Tn(iRow, iTn+1), 1);
              end

            end

          case 2 %ACTNode on==availble for XB binding
            % if the ACTNode is not strongly bound there will be a possible
            % transition from the on to the off state.
            % TFstate(column 2) takes care of the XB binding.
            % unoccupied==0, occupied and weak==1, occupied and
            % strong==2.

            % This will go for the actin node along the Tm as well
            % so we must check both.
            if (~ ( (TFstate(ii, 2)==2)||(TFstate(ii+2, 2)==2)) )
              % compute forward transition likelihoods
              % based on the type of cooperative
              % simululation we are doing. Below looks at
              % TF32 rates--Ru deActivation (p21, prob)
              switch CooperativeType
                case 1 % 1: RURU Only
                  % Again, both neighboring RU
                  inextTn=Tn(iRow, iTn+1);
                  % Again actin [inextTn:2:inextTn+2]
                  % represent all the nodes within the
                  % previous RU.--

                  if (max(TFstate(inextTn:2:inextTn+2, 1)==2)>0)
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 2  %% 2: XB3RU Only
                  inextTn=Tn(iRow, iTn+1);
                  if (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0)
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 3    %% 3: Both RURU or XB3RU
                  %% SEPT MOD bonking on: iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) )
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 4 %% 4: XB2RU Only
                  inextTn=Tn(iRow, iTn+1);
                  if (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0)
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    p21=dt*RuOff*XB2TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 5 %% 5: Both RURU or XB2RU
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 6 %% 6: Both XB3RU or XB2RU
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB2TF32*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 7 %% 7: All three possible;
                  %RURU, XB3RU, XB2RU
                  inextTn=Tn(iRow, iTn+1);
                  if ( (max(TFstate(inextTn:2:inextTn+2, 1)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==2)>0) || (max(TFstate(inextTn:2:inextTn+2, 2)==1)>0) )
                    % XB2TF32=1; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32=1; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32=1; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                otherwise %Not acting on any cooperativity.
                  % No molestarse
                  p21=dt*RuOff; % Transition back to the intermediate state
              end %% End switch CoopType influence pij

              p20=dt*CaOff;  % Cyclic transition back to ground state
              % Move p21 into the otherwise above
              %p21=dt*RuOff; % Transition back to the
              %intermediate state
              gx = rand;
              if (gx <= p20)
                %If the transition is made, then the XB will be knocked off
                % and the position of a bound pair will be set
                % to 0.
                TFstate(ii,:)=0;
                %% MOD March 11, 2008
                %% Here we need to extend our check on whether the next Tn C
                %% over can control itself, or is a knocked out Tn
                if TnKO(iRow, iTn+1) == 0
                  %% Then it is a knockout and should be set at 0, as well as
                  %% the TFstate(ii+2,1) in the middle
                  TFstate(ii+2,:)=0;
                  TFstate(ii+4,:)=0;
                else
                  %% Otherwise, it controls itself and we have to querry
                  %% using of information from the earlier 56 nm span code:
                  TFstate(ii+2,1) = TFstate(Tn(iRow, iTn+1), 1);

                  %% And therefore the following line is really referencing
                  %% itself which seems really inefficient.
                  TFstate(ii+4,1) = TFstate(Tn(iRow, iTn+1), 1);

                  %% Now if either of these sites remains 2 (activated)
                  %% do nothing; otherwise, the XB comes off
                  if TFstate(ii+2,1)<2 && TFstate(ii+4,1)<2
                    TFstate(ii+2,2)=0;  TFstate(ii+2,3)=0;
                    TFstate(ii+4,2)=0;  TFstate(ii+4,3)=0;
                  end
                end

              elseif (gx >= 1.0 - p21)
                TFstate(ii,1)=1;
                %If the transition is made, then the XB will be knocked off
                TFstate(ii, 2)=0;
                TFstate(ii, 3)=0;
                %% MOd March 11, 2008
                if TnKO(iRow, iTn+1) == 0
                  %% Then it is a knockout and should be set at 1, as well as
                  %% the TFstate(ii+2,1) in the middle.  Also the XB comes
                  %% off
                  TFstate(ii+2,1)=1;TFstate(ii+2,2)=0; TFstate(ii+2,3)=0;
                  TFstate(ii+4,1)=1;TFstate(ii+4,2)=0; TFstate(ii+4,3)=0;
                else
                  %% June 25, 2007
                  %% MOD:
                  % TFstate(ii+2, 1)=1; % Tm Span
                  % TFstate(ii+2, 2)=0;% Tm Span on each as well.
                  % TFstate(ii+2, 3)=0;
                  %% So querry our State at the next Tn
                  %% TFstate(Tn(iRow, iTn+1), 1);
                  %% And modify the shared node accordingly
                  TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                  TFstate(ii+4,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                  %% Now if either site remains 2 (activated) do nothing
                  %% otherwise, the site will be 1 and XB comes off
                  %% But we can just querry the other Tn
                  if (TFstate(ii+4,1)==1)
                    TFstate(ii+2,2)=0; TFstate(ii+2,3)=0;
                    TFstate(ii+4,2)=0; TFstate(ii+4,3)=0;
                  end
                end
              end  %% End here

            end % end if on whether there is a strongly bound XB.
        end  % end switch Note: it skips over handel nodes: == -1

      elseif (mod(iTn,(NTn/2))==0) %if iTn=23 or 46
        % Then we are at the free
        % end of the filament and the Tn of choice only
        % controls itself.  mod(i*23, 23)
        % Each Ca state is controlled by itself, but different from the
        % TmSpan of 37 nm before where there is basically no Tm, only
        % Tn on the actin node.  Here however, we will worry about the
        % actin node just to the left of the Tn in question here.
        % Which I believe is going to come up as ii-2 on the
        % TFstate(ii-2,1)=? lines below.

        % Now we need to check the binding:
        % clear p01 p02 p12 p10 p20 p21 %PERFMOD
        switch TFstate(ii, 1)
          case 0 %ACTNode off
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF12 rates--Ca2+binding (p01 in prob.)
            switch CooperativeType
              case 1 % 1: RURU Only
                % At the free end of the actin filament
                % we will only have to search at a single
                % neighboring regulatory unit. The previous
                % Tn site will be found, and then querried
                % to see whether any binding occurs.
                iprevTn=Tn(iRow, iTn-1);
                % Now actin [iprevTn:2:iprevTn+2] should
                % represent all the nodes within the
                % previous RU.  And we can look at either
                % column 1 (thin fil state)  or column 2
                % (XB state) of TFstate to querry
                % cooperativity.--

                if (max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0)
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 2  %% 2: XB3RU Only
                iprevTn=Tn(iRow, iTn-1);
                if (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0)
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 3    %% 3: Both RURU or XB3RU
                iprevTn=Tn(iRow, iTn-1);
                if ((max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 4 %% 4: XB2RU Only
                iprevTn=Tn(iRow, iTn-1);
                if (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0)
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 5 %% 5: Both RURU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                if ((max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 6 %% 6: Both XB3RU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                if ((max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                iprevTn=Tn(iRow, iTn-1);
                if ( (max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p01 = dt*kon*Ca(iCa);
            end %% End switch CoopType influence pij

            %% p02 is reverse transition, un affected
            %% by cooperativity.
            p02 = dt*CaOn*Ca(iCa);
            gx = rand;  % seed this puppy ? -- seeded at top of program
            % If the rand is within the transition prob.
            % Then move to the second state
            % computer forward and reverse transition likelihoods

            if (gx <= p01)
              TFstate(ii,1)=1;  % The transition to the next state
              %% June 25, 2007
              %% MOD: TFstate(ii-2,1)=? query now that TmSpan is
              %% 56 nm.
              %% So querry our State at the prev Tn
              %% iprevTn=Tn(iRow, iTn-1);
              %% And modify the shared node accordingly
              TFstate(ii-2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn-1), 1) );

            elseif (gx >= 1.0 - p02)
              TFstate(ii, 1)=2;  % The reverse transition, which should rarely happen.
              % And the TmSpan
              TFstate(ii-2, 1)=2;
            end

          case 1 %ACTNode is in the intermediate position
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF23 rates--Ru Activation (p12, prob) &
            % TF21 rates--Ca dissociation (p10, prob)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Again, free end, only search 1 prevTn
                iprevTn=Tn(iRow, iTn-1);
                % Again actin [iprevTn:2:iprevTn+2]
                % represent all the nodes within the
                % previous RU.--

                if (max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0)
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*TF3TF23;
                  p10 = dt*koff*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 2  %% 2: XB3RU Only
                iprevTn=Tn(iRow, iTn-1);
                if (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0)
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB3TF23;
                  p10 = dt*koff*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 3    %% 3: Both RURU or XB3RU
                iprevTn=Tn(iRow, iTn-1);
                if ((max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) )
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 4 %% 4: XB2RU Only
                iprevTn=Tn(iRow, iTn-1);
                if (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0)
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  p12 = dt*RuOn*XB2TF23;
                  p10 = dt*koff*XB2TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 5 %% 5: Both RURU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                if ((max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 6 %% 6: Both XB3RU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                if ((max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB2TF23*XB3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                iprevTn=Tn(iRow, iTn-1);
                if ( (max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p12 = dt*RuOn;
                p10 = dt*koff;
            end %% End switch CoopType influence pij

            % Port the normal cases to otherwise of switch above
            %p12 = dt*RuOn;
            %p10 = dt*koff;
            gx = rand;  % seed this puppy ? -- seeded at top of program
            if (gx <= p12)
              TFstate(ii,1)=2;  % The transition to the next state, which makes
              % it availble for binding.  And for the TmSpan too.
              TFstate(ii-2,1)=2;


            elseif (gx >= 1.0 - p10)
              TFstate(ii, 1)=0;  % The reverse transition.
              % This also will transition the next actin
              % on the particular helix: TFstate(ii-2,1)=0;
              %% June 26, 2007
              %% Modify as above with
              %% But use the Coop_2aTo0 function
              %% The Coop_2aTo0 function simply works a switch statement on
              %% the state of the neighboring Tn. Therefore, we can collapse
              %% that into a single line by querying the neighboring Tn
              %% directly.  This only works for To0 (to zero) where the 0 is
              %% as low as we can go, and anything else helps out activation.
              TFstate(ii-2,1) = TFstate(Tn(iRow, iTn-1), 1);
            end

          case 2 %ACTNode on==availble for XB binding
            % if the ACTNode is not strongly bound there will be a possible
            % transition from the on to the off state.
            % TFstate(column 2) takes care of the XB binding.
            % unoccupied==0, occupied and weak==1, occupied and
            % strong==2.
            %% MOD June 28, 2007  Here we will modify the
            %% condition here.  Therefore, if a XB is bound within
            %% the RU (TmSpan) then we cannot have the RU turn off.
            if ~( (TFstate(ii, 2)==2) || (TFstate(ii-2, 2)==2) )
              % compute forward transition likelihoods
              % based on the type of cooperative
              % simululation we are doing. Below looks at
              % TF32 rates--Ru deActivation (p21, prob)
              switch CooperativeType
                case 1 % 1: RURU Only
                  % Again, free end, only search 1 prevTn
                  iprevTn=Tn(iRow, iTn-1);
                  % Again actin [iprevTn:2:iprevTn+2]
                  % represent all the nodes within the
                  % previous RU.--

                  if (max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0)
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 2  %% 2: XB3RU Only
                  iprevTn=Tn(iRow, iTn-1);
                  if (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0)
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 3    %% 3: Both RURU or XB3RU
                  iprevTn=Tn(iRow, iTn-1);
                  if ((max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) )
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 4 %% 4: XB2RU Only
                  iprevTn=Tn(iRow, iTn-1);
                  if (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0)
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    p21=dt*RuOff*XB2TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 5 %% 5: Both RURU or XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  if ((max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 6 %% 6: Both XB3RU or XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  if ((max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB2TF32*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 7 %% 7: All three possible;
                  %RURU, XB3RU, XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  if ( (max(TFstate(iprevTn:2:iprevTn+2, 1)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==2)>0) || (max(TFstate(iprevTn:2:iprevTn+2, 2)==1)>0) )
                    % XB2TF32=1; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32=1; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32=1; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                otherwise %Not acting on any cooperativity.
                  % No molestarse
                  p21=dt*RuOff; % Transition back to the intermediate state
              end %% End switch CoopType influence pij

              p20=dt*CaOff;  % Cyclic transition back to ground state
              % Move p21 into the otherwise above
              %p21=dt*RuOff; % Transition back to the intermediate state
              gx = rand;

              if (gx <= p20)
                %If the transition is made, then the XB will be knocked off
                % and the position of a bound pair will be set
                % to 0.
                TFstate(ii,:)=0;
                %% TFstate(ii-2,:)=0; % And also for Tm span.
                %% MOD on June 26, 2007
                %% The Coop_2aTo0 function simply works a switch statement on
                %% the state of the neighboring Tn. Therefore, we can collapse
                %% that into a single line by querying the neighboring Tn
                %% directly.  This only works for To0 (to zero) where the 0 is
                %% as low as we can go, and anything else helps out activation.
                TFstate(ii-2,1) = TFstate(Tn(iRow, iTn-1), 1);
                %% Now if this site remains 2 (activated) do nothing
                %% otherwise, the XB comes off
                if TFstate(ii-2,1)<2
                  TFstate(ii-2,2)=0;
                  TFstate(ii-2,3)=0;
                end

              elseif (gx >= 1.0 - p21)
                TFstate(ii,1)=1;
                %If the transition is made, then the XB will be knocked off
                TFstate(ii, 2)=0;
                TFstate(ii, 3)=0;
                %% So querry our State at the next Tn
                %% TFstate(Tn(iRow, iTn-1), 1);
                %% And modify the shared node accordingly
                TFstate(ii-2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn-1), 1) );
                %% Now if this site remains 2 (activated) do nothing
                %% otherwise, the site will be 1 and XB comes off
                if TFstate(ii-2,1)==1
                  TFstate(ii-2,2)=0;
                  TFstate(ii-2,3)=0;
                end

              end
            end
        end % end Switch on Cooperativity

      elseif  (mod(iTn,(NTn/2))==22) %if iTn=22 or 45
        % If we are not at the either end, we will have Tm
        % and that will be treated a little differently.
        % Each Tn will controll itself (Super Actin) and
        % 1 actin further down the filament from it.
        % However, when looking towards the free end of the
        % filament for cooperativity, there is not a full
        % RU, so we can only look at the next Tn.
        switch TFstate(ii, 1)
          case 0 %ACTNode off
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF12 rates--Ca2+binding (p01 in prob.)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Not at the free end of the actin filament
                % we have to search at a single
                % neighboring regulatory unit to the left and
                % right of the current RU. The previous
                % Tn site and the next Tn site will be found.
                % Then querried
                % to see whether any binding occurs.
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                % Now actin [iprevTn:2:iprevTn+2, inextTn] should
                % represent all the nodes within the
                % two neighboring RUs.  And we can look at either
                % column 1 (thin fil state)  or column 2
                % (XB state) of TFstate to querry
                % cooperativity.--

                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0)
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 2  %% 2: XB3RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0)
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 3    %% 3: Both RURU or XB3RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 4 %% 4: XB2RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0)
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 5 %% 5: Both RURU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 6 %% 6: Both XB3RU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ( (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p01 = dt*kon*Ca(iCa);
            end %% End switch CoopType influence pij

            %% p02 is reverse transition, un affected
            %% by cooperativity.
            p02 = dt*CaOn*Ca(iCa);
            gx = rand;  % seed this puppy ? -- seeded at top of program
            % If the rand is within the transition prob.
            % Then move to the second state
            % computer forward and reverse transition likelihoods

            if (gx <= p01)
              TFstate(ii,1)=1;  % The transition to the next state
              % This also will transition the next actin on the particular
              % helix
              TFstate(ii-2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn-1), 1) );
              %% MOD: March 11, 2008
              %% For the added span, about what to do at the next TnC
              %% down the pike.  So, if it is one of the TnKO, then it will
              %% be controlled only by the current Tn.  However, if it is
              %% it's individual (not a KO), then the same rule applies as
              %% the TFstate(ii+1, 1) just above.
              if TnKO(iRow, iTn+1) == 0
                %% Then it is a knockout and should be set at 1
                TFstate(ii+4,1)=1;
                TFstate(ii+2,1)=1;
              else
                %% Otherwise, it controls itself and we have to querry
                TFstate(ii+4,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
              end


            elseif (gx >= 1.0 - p02)
              TFstate(ii, 1)=2;  % The reverse transition, which should rarely happen.
              % This also will transition the next actin on the particular
              % helix (left and right)
              TFstate(ii+2,1)=2;
              TFstate(ii-2,1)=2;
              TFstate(ii+4,1)=2;
            end

          case 1 %ACTNode is in the intermediate position
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF23 rates--Ru Activation (p12, prob) &
            % TF21 rates--Ca dissociation (p10, prob)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Again, both neighboring RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                % Again actin [iprevTn:2:iprevTn+2, inextTn]
                % represent all the nodes within the
                % previous RU.--

                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0)
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*TF3TF23;
                  p10 = dt*koff*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 2  %% 2: XB3RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0)
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB3TF23;
                  p10 = dt*koff*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 3    %% 3: Both RURU or XB3RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) )
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 4 %% 4: XB2RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0)
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  p12 = dt*RuOn*XB2TF23;
                  p10 = dt*koff*XB2TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 5 %% 5: Both RURU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 6 %% 6: Both XB3RU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB2TF23*XB3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ( (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p12 = dt*RuOn;
                p10 = dt*koff;
            end %% End switch CoopType influence pij

            gx = rand;  % seed this puppy ? -- seeded at top of program
            if (gx <= p12)
              TFstate(ii,1)=2;  % The transition to the next state, which makes
              % it availble for binding.
              % This also will transition the next actin on the particular
              % helix
              TFstate(ii+2,1)=2;
              TFstate(ii-2,1)=2;
              TFstate(ii+4,1)=2;

            elseif (gx >= 1.0 - p10)
              TFstate(ii, 1)=0;  % The reverse transition.
              %% June 26, 2007
              %% Modify as above with
              %% But use the Coop_2aTo0 function
              %% The Coop_2aTo0 function simply works a switch statement on
              %% the state of the neighboring Tn. Therefore, we can collapse
              %% that into a single line by querying the neighboring Tn
              %% directly.  This only works for To0 (to zero) where the 0 is
              %% as low as we can go, and anything else helps out activation.
              TFstate(ii-2,1) = TFstate(Tn(iRow, iTn-1), 1);
              if TnKO(iRow, iTn+1) == 0 %% Then it is a knockout
                %% and we transition to 1
                TFstate(ii+2,1) = 0;
                TFstate(ii+4,1) = 0;
              else %% we need to check on the state of the neighbor
                TFstate(ii+2,1) = TFstate(Tn(iRow, iTn+1), 1);
                TFstate(ii+4,1) = TFstate(Tn(iRow, iTn+1), 1);
              end
            end

          case 2 %ACTNode on==availble for XB binding
            % if the ACTNode is not strongly bound there will be a possible
            % transition from the on to the off state.
            % TFstate(column 2) takes care of the XB binding.
            % unoccupied==0, occupied and weak==1, occupied and
            % strong==2.

            % This will go for the actin node along the Tm as well
            % so we must check both left and right for the new,
            % larger RU.
            if (~ ( (TFstate(ii-2, 2)==2)||(TFstate(ii, 2)==2)||(TFstate(ii+2, 2)==2)) )
              % compute forward transition likelihoods
              % based on the type of cooperative
              % simululation we are doing. Below looks at
              % TF32 rates--Ru deActivation (p21, prob)
              switch CooperativeType
                case 1 % 1: RURU Only
                  % Again, both neighboring RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  % Again actin [[iprevTn:2:iprevTn+2, inextTn]]
                  % represent all the nodes within the
                  % previous RU.--

                  if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0)
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 2  %% 2: XB3RU Only
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0)
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 3    %% 3: Both RURU or XB3RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) )
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 4 %% 4: XB2RU Only
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0)
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    p21=dt*RuOff*XB2TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 5 %% 5: Both RURU or XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 6 %% 6: Both XB3RU or XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB2TF32*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 7 %% 7: All three possible;
                  %RURU, XB3RU, XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ( (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn], 2)==1)>0) )
                    % XB2TF32=1; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32=1; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32=1; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                otherwise %Not acting on any cooperativity.
                  % No molestarse
                  p21=dt*RuOff; % Transition back to the intermediate state
              end %% End switch CoopType influence pij

              p20=dt*CaOff;  % Cyclic transition back to ground state
              % Move p21 into the otherwise above
              %p21=dt*RuOff; % Transition back to the
              %intermediate state
              gx = rand;
              if (gx <= p20)
                %If the transition is made, then the XB will be knocked off
                % and the position of a bound pair will be set
                % to 0.
                TFstate(ii,:)=0;
                %% MOD March 11, 2008
                %% Here we need to extend our check on whether the next Tn C
                %% over can control itself, or is a knocked out Tn
                if TnKO(iRow, iTn+1) == 0
                  %% Then it is a knockout and should be set at 0, as well as
                  %% the TFstate(ii+2,1) in the middle
                  TFstate(ii+2,:)=0;
                  TFstate(ii+4,:)=0;
                else
                  %% Otherwise, it controls itself and we have to querry
                  %% using of information from the earlier 56 nm span code:
                  TFstate(ii+2,1) = TFstate(Tn(iRow, iTn+1), 1);
                  %% And therefore the following line is really referencing
                  %% itself which seems really inefficient.
                  TFstate(ii+4,1) = TFstate(Tn(iRow, iTn+1), 1);
                  %% Now if either of these sites remains 2 (activated)
                  %% do nothing; otherwise, the XB comes off
                  if TFstate(ii+2,1)<2 && TFstate(ii+4,1)<2
                    TFstate(ii+2,2)=0;  TFstate(ii+2,3)=0;
                    TFstate(ii+4,2)=0;  TFstate(ii+4,3)=0;
                  end
                end
                %% Now the guy to the left
                TFstate(ii-2,1) = TFstate(Tn(iRow, iTn-1), 1);
                %% Now if this site remains 2 (activated) do nothing
                %% otherwise, the XB comes off
                if TFstate(ii-2,1)<2
                  TFstate(ii-2,2)=0;
                  TFstate(ii-2,3)=0;
                end

              elseif (gx >= 1.0 - p21)
                TFstate(ii,1)=1;
                %If the transition is made, then the XB will be knocked off
                TFstate(ii, 2)=0;
                TFstate(ii, 3)=0;
                %% MOd March 11, 2008
                if TnKO(iRow, iTn+1) == 0
                  %% Then it is a knockout and should be set at 1, as well as
                  %% the TFstate(ii+2,1) in the middle.  Also the XB comes
                  %% off
                  TFstate(ii+2,1)=1;TFstate(ii+2,2)=0; TFstate(ii+2,3)=0;
                  TFstate(ii+4,1)=1;TFstate(ii+4,2)=0; TFstate(ii+4,3)=0;
                else
                  %% June 25, 2007
                  %% MOD:
                  % TFstate(ii+2, 1)=1; % Tm Span
                  % TFstate(ii+2, 2)=0;% Tm Span on each as well.
                  % TFstate(ii+2, 3)=0;
                  %% So querry our State at the next Tn
                  %% TFstate(Tn(iRow, iTn+1), 1);
                  %% And modify the shared node accordingly
                  TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                  TFstate(ii+4,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                  %% Now if either site remains 2 (activated) do nothing
                  %% otherwise, the site will be 1 and XB comes off
                  %% But we can just querry the other Tn
                  if (TFstate(ii+4,1)==1)
                    TFstate(ii+2,2)=0; TFstate(ii+2,3)=0;
                    TFstate(ii+4,2)=0; TFstate(ii+4,3)=0;
                  end
                end
                %% Now do the next guy to the left
                TFstate(ii-2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn-1), 1) );
                %% Now if this site remains 2 (activated) do nothing
                %% otherwise, the site will be 1 and XB comes off
                if TFstate(ii-2,1)==1
                  TFstate(ii-2,2)=0;
                  TFstate(ii-2,3)=0;
                end
              end

            end %End if on whether strongly bound
        end % End the Switch on TFState

      else % If we are not at the either end, we will have Tm
        % and that will be treated a little differently.
        % Each Tn will controll itself (Super Actin) and
        % 1 actin further down the filament from it.
        % clear p01 p02 p12 p10 p20 p21 %PERFMOD
        switch TFstate(ii, 1)
          case 0 %ACTNode off
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF12 rates--Ca2+binding (p01 in prob.)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Not at the free end of the actin filament
                % we have to search at a single
                % neighboring regulatory unit to the left and
                % right of the current RU. The previous
                % Tn site and the next Tn site will be found.
                % Then querried
                % to see whether any binding occurs.
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                % Now actin [iprevTn:2:iprevTn+2, inextTn:2:inextTn+2] should
                % represent all the nodes within the
                % two neighboring RUs.  And we can look at either
                % column 1 (thin fil state)  or column 2
                % (XB state) of TFstate to querry
                % cooperativity.--

                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0)
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 2  %% 2: XB3RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0)
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 3    %% 3: Both RURU or XB3RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 4 %% 4: XB2RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0)
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 5 %% 5: Both RURU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                  % XB2TF12; %Feedback on Ca
                  % binding from pre powerstroke state.
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  p01 = dt*kon*Ca(iCa)*TF3TF12*XB2TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 6 %% 6: Both XB3RU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ( (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                  % TF3TF12; %Feedback on Ca
                  % binding from neighboring
                  % actiavted RU.
                  % XB3TF12; %Feedback on Ca
                  % binding from post powerstroke state.
                  % XB2TF12; %Feedback on Ca
                  % binding from pre
                  % powerstroke state.
                  p01 = dt*kon*Ca(iCa)*XB2TF12*XB3TF12*TF3TF12;
                else
                  p01 = dt*kon*Ca(iCa);
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p01 = dt*kon*Ca(iCa);
            end %% End switch CoopType influence pij

            %% p02 is reverse transition, un affected
            %% by cooperativity.
            p02 = dt*CaOn*Ca(iCa);
            gx = rand;  % seed this puppy ? -- seeded at top of program
            % If the rand is within the transition prob.
            % Then move to the second state
            % computer forward and reverse transition likelihoods
            if (gx <= p01)
              TFstate(ii,1)=1;  % The transition to the next state
              % This also will transition the next actin on the particular
              % helix
              TFstate(ii-2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn-1), 1) );
              %% MOD: March 11, 2008
              %% For the added span, about what to do at the next TnC
              %% down the pike.  So, if it is one of the TnKO, then it will
              %% be controlled only by the current Tn.  However, if it is
              %% it's individual (not a KO), then the same rule applies as
              %% the TFstate(ii+1, 1) just above.
              if TnKO(iRow, iTn+1) == 0
                %% Then it is a knockout and should be set at 1
                TFstate(ii+4,1)=1;
                TFstate(ii+2,1)=1;
              else
                %% Otherwise, it controls itself and we have to querry
                TFstate(ii+4,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
              end


            elseif (gx >= 1.0 - p02)
              TFstate(ii, 1)=2;  % The reverse transition, which should rarely happen.
              % This also will transition the next actin on the particular
              % helix (left and right)
              TFstate(ii+2,1)=2;
              TFstate(ii-2,1)=2;
              TFstate(ii+4,1)=2;
            end

          case 1 %ACTNode is in the intermediate position
            % compute forward transition likelihoods
            % based on the type of cooperative
            % simululation we are doing. Below looks at
            % TF23 rates--Ru Activation (p12, prob) &
            % TF21 rates--Ca dissociation (p10, prob)
            switch CooperativeType
              case 1 % 1: RURU Only
                % Again, both neighboring RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                % Again actin [iprevTn:2:iprevTn+2, inextTn:2:inextTn+2]
                % represent all the nodes within the
                % previous RU.--

                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0)
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*TF3TF23;
                  p10 = dt*koff*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 2  %% 2: XB3RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0)
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB3TF23;
                  p10 = dt*koff*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 3    %% 3: Both RURU or XB3RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) )
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 4 %% 4: XB2RU Only
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0)
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  p12 = dt*RuOn*XB2TF23;
                  p10 = dt*koff*XB2TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 5 %% 5: Both RURU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 6 %% 6: Both XB3RU or XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  p12 = dt*RuOn*XB2TF23*XB3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              case 7 %% 7: All three possible;
                %RURU, XB3RU, XB2RU
                iprevTn=Tn(iRow, iTn-1);
                inextTn=Tn(iRow, iTn+1);
                if ( (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                  % XB2TF23; %Feedback on RU activation from pre powerstroke state.
                  % XB3TF23; %Feedback on RU actiavtion from post powerstroke state.
                  % TF3TF23; %Feedback on RU activation from neighboring activated RU.
                  % XB2TF21; %Feedback on Ca release from pre powerstroke state.
                  % XB3TF21; %Feedback on Ca release from post powerstroke state.
                  % TF3TF21; %Feedback on Ca release from neighboring actiavted RU.
                  p12 = dt*RuOn*XB2TF23*XB3TF23*TF3TF23;
                  p10 = dt*koff*XB2TF21*XB3TF21*TF3TF21;
                else
                  p12 = dt*RuOn;
                  p10 = dt*koff;
                end
              otherwise %Not acting on any cooperativity.
                % No molestarse
                p12 = dt*RuOn;
                p10 = dt*koff;
            end %% End switch CoopType influence pij

            % Port the normal cases to otherwise of switch above
            %p12 = dt*RuOn;
            %p10 = dt*koff;
            gx = rand;  % seed this puppy ? -- seeded at top of program
            if (gx <= p12)
              TFstate(ii,1)=2;  % The transition to the next state, which makes
              % it availble for binding.
              % This also will transition the next actin on the particular
              % helix
              TFstate(ii+2,1)=2;
              TFstate(ii-2,1)=2;
              TFstate(ii+4,1)=2;

            elseif (gx >= 1.0 - p10)
              TFstate(ii, 1)=0;  % The reverse transition.
              %% June 26, 2007
              %% Modify as above with
              %% But use the Coop_2aTo0 function
              %% The Coop_2aTo0 function simply works a switch statement on
              %% the state of the neighboring Tn. Therefore, we can collapse
              %% that into a single line by querying the neighboring Tn
              %% directly.  This only works for To0 (to zero) where the 0 is
              %% as low as we can go, and anything else helps out activation.
              TFstate(ii-2,1) = TFstate(Tn(iRow, iTn-1), 1);
              if TnKO(iRow, iTn+1) == 0 %% Then it is a knockout
                %% and we transition to 1
                TFstate(ii+2,1) = 0;
                TFstate(ii+4,1) = 0;
              else %% we need to check on the state of the neighbor
                TFstate(ii+2,1) = TFstate(Tn(iRow, iTn+1), 1);
                TFstate(ii+4,1) = TFstate(Tn(iRow, iTn+1), 1);
              end
            end
            
          case 2 %ACTNode on==availble for XB binding
            % if the ACTNode is not strongly bound there will be a possible
            % transition from the on to the off state.
            % TFstate(column 2) takes care of the XB binding.
            % unoccupied==0, occupied and weak==1, occupied and
            % strong==2.

            % This will go for the actin node along the Tm as well
            % so we must check both left and right for the new,
            % larger RU.
            if (~ ( (TFstate(ii-2, 2)==2)||(TFstate(ii, 2)==2)||(TFstate(ii+2, 2)==2)) )
              % compute forward transition likelihoods
              % based on the type of cooperative
              % simululation we are doing. Below looks at
              % TF32 rates--Ru deActivation (p21, prob)
              switch CooperativeType
                case 1 % 1: RURU Only
                  % Again, both neighboring RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  % Again actin [[iprevTn:2:iprevTn+2, inextTn:2:inextTn+2]]
                  % represent all the nodes within the
                  % previous RU.--

                  if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0)
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 2  %% 2: XB3RU Only
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0)
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 3    %% 3: Both RURU or XB3RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) )
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 4 %% 4: XB2RU Only
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0)
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    p21=dt*RuOff*XB2TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 5 %% 5: Both RURU or XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % TF3TF32; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 6 %% 6: Both XB3RU or XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ((max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                    % XB2TF32; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32; %Feedback on RU deactiavtion from post powerstroke state.
                    p21=dt*RuOff*XB2TF32*XB3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                case 7 %% 7: All three possible;
                  %RURU, XB3RU, XB2RU
                  iprevTn=Tn(iRow, iTn-1);
                  inextTn=Tn(iRow, iTn+1);
                  if ( (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 1)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==2)>0) || (max(TFstate([iprevTn:2:iprevTn+2, inextTn:2:inextTn+2], 2)==1)>0) )
                    % XB2TF32=1; %Feedback on RU deactivation from pre powerstroke state.
                    % XB3TF32=1; %Feedback on RU deactiavtion from post powerstroke state.
                    % TF3TF32=1; %Feedback on RU deactivation from neighboring activated RU.
                    p21=dt*RuOff*XB2TF32*XB3TF32*TF3TF32; % Transition back to the intermediate state
                  else
                    p21=dt*RuOff; % Transition back to the intermediate state
                  end
                otherwise %Not acting on any cooperativity.
                  % No molestarse
                  p21=dt*RuOff; % Transition back to the intermediate state
              end %% End switch CoopType influence pij

              p20=dt*CaOff;  % Cyclic transition back to ground state
              % Move p21 into the otherwise above
              %p21=dt*RuOff; % Transition back to the
              %intermediate state
              gx = rand;
              if (gx <= p20)
                %If the transition is made, then the XB will be knocked off
                % and the position of a bound pair will be set
                % to 0.
                TFstate(ii,:)=0;
                %% MOD March 11, 2008
                %% Here we need to extend our check on whether the next Tn C
                %% over can control itself, or is a knocked out Tn
                if TnKO(iRow, iTn+1) == 0
                  %% Then it is a knockout and should be set at 0, as well as
                  %% the TFstate(ii+2,1) in the middle
                  TFstate(ii+2,:)=0;
                  TFstate(ii+4,:)=0;
                else
                  %% Otherwise, it controls itself and we have to querry
                  %% using of information from the earlier 56 nm span code:
                  TFstate(ii+2,1) = TFstate(Tn(iRow, iTn+1), 1);
                  %% And therefore the following line is really referencing
                  %% itself which seems really inefficient.
                  TFstate(ii+4,1) = TFstate(Tn(iRow, iTn+1), 1);
                  %% Now if either of these sites remains 2 (activated)
                  %% do nothing; otherwise, the XB comes off
                  if TFstate(ii+2,1)<2 && TFstate(ii+4,1)<2
                    TFstate(ii+2,2)=0;  TFstate(ii+2,3)=0;
                    TFstate(ii+4,2)=0;  TFstate(ii+4,3)=0;
                  end
                end
                %% Now the guy to the left
                TFstate(ii-2,1) = TFstate(Tn(iRow, iTn-1), 1);
                %% Now if this site remains 2 (activated) do nothing
                %% otherwise, the XB comes off
                if TFstate(ii-2,1)<2
                  TFstate(ii-2,2)=0;
                  TFstate(ii-2,3)=0;
                end

              elseif (gx >= 1.0 - p21)
                TFstate(ii,1)=1;
                %If the transition is made, then the XB will be knocked off
                TFstate(ii, 2)=0;
                TFstate(ii, 3)=0;
                %% MOd March 11, 2008
                if TnKO(iRow, iTn+1) == 0
                  %% Then it is a knockout and should be set at 1, as well as
                  %% the TFstate(ii+2,1) in the middle.  Also the XB comes
                  %% off
                  TFstate(ii+2,1)=1;TFstate(ii+2,2)=0; TFstate(ii+2,3)=0;
                  TFstate(ii+4,1)=1;TFstate(ii+4,2)=0; TFstate(ii+4,3)=0;
                else
                  %% June 25, 2007
                  %% MOD:
                  % TFstate(ii+2, 1)=1; % Tm Span
                  % TFstate(ii+2, 2)=0;% Tm Span on each as well.
                  % TFstate(ii+2, 3)=0;
                  %% So querry our State at the next Tn
                  %% TFstate(Tn(iRow, iTn+1), 1);
                  %% And modify the shared node accordingly
                  TFstate(ii+2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                  TFstate(ii+4,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn+1), 1) );
                  %% Now if either site remains 2 (activated) do nothing
                  %% otherwise, the site will be 1 and XB comes off
                  %% But we can just querry the other Tn
                  if (TFstate(ii+4,1)==1)
                    TFstate(ii+2,2)=0; TFstate(ii+2,3)=0;
                    TFstate(ii+4,2)=0; TFstate(ii+4,3)=0;
                  end
                end
                %% Now do the next guy to the left
                TFstate(ii-2,1) = Coop_2aTo1( TFstate(Tn(iRow, iTn-1), 1) );
                %% Now if this site remains 2 (activated) do nothing
                %% otherwise, the site will be 1 and XB comes off
                if TFstate(ii-2,1)==1
                  TFstate(ii-2,2)=0;
                  TFstate(ii-2,3)=0;
                end
              end
              
            end
        end  % end switch Note: it skips over handel nodes: == -1
        % Because they are not in the Tn matrix, and we
        % only access the Tn associated with the given
        % actin node.  Then the span is set by how far in
        % front it looks.
      end % End of the if/else loop on free end.  Where we use the MOD opperator to check
      % and see if each Tn is one of the end of the 23 at the
      % free end.
    end % End the If on whether we can go into this guy or not based on TnKO
  end % end loop over the Tn columns.  Where the columns of each row is given in Tn
  % by the Tn on 1 helix, then followed by the Tn on the
  % other helix.  Together this makes up all the Tn on one of
  % the 8 filaments. (Note: there are 46 col. in Tn matrix)
end % end loop over the Tn rows.  Each row of the Tn matrix satisfys one of the
% Thin filaments.  There are eight rows.


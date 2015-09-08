% Main function to run the dynamic model
% Typical usage:
% options = dynamic_xb_prep('mouse');
% ca = your ca transient
% [t,y] = ode15s( @(t,y) make_dy(t,y,interp1(ca_time,ca,mod(t,ca_time(end)),options), 0:1000, options.y0);
%
% dydt = make_dy(t,y,ca,options)
% t: time, not used
% y: current ODE variables
% ca: Ca2+ in micromolar
% options: structure with options and model parameters, use set_options function to obtain this.
% return value dydt: derivatives for the ODE system
function dydt = make_dy(t,y,ca,options)
  n    = options.n;
  nxb  = options.num_xb;
  eps  = 1e-16;
  % STATES
  Tmv     = y(1:options.num_tmxb); 
  TmXb    = zeros(n+1,nxb+1);
  tmxb_ix = sub2ind(size(TmXb),options.state_tmi,options.state_xbj);
  TmXb( tmxb_ix ) = Tmv; % vector -> regular grid of states
  
  TnCB    = y(options.num_tmxb+1);
  TnCC    = y(options.num_tmxb+2);
  TnITnC  = y(options.num_tmxb+3);

  C = options.frac_c' * Tmv;
  B = 1 - C;
  TnI = C - TnITnC;
  
  % Tm kinetics
  TnIFrac =  max(0, TnI ./ max(C,eps));
  
  kbc = options.kbc;
  kcb = options.kcb_no_tni .* TnIFrac;
  xb_off = options.xb_off;
  xb_on  = options.xb_on;
  
  dTm = zeros(n+1,nxb+1);
  tm_bc = kbc(1:n,:)   .* TmXb(1:n  ,:);
  tm_cb = kcb(2:n+1,:) .* TmXb(2:n+1,:);
    
  dTm(2:n+1,:) =              tm_bc - tm_cb; % transition to/from previous tm
  dTm(1:n,:)   = dTm(1:n,:) + tm_cb - tm_bc; % transition to/from next  tm

  assert(~isnan(sum(dTm(tmxb_ix))));
  
  dTm(:,2:nxb+1) = dTm(:,2:nxb+1) + xb_on(:,1:nxb) .* TmXb(:,1:nxb) - xb_off(:,2:nxb+1) .* TmXb(:,2:nxb+1); % xb kinetics
  dTm(:,1:nxb)   = dTm(:,1:nxb)   - xb_on(:,1:nxb) .* TmXb(:,1:nxb) + xb_off(:,2:nxb+1) .* TmXb(:,2:nxb+1);
  
  assert(~isnan(sum(dTm(:))));
  
  % TnC kinetics
  dTnITnC  = options.k_Ion * TnCC - options.k_Ioff * TnITnC;
  
  Jbc = (1/n) * sum(tm_bc(:)) * TnCB/max(B,eps);
  Jcb = (1/n) * sum(tm_cb(:)) * TnCC/max(TnI,eps);
  
  dTnCB  = options.k_Con * ca * (B - TnCB)          - options.k_Coff * TnCB           + Jcb - Jbc;
  dTnCC  = options.k_Con * ca * (C - TnCC - TnITnC) - options.k_Coff * TnCC - dTnITnC + Jbc - Jcb;

  dydt = [dTm(tmxb_ix); dTnCB; dTnCC; dTnITnC];
 
  assert(~isnan(sum(dydt)));




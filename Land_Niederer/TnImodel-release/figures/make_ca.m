% up to tpt: hmt with tau=tpt
% tpt to rt50: solve for tau
% rt50 to rt90: same tau, stretch time
% tpt+2*rt90 to end: linear
function ca = make_ca(dia_ca,delta_ca,tpt,rt50,rt90)
  
  f = @(t,tau) dia_ca + delta_ca * t/tau .* exp(1 - t/tau);

 
  ca1 = f(0:tpt,tpt);

  
  tau1 = fsolve( @(tau) f(tau+rt50,tau) - (dia_ca + delta_ca/2), tpt );

  ca2 = f(tau1+1:tau1+rt50,tau1);

  t5090 = tau1+rt50:tau1+1000;
  ca5090 = f(t5090,tau1);
  
  dt = find(ca5090 <= dia_ca + 0.1*delta_ca,1);
  
  ca3 = interp1(tau1+rt50 + (t5090-tau1-rt50) * (rt90-rt50)/dt,ca5090, tau1+rt50+1:tau1+rt90*2);

  ca4 = interp1( [length([ca1 ca2 ca3]) 1001], [ca3(end) dia_ca], length([ca1 ca2 ca3])+1:1000);
  
  ca = [ca1 ca2 ca3 ca4];
  
  
  ca(end+1) = ca(1);
  
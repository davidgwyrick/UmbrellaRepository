function [] = SimpleForceTime(dt, MFvec, AFvec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

NSTEPS = size(MFvec, 2);
t = (dt*(1:NSTEPS))-dt;

plot(t, MFvec, 'k-'), hold on %%%Plot Thick filament force vs. Time
plot(t, AFvec, 'b-') %%%Plot Thin filament force vs. Time
legend('Thick', 'Thin', 'location', 'southeast')
ylabel('Force (pN)')
xlabel('Time (s)')

end


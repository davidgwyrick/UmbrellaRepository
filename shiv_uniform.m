matlabpool(4)

Dir = '.';

Sets = [531:535; 571:575; 586:590];
kx = [5, 1, 5];
half_sl = [1300, 1150, 1150];
kt = [342, 247, 247];

TnFraction=[0.9,0.8,0.7,0.6,0.5];

for iSet = 1:length(Sets)
    for iTnFraction = 1:length(TnFraction)
        Launch(Dir, Sets(iSet, iTnFraction), 0, kx(iSet), 1, TnFraction(iTnFraction), half_sl(iSet), kt(iSet));
    end
end
function [ThinFMeans] = Plot50(SetNumber)

initfn = ['DataFiles' filesep 'Set' num2str(SetNumber) filesep];

pCaRange=[
    4.0000
    4.5
    5.0
    %    5.25
    5.5
    5.7
    5.8
    5.9
    6.0
    6.1
    6.25
    6.5
    7.0
    7.5
    8.0
    %    9.0000
    ];

pCaLength = length(pCaRange);

ThinFMeans = zeros(1, pCaLength);

for i=1:pCaLength
    data = importdata([initfn 'SSData_pCa_' sprintf('%3.2f', pCaRange(i)) '.txt']);
    ThinFMeans(1, i) = mean(data.data(:, 5));
end

plot(pCaRange, ThinFMeans)
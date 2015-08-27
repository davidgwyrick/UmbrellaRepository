function [mat] = old_randmat(s1,s2,frac)

mat=rand(s1, s2);
[iTest, jTest]=find(mat<=frac);
mat=zeros(s1, s2);
for iii=1:length(iTest)
    %% Now the matrix is zeros, update with the Tn that will actually
    %% exist. As a 1;
    mat(iTest(iii,1), jTest(iii,1) )=1;
end
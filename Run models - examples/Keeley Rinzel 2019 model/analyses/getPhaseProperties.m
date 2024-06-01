function [Mean, Std, rayleighTest, C, Stat] = getPhaseProperties(R)
    C = []; 
    PP = round(R);
    for i = 1:100
        CC = ones(1,PP(i))*3.6*(i-1);
        C = cat(2,C,CC);
    end
    C = C;
    Mean = mean(C);
    Std = std(C);
    C = C/180;
    rayleighTest = [];
    Stat = [];
%     Stat = circ_stats(C);
%     Mean = Stat.mean;
%     Std = Stat.std;
%     rayleighTest = circ_rtest(C);  
end
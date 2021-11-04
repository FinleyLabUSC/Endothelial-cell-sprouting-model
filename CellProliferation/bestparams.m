close all;
clear;

load('r_p.mat');
load('r_s.mat');
load('r_g.mat');
load('r_p_V.mat');
load('r_s_V.mat');
load('r_g_V.mat');

load('results.mat');



 
 GoodFits = [3 21 24 44 49 50 62 66 67 70 82  91 111  126 162 181 183 186 191 195 196]; 

for numfit = 1:size(GoodFits,2)  
    
    besterrors(numfit,1) = Error(GoodFits(numfit),1);
    bestfits(numfit,:) = Fits(GoodFits(numfit),:);
    best_r_p(:,numfit) = r_p(:,GoodFits(numfit));
    best_r_s(:,numfit) = r_s(:,GoodFits(numfit));
    best_r_g(:,numfit) = r_g(:,GoodFits(numfit));
    best_r_p_V(:,numfit) = r_p_V(:,GoodFits(numfit));
    best_r_s_V(:,numfit) = r_s_V(:,GoodFits(numfit));
    best_r_g_V(:,numfit) = r_g_V(:,GoodFits(numfit));
    
end

save('besterrors.mat','besterrors');
save('bestfits.mat','bestfits');
save('best_r_p.mat','best_r_p');
save('best_r_s.mat','best_r_s');
save('best_r_g.mat','best_r_g');
save('best_r_p_V.mat','best_r_p_V');
save('best_r_s_V.mat','best_r_s_V');
save('best_r_g_V.mat','best_r_g_V');
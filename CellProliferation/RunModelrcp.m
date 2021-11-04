clear all
close all
type = 'weighted';

% Best fits #4,
load('initvalue.mat','initvalue');
load('params.mat','params');

load('results.mat');


    
for numfit = 1:4
    numfit
POS = Fits(numfit,:);

initvalue = initvalue(:,1);

%--- Assign the baseline parameter values 
params = params(:,1);

pos =POS;

EC50_pERK = pos(	1	);
EC50_pAkt = pos(	2	);
Vmax_p = pos(	3	);
Vmax_s = (28/50)*Vmax_p;
n_g = pos(	4	);
basal_rg = pos(	5	);



    options = odeset('RelTol',1e-9,'AbsTol',1e-12);
    Ncell = 4e4;
    Avogadro = 6.02214e23;
%     FGFdose_data = West_FGF(:,1); %ng/ml
    Time = [0:60:120*60];
    FGFdose = [0:0.01:0.1 0.2 0.3 0.4 0.5 1:30];
  
    

    r_g(1,numfit) = basal_rg;
    
     % simulate the model
    for i = 2:size(FGFdose,2)      
    FGFconc = FGFdose(i)*1e-9*1e3*1.24e-4/23e3*Avogadro/Ncell;
    % Specify initial values
    initvalue(2,1) = FGFconc; %s2 pape
    initvalue(50,1) = 0; %VEGF
            
    [~, predConc_F] = ode15s(@coreFile_F_V,Time,initvalue,options,params);
    pERK_F(:,i) = sum(predConc_F(:,[7,16,18,19,23,24,31,32]),2); 
    maxERK_F(i,1) = max(pERK_F(:,i));
    
    pAkt_F(:,i) = sum(predConc_F(:,[84,86,87,88,90,91,93]),2); 
    maxAkt_F(i,1) = max(pAkt_F(:,i));
    
    r_p(i,numfit) = Vmax_p/(1+(EC50_pERK/maxERK_F(i,1))^n_g);
    r_s(i,numfit) = Vmax_s/(1+(EC50_pAkt/maxAkt_F(i,1))^n_g);

    r_g(i,numfit) = basal_rg + r_p(i,numfit) + r_s(i,numfit);
    

    end
    

    % Specify initial values
VEGFdose = [0:0.1:1 2:1:10];

    r_g_V(1,numfit) = basal_rg;
    
   for m = 2:size(VEGFdose,2)
    VEGFconc = (7430000/50)*VEGFdose(m);
    initvalue(2,1) = 0; % s34 FGF     
    initvalue(50,1) = VEGFconc;
    % simulate the model
    [~, predConcV] = ode15s(@coreFile_F_V,Time,initvalue,options,params);
    pAkt_V(:,m) = sum(predConcV(:,[84,86,87,88,90,91,93]),2);
    maxpAkt_V(m,1) = max(pAkt_V(:,m));
    
    pERK_V(:,m) = sum(predConcV(:,[7,16,18,19,23,24,31,32]),2);
    maxpERK_V(m,1) = max(pERK_V(:,m));
 
    r_p_V(m,numfit) = Vmax_p/(1+(EC50_pERK/maxpERK_V(m,1))^n_g);
    r_s_V(m,numfit) = Vmax_s/(1+(EC50_pAkt/maxpAkt_V(m,1))^n_g);
    
            
    r_g_V(m,numfit) = basal_rg + r_p_V(m,numfit) + r_s_V(m,numfit);
 

    end
    




save('r_p.mat','r_p');
save('r_s.mat','r_s');
save('r_g.mat','r_g');
save('r_p_V.mat','r_p_V');
save('r_s_V.mat','r_s_V');
save('r_g_V.mat','r_g_V');

end


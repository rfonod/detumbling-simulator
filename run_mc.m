clc, clear

N     = 500; % Number of MC runs
N_orb = 32;  % Number of simulated orbits

data = zeros(25+3*N_orb-1,N); % Data array pre-allocation

tic
parfor_progress(N);
parfor i=1:N
    data_tmp = main_new(i); 
    data(:,i) = data_tmp;
  
    parfor_progress;
end
parfor_progress(0);
toc

cd('MC_Results')
save data_350_32_1_new3_high data
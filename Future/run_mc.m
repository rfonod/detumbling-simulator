clc, clear

N     = 1000; % Number of MC runs
N_orb = 16;  % Number of simulated orbits

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

disp('Dont forget to save the data')

cd('MC_Results')
save data_350_16_1_new data
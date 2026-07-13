%--------------------------------------------------------------------------
%                          MONTE CARLO CAMPAIGN
%   Runs main.m N times with dispersed satellite, actuator, and sensor
%   parameters. Each run returns one packed column of results (see main.m).
%   Requires the Parallel Computing Toolbox for the parfor loop.
%--------------------------------------------------------------------------
clc, clear

N        = 500;    % number of MC runs [-]
N_orb    = 32;     % number of simulated orbits per run [-]
out_name = 'data_350_32_1_new3_high'; % result file, saved into MC_Results/

n_rows = 24 + 3*N_orb;   % height of the packed result column (see main.m)
data   = zeros(n_rows,N);

tic
parfor_progress(N);
parfor i = 1:N
    data(:,i) = main(i,N_orb);

    parfor_progress;
end
parfor_progress(0);
toc

save(fullfile('MC_Results',out_name),'data')
fprintf('Saved %d runs to MC_Results/%s.mat\n',N,out_name);

clear,clc
set(0,'defaultAxesFontSize',15)

%% Simulation Data
ver     = 'new3';          % data version (suffix) /new3 - IAC paper/
N_orb   = 32;              % number of orbits
C_type  = 1;               % control type {0 - modif. B-dot, 1 - weighted B-dot}
alt     = 350e3;           % altitude (LEO) [m]
%alt    = 200e3;          % altitude (LEO) [m]
re      = 6378.137e3;      % Earth (equatorial) radius [m]
a       = re + alt;        % semi-major axis (assuming circular orbit) [m]
mu      = 3.9860044181e14; % Earth's standard gravitational parameter [m^3/s^2]
n       = sqrt(mu/a^3);    % mean motion [rad/s]
T_orb   = 2*pi/n;          % orbital period [s]

cmp_act = 1;              % data comparision active?

%% Load 1st data set and compute some statistics
orbit_str = sprintf('data_%d_%d_%d_%s.mat',alt/1e3,N_orb,C_type,ver);

% Load simulation results and extract data
load(orbit_str)

N_mc = length(data);
data = data(:,~isnan(data(1,:)'));

t_det_w  = data(1,:)/T_orb; % true detumbled time [orbits]
t_det_p  = data(2,:)/T_orb; % algorithm estimated detumbled time [orbits]
mass     = data(3,:);       % mass [kg]
Id       = data(4:6,:);     % mass moment of inertia [kg.m^2]
m_max    = data(7:9,:);     % maximum magnetic dipole moment [A.m^2]
m_res    = data(10:12,:);   % maximum residual mag. dipole moment [A.m^2]
p_tmb    = data(13:15,:);   % initial tumbling parameter [T/s]
w0       = rad2deg(data(16:18,1)); % initial angular rate [deg/s]
t_on_w   = data(19:21,:);          % sum of t_on (per axis) until t_det_w [s]
t_on_p   = data(22:24,:);          % sum of t_on (per axis) until t_det_p [s]

if size(data,1) == 56
    data_ver = 0;
    t_on_orb = data(25:25+N_orb-1,:);  % sum of t_on (all axis) per orbit [s]
else
    data_ver = 1;
    t_on_orb = nan(N_orb,3,size(data,2));
    for i = 1:N_orb
        t_on_orb(i,:,:)  = data(25+(i-1)*3:25+i*3-1,:);
    end
end

% Compute some statistics
% in orbits
t_det_w_mean = mean(t_det_w);
t_det_w_med  = median(t_det_w);
t_det_w_std  = std(t_det_w);

% in hours
t_err_mean = mean(T_orb*(t_det_p-t_det_w)/3600);
t_err_med  = median(T_orb*(t_det_p-t_det_w)/3600);
t_err_std  = std(T_orb*(t_det_p-t_det_w)/3600);


% in seconds
t_on_w_mean = mean(t_on_w,2);
t_on_w_med  = median(t_on_w,2);
t_on_w_std  = std(t_on_w,0,2);
t_on_w_mean_sum  = sum(t_on_w_mean);
t_on_w_std_sum = std(sum(t_on_w,1));


% in seconds
t_on_p_mean = mean(t_on_p,2);
t_on_p_med  = median(t_on_p,2);
t_on_p_std  = std(t_on_p,0,22);

% in seconds
if size(data,1) == 56
    t_on_orb_sum_mean = nanmean(t_on_orb,2);
    t_on_orb_sum_std  = nanstd(t_on_orb,0,2);
else
    t_on_orb_sum_mean = nanmean(t_on_orb,3);
    t_on_orb_sum_std  = nanstd(t_on_orb,0,3);
end

%% Labels
x_text_orb = sprintf('Number of orbits (1 orbit = %g hours; 1 day = %g orbits)',...
    round(T_orb/3600,2),round(24/(T_orb/3600),2));
x_text_hrs = sprintf('Number of hours (1 orbit = %g hours)',...
    round(T_orb/3600,2));
y_text = sprintf('Frequency (N_{MC} = %g)',N_mc);

%% Load 2nd data set
if cmp_act
    orbit_str = sprintf('data_%d_%d_%d_%s.mat',alt/1e3,N_orb,~C_type,ver);
    
    % Load simulation results and extract data
    load(orbit_str)
    
    data = data(:,~isnan(data(1,:)'));
    
    t_det_w2  = data(1,:)/T_orb; % true detumbled time [orbits]
    t_det_p2  = data(2,:)/T_orb; % algorithm estimated detumbled time [orbits]
    mass2     = data(3,:);       % mass [kg]
    Id2       = data(4:6,:);     % mass moment of inertia [kg.m^2]
    m_max2    = data(7:9,:);     % maximum magnetic dipole moment [A.m^2]
    m_res2    = data(10:12,:);   % maximum residual mag. dipole moment [A.m^2]
    p_tmb2    = data(13:15,:);   % initial tumbling parameter [T/s]
    w02       = rad2deg(data(16:18,1)); % initial angular rate [deg/s]
    t_on_w2   = data(19:21,:);          % sum of t_on (per axis) until t_det_w [s]
    t_on_p2   = data(22:24,:);          % sum of t_on (per axis) until t_det_p [s]
    
    if size(data,1) == 56
        data_ver2 = 0;
        t_on_orb2 = data(25:25+N_orb-1,:);  % sum of t_on (all axis) per orbit [s]
    else
        data_ver2 = 1;
        t_on_orb2 = nan(N_orb,3,size(data,2));
        for i = 1:N_orb
            t_on_orb2(i,:,:)  = data(25+(i-1)*3:25+i*3-1,:);
        end
    end
    
    % Compute some statistics
    % in orbits
    t_det_w_mean2 = mean(t_det_w2);
    t_det_w_med2  = median(t_det_w2);
    t_det_w_std2  = std(t_det_w2);
    
    % in hours
    t_err_mean2 = mean(T_orb*(t_det_p2-t_det_w2)/3600);
    t_err_med2  = median(T_orb*(t_det_p2-t_det_w2)/3600);
    t_err_std2  = std(T_orb*(t_det_p2-t_det_w2)/3600);
    
    % in seconds
    t_on_w_mean2 = mean(t_on_w2,2);
    t_on_w_med2  = median(t_on_w2,2);
    t_on_w_std2  = std(t_on_w2,0,2);
    t_on_w_mean_sum2  = sum(t_on_w_mean2);
    t_on_w_std_sum2 = std(sum(t_on_w2,1));
    
    % in seconds
    t_on_p_mean2 = mean(t_on_p2,2);
    t_on_p_med2  = median(t_on_p2,2);
    t_on_p_std2  = std(t_on_p2,0,22);
    
    % in seconds
    if ~data_ver2
        t_on_orb_sum_mean2 = nanmean(t_on_orb2,2);
        t_on_orb_sum_std2  = nanstd(t_on_orb2,0,2);
    else
        t_on_orb_sum_mean2 = nanmean(t_on_orb2,3);
        t_on_orb_sum_std2  = nanstd(t_on_orb2,0,3);
    end
end

%% Tabelized data
if cmp_act
    disp('Detumbling time in [hours]')
    disp(['   Mean (1)   ','Std (1)  ','Mean (2)   ','Std (2)   ', '% difference'])
    disp([T_orb*[t_det_w_mean,t_det_w_std,t_det_w_mean2,t_det_w_std2]/3600,...
        [100*(t_det_w_mean-t_det_w_mean2)/t_det_w_mean2]])
    
    disp('Power consumption (per axis) in [hours]')
    disp(['   Mean (1)   ','Std (1)  ','Mean (2)   ','Std (2)   ', '% difference'])
    disp([[t_on_w_mean,t_on_w_std,t_on_w_mean2,t_on_w_std2]/3600,...
        [100*(t_on_w_mean-t_on_w_mean2)./t_on_w_mean2]])
    
    disp(['   Mean (1)   ','Std (1)  ','Mean (2)   ','Std (2)   ', '% difference'])
    disp([[t_on_w_mean_sum,t_on_w_std_sum,t_on_w_mean_sum2,t_on_w_std_sum2]/3600,...
        [100*(t_on_w_mean_sum-t_on_w_mean_sum2)/t_on_w_mean_sum2]])
else
    disp('Detumbling time in [orbits]')
    disp(['   Mean (1)   ','Std (1)  '])
    disp([t_det_w_mean,t_det_w_std])
    
    disp('Power consumption (per axis) in [hours]')
    disp(['   Mean (1)   ','Std (1)  '])
    disp([t_on_w_mean,t_on_w_std]/3600)
end

%% Plot 0 - Parameter - Histograms
figure

subplot(2,3,1)
histfit(mass), grid minor
xlabel('Mass [kg]')
ylabel('Frequency')
subplot(2,3,2)
histfit(m_max(1,:)), grid minor
xlabel('Max. dipole [A.m^2]')
%ylabel('Frequency')
subplot(2,3,3)
histfit(vecnorm(m_max)), grid minor
xlabel('Res. mag. dipole (magnitude) [A.m^2]')
%ylabel('Frequency')


subplot(2,3,4)
histfit(Id(1,:)), grid minor
xlabel('$I_x$ [kg.m$^2$]','interpreter','latex')
ylabel('Frequency')
subplot(2,3,5)
histfit(Id(1,:)), grid minor
xlabel('$I_y$ [kg.m$^2$]','interpreter','latex')
%ylabel('Frequency')
subplot(2,3,6)
histfit(Id(1,:)), grid minor
xlabel('$I_z$ [kg.m$^2$]','interpreter','latex')
%ylabel('Frequency')

%% Plot 1 - True detumbling time (Orbits)
figure
histfit(t_det_w,25), grid minor, hold on
%title('Detumbling time from 180 deg/s (in all 3 axes)')
xlabel(x_text_orb,'fontsize',15)
%ylabel(y_text)
ylabel('Frequency','fontsize',15)
legend('Detumbling time','Gaussian fit')

text(.5,.8,'Median: ','Units','normalized','FontSize',14,'fontweight','bold')
text(.5,.75,'Mean: ','Units','normalized','FontSize',14,'fontweight','bold')
text(.5,.7,'Std: ','Units','normalized','FontSize',14,'fontweight','bold')


text(.64,.8,sprintf('%g orbits / %g hours',round(t_det_w_med,1),...
    round(t_det_w_med*T_orb/3600,1)),'Units','normalized','FontSize',14)
text(.64,.75,sprintf('%g orbits / %g hours',round(t_det_w_mean,1),...
    round(t_det_w_mean*T_orb/3600,1)),'Units','normalized','FontSize',14)
text(.66,.7,sprintf('%g orbits / %g hours',round(t_det_w_std,1),...
    round(t_det_w_std*T_orb/3600,1)),'Units','normalized','FontSize',14)
% if savePlots
    fig_text=sprintf(['hist.eps']);
    print('-depsc2',fig_text)
% end

return
%% Plot 2 - Alg. estimated detumbling time
figure
histfit(T_orb*(t_det_p-t_det_w)/3600), grid minor
title('Algorithm detumbling confirmation time delay')
xlabel('Hours')
ylabel(y_text)
legend('Alg. confirmation time delay','Gaussian fit')

text(.03,.95,'Mean: ','Units','normalized','FontSize',14,'fontweight','bold')
text(.03,.9,'Std: ','Units','normalized','FontSize',14,'fontweight','bold')
text(.03,.85,'Median: ','Units','normalized','FontSize',14,'fontweight','bold')

text(.17,.95,sprintf('%g hours',round(t_err_mean,1)),'Units','normalized','FontSize',14)
text(.17,.90,sprintf('%g hours',round(t_err_std,1)),'Units','normalized','FontSize',14)
text(.17,.85,sprintf('%g hours',round(t_err_med,1)),'Units','normalized','FontSize',14)

annotation(gcf,'textbox',...
    [0.151785714285714 0.682333333333338 0.333482142857143 0.0595238095238095],...
    'String',{'Confirmation window: 1 hour!'},...
    'LineWidth',1,...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Helvetica Neue',...
    'BackgroundColor',[0.301960784313725 0.749019607843137 0.929411764705882]);

%% Plot 3 - Correlations - Inertia - Separate terms

% I_x
figure

subplot(2,3,1)
%corrplot([Id(1,:)',t_det_w']), grid minor
scatter(Id(1,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(Id,2),1),Id(1,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
%title('Detumbling time correlation with I_x')
ylabel('$t_{det}$ [Orbits]','interpreter','latex')
xlabel('$I_x$ [kg.m$^2$]','interpreter','latex')
legend('Samples','Correlation')

R = corr(Id(1,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

% I_y
%figure
subplot(2,3,2)
%corrplot([Id(2,:)',t_det_w']), grid minor
scatter(Id(2,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(Id,2),1),Id(2,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
%title('Detumbling time correlation with I_y')
%ylabel('$t_{det}$ [Orbits]','interpreter','latex')
xlabel('$I_y$ [kg.m$^2$]','interpreter','latex')
%legend('Samples','Correlation')

R = corr(Id(2,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

% I_z
%figure
subplot(2,3,3)
%corrplot([Id(3,:)',t_det_w']), grid minor
scatter(Id(3,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(Id,2),1),Id(3,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
%title('Detumbling time correlation with I_z')
%ylabel('$t_{det}$ [Orbits]','interpreter','latex')
xlabel('$I_z$ [kg.m$^2$]','interpreter','latex')
%legend('Samples','Correlation')

R = corr(Id(3,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')


for i=1:3
    subplot(2,3,3+i)
    
    scatter(m_max(i,:)',t_det_w','filled'), grid minor, hold on
    kk = [ones(size(Id,2),1),m_max(i,:)']\t_det_w';
    x_ = xlim;
    xx = linspace(x_(1),x_(2),1e3);
    yy = kk(1) + kk(2)*xx;
    plot(xx,yy,'r','linewidth',1.5)
    %title('Detumbling time correlation with m_{max}')
    xlabel('$\bar m_x$ [A.m$^2$]','interpreter','latex')
    %legend('Samples','Correlation')
    
    R = corr(m_max(i,:)',t_det_w');
    text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')
end
subplot(2,3,4)
ylabel('$t_{det}$ [Orbits]','interpreter','latex')
subplot(2,3,5)
ylabel('$t_{det}$ [Orbits]','interpreter','latex')
xlabel('$\bar m_y$ [A.m$^2$]','interpreter','latex')
subplot(2,3,6)
xlabel('$\bar m_z$ [A.m$^2$]','interpreter','latex')
%% Plot 3 - Correlations - Inertia - All

% I_all
figure
%corrplot([vecnorm(Id)',t_det_w']), grid minor
scatter(vecnorm(Id)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(Id,2),1),vecnorm(Id)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with inertia magnitude')
ylabel('Detumbling time [Orbits]')
xlabel('Euclidian norm of the inertia [kg.m^2]')
legend('Samples','Correlation')

R = corr(vecnorm(Id)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

%% Plot 3b - Correlations - Mass
figure
%corrplot([mass',t_det_w']), grid minor
scatter(mass',t_det_w','filled'), grid minor, hold on
kk = [ones(size(mass,2),1),mass(1,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with mass')
ylabel('Detumbling time [Orbits]')
xlabel('Mass [kg]')
legend('Samples','Correlation')

R = corr(mass(1,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

%% Plot 4 - Correlations - Max dipole
figure
%corrplot([m_max(1,:)',t_det_w']), grid minor
scatter(m_max(1,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(Id,2),1),m_max(1,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with m_{max}')
ylabel('Detumbling time [Orbits]')
xlabel('Maximal magnetic dipole moment [A.m^2]')
legend('Samples','Correlation')

R = corr(m_max(1,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

%% Plot 5 - Correlation - Residual dipole magnitude - separate
% m_res_x
figure

subplot(1,3,1)
%corrplot([m_res(1,:)',t_det_w']), grid minor
scatter(m_res(1,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(m_res,2),1),m_res(1,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with m_{res}^x')
ylabel('Detumbling time [Orbits]')
xlabel('Residual mag. dipole along X-axis [A.m^2]')
legend('Samples','Correlation')

R = corr(m_res(1,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

% m_res_y
%figure
subplot(1,3,2)
%corrplot([m_res(2,:)',t_det_w']), grid minor
scatter(m_res(2,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(m_res,2),1),m_res(2,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with m_{res}^y')
ylabel('Detumbling time [Orbits]')
xlabel('Residual mag. dipole along y-axis [A.m^2]')
legend('Samples','Correlation')

R = corr(m_res(2,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

% m_res_z
%figure
subplot(1,3,3)
%corrplot([m_res(3,:)',t_det_w']), grid minor
scatter(m_res(3,:)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(m_res,2),1),m_res(3,:)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with m_{res}^z')
ylabel('Detumbling time [Orbits]')
xlabel('Residual mag. dipole along z-axis [A.m^2]')
legend('Samples','Correlation')

R = corr(m_res(3,:)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')


%% Plot 6 - Correlation - Residual dipole magnitude - all
figure
%corrplot([vecnorm(m_res)',t_det_w']), grid minor
scatter(vecnorm(m_res)',t_det_w','filled'), grid minor, hold on
kk = [ones(size(m_res,2),1),vecnorm(m_res)']\t_det_w';
x_ = xlim;
xx = linspace(x_(1),x_(2),1e3);
yy = kk(1) + kk(2)*xx;
plot(xx,yy,'r','linewidth',1.5)
title('Detumbling time correlation with ||m_{res}||')
ylabel('Detumbling time [Orbits]')
xlabel('Euclidian norm of the residual mag. dipole moment [A.m^2]')
legend('Samples','Correlation')

R = corr(vecnorm(m_res)',t_det_w');
text(.09,.9,sprintf('R = %g',round(R,2)),'Units','normalized','FontSize',14,'fontweight','bold')

%% Plot 7 - Average consumption
x_text = sprintf('Number of orbits (1 orbit = %g hours; 1 day = %g orbits)',...
    round(T_orb/3600,2),round(24/(T_orb/3600),2));

if ~data_ver
    figure
    
    lfz = 16;
    plot(1:N_orb,t_on_orb_sum_mean','bo-','markeredgecolor','b',...
        'markerfacecolor','r','markersize',5),grid,hold on
    plot(1:N_orb,t_on_orb_sum_mean'+[1;-1]*t_on_orb_sum_std','ks-.',...
        'markeredgecolor','k','markersize',3)
    
    h=legend('mean','$1\sigma$ envelope');
    set(h,'FontSize',lfz,'interpreter','latex','location','SouthWest')
    ylabel('Relative intercept angle error [deg]')
else
    figure
    
    lfz = 15;
    for i=1:3
        subplot(3,1,i)
        plot(1:N_orb,t_on_orb_sum_mean(:,i)'/3600,'ko-','markeredgecolor','b',...
            'markerfacecolor','r','markersize',5),grid minor,hold on
        plot(1:N_orb,t_on_orb_sum_mean2(:,i)'/3600,'ks-','markeredgecolor','b',...
            'markerfacecolor','k','markersize',5)
        %plot(1:N_orb,t_on_orb_sum_mean(:,i)'/3600+[1;-1]*t_on_orb_sum_std(:,i)'/3600,'ks-.',...
        %    'markeredgecolor','k','markersize',3)
        if i==1
        h=legend('Weighted B-dot','B-dot of G.Avanzini & F.Giulietti');
        set(h,'FontSize',lfz) %,'interpreter','latex','location','SouthWest')
        end
    end
    subplot(3,1,1)
    ylabel('$\mu(T_{orb;x}^{on})$ [hours]','interpreter','latex')
    subplot(3,1,2)
    ylabel('$\mu(T_{orb;y}^{on})$ [hours]','interpreter','latex')
    ylim([0.6,1])
    subplot(3,1,3)
    ylabel('$\mu(T_{orb;z}^{on})$ [hours]','interpreter','latex')
    ylim([0.7,1])
    xlabel(x_text)
    
% if savePlots
%     fig_text=sprintf(['power_orb.eps']);
%     print('-depsc2',fig_text)
%     %Process_AtendHeader(fig_text,fig_text);
% end
    
end


add_text = '';
set(0,'defaultAxesFontSize',13)
x_text = sprintf('Number of orbits (1 orbit = %g hours; 1 day = %g orbits)',...
    round(T_orb/3600,2),round(24/(T_orb/3600),2));

%% Effect of magnetometer errors on the measurements
figure
nn = 1; % number of orbits considered for this plot 
t_mk = floor(nn*t_len/N_orb);
plot3(B_B(1,1:t_mk),B_B(2,1:t_mk),B_B(3,1:t_mk),'r.','markersize',4),hold,grid
plot3(B_B_meas(1,1:t_mk),B_B_meas(2,1:t_mk),B_B_meas(3,1:t_mk),'k.','markersize',4)
xlabel('$b_{E;x}^B$ [$\mu$T]','interpreter','latex')
ylabel('$b_{E;y}^B$ [$\mu$T]','interpreter','latex')
zlabel('$b_{E;z}^B$ [$\mu$T]','interpreter','latex')
h=legend('$\mathbf{b}_E^B$','$\mathbf{z}_{raw}^{(1)}$');
set(h,'interpreter','latex')

if savePlots
    fig_text=sprintf(['b_circ',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

%% Delivered magnetic dipole moment
% figure
% subplot(3,1,1)
% plot(t/T_orb,M_del(1,:)), grid on
% xlabel(x_text)
% ylabel('$m^{del}_x$','interpreter','latex')
% subplot(3,1,2)
% plot(t/T_orb,M_del(2,:)), grid on
% xlabel(x_text)
% ylabel('$m^{del}_y$','interpreter','latex')
% subplot(3,1,3)
% plot(t/T_orb,M_del(3,:)), grid on
% xlabel(x_text)
% ylabel('$m^{del}_z$','interpreter','latex')
% 
% if savePlots
%     fig_text=sprintf(['mdel',add_text,'.eps']);
%     print('-depsc2',fig_text)
%     %Process_AtendHeader(fig_text,fig_text);
% end

%% Commanded magnetorquer activation time in %
figure
subplot(7,1,1:2)
plot(t/T_orb,100*S_on(1,:).*T_on(1,:)/delta/T_c), grid on
xlabel(x_text)
ylabel('$t^{on}_x$ [\% cycle]','interpreter','latex')
legend(sprintf('INT = %d',sum(abs(T_on(1,:)))));
subplot(7,1,3:4)
plot(t/T_orb,100*S_on(2,:).*T_on(2,:)/delta/T_c), grid on
xlabel(x_text)
ylabel('$t^{on}_y$ [\% cycle]','interpreter','latex')
legend(sprintf('INT = %d',sum(abs(T_on(2,:)))));
subplot(7,1,5:6)
plot(t/T_orb,100*S_on(3,:).*T_on(3,:)/delta/T_c), grid on
xlabel(x_text)
ylabel('$t^{on}_z$ [\% cycle]','interpreter','latex')
legend(sprintf('INT = %d',sum(abs(T_on(3,:)))));
subplot(7,1,7)
plot(1:N_orb,t_on_orb,'*-'), grid on
xlabel(x_text)
ylabel('$\sum t^{on}$ ','interpreter','latex')
legend(sprintf('Sum of t_on per orbit basis',sum(abs(T_on(3,:)))));
if savePlots
    fig_text=sprintf(['ton',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

%% Angular Speed Only
figure
plot(t/T_orb,W*r2d,'linewidth',1), grid on
ylabel('Angular rates [deg/s]')
h=legend('$\omega_x$','$\omega_y$','$\omega_z$');
set(h,'fontsize',15)
set(h,'interpreter','latex')
xlabel(x_text)
text(.22,.1,sprintf('Detumbled after: %g orbits = %g hrs',round(t_det_w/T_orb,1),round(t_det_w/3600,1)),'Units','normalized','FontSize',15,'fontweight','bold')
text(.22,.05,sprintf('Algorithm confirmed after: %g orbits = %g hrs',round(t_det_p/T_orb,1),round(t_det_p/3600,1)),'Units','normalized','FontSize',15,'fontweight','bold')

if savePlots
    fig_text=sprintf(['w',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

%% Magnetic field in the body frame (true) and B_dot (estimated)
figure
subplot(3,1,1)
plot(t/T_orb,B_B), grid on
xlabel(x_text)
ylabel('b_B (true) [\mu T]')
h=legend('$b_x$','$b_y$','$b_z$');
set(h,'interpreter','latex')
subplot(3,1,2)
plot(t/T_orb,B_B_meas), grid on
xlabel(x_text)
ylabel('b_B (meas) [\mu T]')
h=legend('$b_x$','$b_y$','$b_z$');
set(h,'interpreter','latex')
subplot(3,1,3)
plot(t/T_orb,B_dot), grid on
ylabel('B\_dot (meas) [\mu T/s]')
h=legend('$\dot b_x$','$\dot b_y$','$\dot b_z$');
set(h,'interpreter','latex')
xlabel(x_text)

if savePlots
    fig_text=sprintf(['bbdot',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

%% Tumble parameter 
figure
semilogy(t/T_orb,P_tmb), grid on,hold on
semilogy([t(1)/T_orb t(end)/T_orb],[p_bar_u p_bar_u],'g--','linewidth',1)
semilogy([t(1)/T_orb t(end)/T_orb],[p_bar_l p_bar_l],'k:','linewidth',2)
ylabel('P_{tumble} [1/s]')
h=legend('$p^{tumb}_x$','$p^{tumb}_y$','$p^{tumb}_z$','$\bar p_{upp}$','$\bar p_{low}$');
set(h,'interpreter','latex')
xlabel(x_text)

if savePlots
    fig_text=sprintf(['ptumb',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

%% Quaternion, Angular Speed, and Angular Momentum Magnitude
figure
subplot(10,1,[1:3])
plot(t/T_orb,Q), grid on
xlabel(x_text)
ylabel('Quaternions [-]')
h=legend('$q_1$','$q_2$','$q_3$','$q_4$');
set(h,'interpreter','latex')
subplot(10,1,[5:7])
plot(t/T_orb,W*r2d), grid on
ylabel('Angular rates [deg/s]')
h=legend('$\omega_x$','$\omega_y$','$\omega_z$');
set(h,'interpreter','latex')
xlabel(x_text)

subplot(10,1,[9:10])
plot(t/T_orb,vecnorm(I*W),'k','linewidth',1.5), grid on
ylabel('||H|| [N.m.s]')
xlabel(x_text)

if savePlots
    fig_text=sprintf(['qw',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

%% Ground Track
% Ground station coordinates
lat_gs = 51.9986;
lon_gs = 4.3736;

% Cartesian to spherical transformation
[lon,lati,~] = cart2sph(R_E(1,:),R_E(2,:),R_E(3,:));

figure
load coast % requires Mapping Toolbox
plot(long,lat,'k'), hold on

pa = plot(r2d*lon(1),r2d*lati(1),'b','LineWidth',2);
pb = plot(r2d*lon(1),r2d*lati(1),'k','LineWidth',3);
pc = plot(r2d*lon(1),r2d*lati(1),'g','LineWidth',3);
pd = plot(lon_gs,lat_gs,'yh','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','y');

plot(r2d*lon,r2d*lati,'b.','MarkerSize',3,'MarkerEdgeColor','b');

pe = plot(r2d*lon(1),r2d*lati(1),'ro','MarkerSize',9,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');

axis([-180 180 -90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
%legend([pa pb pc pd pe],'Satellite ground track','Eclipse period',...
%    'Ground station pass','Delft ground station','Start of simulation')
legend([pa pd pe],'Satellite ground track','Delft ground station',...
    'Start of simulation')

if savePlots
    fig_text=sprintf(['orbits',add_text,'.eps']);
    print('-depsc2',fig_text)
    %Process_AtendHeader(fig_text,fig_text);
end

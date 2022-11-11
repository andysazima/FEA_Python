clc
clear
close all
format shortg
tic



%% Inputs

% Physical Inputs
Ri = 10; % Inner radius of the ball
Ro = 20; % Outer radius of the ball
L_tot = Ro - Ri;

% Discretization
N_elements = 100; % Number of elemetns

% Integration parameters
Ngp_M = 3; % Gauss points from 1 to 6
Ngp_K = 3; % Gauss points from 1 to 6

% Analysis Parameters
beta  = 0.0;
gamma = 0.5;
tcr_ratio = 0.99; % ratio of t_cr to use for time step
t_total = 5;      % total time to be calculated

% Force BCs
Po = 1.0;
Fo = 4*pi*Ri^2*Po;
Force_BCs = [1 Fo]; % First DOF with force Fo. No others with a force BC.

% Elastic Constants
E = 100;
v = 0.3;
rho = 0.01;
lam = E*v/(1+v)/(1-2*v);
mu = E/2/(1+v);




%% Calcs

% Finish Discretization
N_nodes = N_elements + 1;

% Initial Conditions set to 0
do = zeros(N_nodes,1);
vo = zeros(N_nodes,1);

% Select lumped mass for Central Difference automatically
if beta == 0
    MassType = 'lumped'; % 'consistent' or 'lumped'
else
    MassType = 'consistent'; % 'consistent' or 'lumped'
end

% Initialize the vector containing the nodal coordinates
x = Ri;
Nodal_Positions = [Ri];
for iN = 2:N_nodes
    x = x + L_tot/N_elements;
    Nodal_Positions(iN,1) = x;
end
elements = zeros(N_elements,2);
for iE=1:N_elements
    st_pt=Nodal_Positions(iE,1);
    end_pt=Nodal_Positions(iE+1,1);
    elements(iE,:)=[st_pt end_pt];
end

% Create constants matrix
D = [lam+2*mu   lam        lam;
     lam        lam+2*mu   lam;
     lam        lam        lam+2*mu];
L_ele = Preprocessing(1,L_tot,N_elements);

% Calc K, C, M, and F
[K, K_loc] = StiffnessMatrix(L_ele,Nodal_Positions,D,Ngp_K);
C = zeros(size(K));
[M, M_loc] = MassMatrix(L_ele,Nodal_Positions,rho,Ngp_M,MassType);
F = ForceVector(K,Force_BCs);

% Stability analysis
EigVals = zeros(size(K_loc));
EigMins = zeros(size(K_loc,3),1);
EigMaxs = zeros(size(K_loc,3),1);
for n = 1:size(K_loc,3)
    [~, EigVals(:,:,n)] = eig(K_loc(:,:,n), M_loc(:,:,n));
    EigMins(n) = min(EigVals(EigVals > 0),[],"all");
    EigMaxs(n) = max(EigVals(EigVals > 0),[],"all");
end
EigMin = min(EigMins,[],"all");
EigMax = max(EigMaxs,[],"all");
Omega_cr = 1 / (sqrt(gamma/2 - beta));
t_cr = Omega_cr / sqrt(EigMax);

% Newmark Beta Method
if beta >= 0.25
    t_cr = 0.05 / N_elements;
    t_cr_text = ['Uncond Stab, use t_cr = ' num2str(t_cr)];
else
    t_cr_text = num2str(t_cr);
end
dt = tcr_ratio*t_cr;
time = 0:dt:t_total;
[d,v,a,algo_name] = NewmarkBetaBall(M,C,K,F,do,vo,beta,gamma,dt,t_total,Po,Ri);

% Calculate Strains
B = zeros(3, 2, N_nodes);
eps = zeros(3, size(d,2), N_nodes);
sig = zeros(3, size(d,2), N_nodes);
for iE = 1:N_elements
    X = elements(iE,:);
    [~,B(:,:,iE),~] = NBJvec(1,X(1),X);
    eps(:,:,iE) = B(:,:,iE) * d([iE iE+1],:);
end
[~,B(:,:,end),~] = NBJvec(1,X(2),X);
eps(:,:,end) = B(:,:,end) * d([end-1 end],:);
for iN = 1:N_nodes
    sig(:,:,iN) = D * eps(:,:,iN);
end

time_elapsed = toc;




%% Plots

plot_d = 1;
plot_sig_rr = 1;
plot_sig_tt = 1;

if plot_d == 1
% plot displacement for inner, mid, and outer point
figure(1)
yline(0,LineWidth=1); grid on; hold on
plot(time, d(1,:),            color='#004488')
plot(time, d(round(end/2),:), color='#DDAA33')
plot(time, d(end,:),          color='#BB5566')
xlabel('Time [T]')
ylabel('Displacement [L]')
title([num2str(N_elements) ' Elements, ' ...
       algo_name ', ' ...
       '$\Delta t=$ '  num2str(tcr_ratio) '$*t_{cr}$' ' = '  num2str(dt)])
legend('','Inner Point', 'Midpoint', 'Outer Point', Location='southoutside', Orientation='horizontal')
end

if plot_sig_rr == 1
% plot radial stress for inner, mid, and outer point
figure(2)
yline(0,LineWidth=1); grid on; hold on
plot(time, sig(1,:,1),           color='#004488')
plot(time, sig(1,:,round(end/2)),color='#DDAA33')
plot(time, sig(1,:,end),         color='#BB5566')
xlabel('Time [T]')
ylabel('Radial Stress, $\sigma_{rr}$ [F/L$^2$]')
title([num2str(N_elements) ' Elements, ' ...
       algo_name ', ' ...
       '$\Delta t=$ '  num2str(tcr_ratio) '$*t_{cr}$' ' = '  num2str(dt)])
legend('','Inner Point', 'Midpoint', 'Outer Point', Location='southoutside', Orientation='horizontal')
end

if plot_sig_tt == 1
% plot radial stress for inner, mid, and outer point
figure(3)
yline(0,LineWidth=1); grid on; hold on
plot(time, sig(2,:,1),           color='#004488')
plot(time, sig(2,:,round(end/2)),color='#DDAA33')
plot(time, sig(2,:,end),         color='#BB5566')
xlabel('Time [T]')
ylabel('Circumferential Stress, $\sigma_{\theta \theta}$ [F/L$^2$]')
title([num2str(N_elements) ' Elements, ' ...
       algo_name ', ' ...
       '$\Delta t=$ '  num2str(tcr_ratio) '$*t_{cr}$' ' = '  num2str(dt)])
legend('','Inner Point', 'Midpoint', 'Outer Point', Location='southoutside', Orientation='horizontal')
end

% % Animation of displacement in cross section
% for t = 1:length(time)
%     plot(9,0); hold on
%     plot(21,0); hold off
%     for ele = 1:length(Nodal_Positions)
%         xline((Nodal_Positions(ele)+d(ele,t)),LineWidth=3)
%     end
%     pause(0.00002)
% end





%% Commandline text

disp('––––––––––––––––––––––––––––––––––––')
disp('Hollow Sphere FEA')
disp('––––––––––––––––––––––––––––––––––––')
disp(' ')
disp('Mesh Discretization:')
fprintf('    Number of Elements:       %d\n', N_elements)
fprintf('    Number of Nodes:          %d\n', N_nodes)
fprintf('    Element Size, h:          %.5g\n', L_tot/N_elements)
fprintf('    M Gauss Points:           %d\n', Ngp_M)
fprintf('    K Gauss Points:           %d\n', Ngp_K)
fprintf('    Mass Type:                %s\n', MassType)
disp(' ')
disp('Analysis:')
fprintf('    Algorithm:                %s\n', algo_name)
fprintf('    First Natural Frequency:  %.5g Hz\n', sqrt(EigMin)/2/pi)
fprintf('    Critical Time Step:       %s\n', t_cr_text)
fprintf('    Time Step Used:           %.3g * t_cr\n', tcr_ratio)
fprintf('    Computation Time:         %.3g seconds\n', time_elapsed)
disp(' ')
disp('.................')
disp('.................')
disp('.................')
disp('ANALYSIS COMPLETE')
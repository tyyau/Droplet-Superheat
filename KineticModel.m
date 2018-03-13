function [MinTempDiff, t_f, t3_f, rout] = KineticModel(fuel, W_Comp, t, r, T, P0, c)
% INPUT
% P0 - Initial ambient pressure (Pa)
% T - Temperature vector (K)

% LOAD PROPERTIES
n = size(W_Comp, 2); %  n = number of components
M_k = fuel.MW; % Kg/mol
Tc_k = fuel.TcVec; % K
Pc_k = fuel.PcVec; % Pa 
rhoc_k = fuel.rhocVec; % Kg/m^3

% TODO : Remove this. Used for Sigma
Patm = 1.01325 * 10^5; % convert bar to Pa
% TODO : Clean this
% Tsl/Tcr ratio
% Tro = ((0.11*P0./Pc_k)+0.89);
P0_ratio = min(1, P0./Pc_k);
Tro = 1-.111*(1-P0_ratio).^0.858; % Law

% Find # of points per time step
for i = 1:length(t)-1
    if t(i) ~= t(i+1)
        pts = i;
        break
    end
end

% Initialize vector arrays
x_Comp = zeros(length(t), n);
T_sl_ps = zeros(length(t), 1);
T_sl_li = zeros(length(t), 1);
T_sl_p0 = zeros(length(t), 1);
T_diff = zeros(length(t), 1);
T_sl_c = zeros(length(t), 1);
J_i = zeros(length(t), 1);
sigma = zeros(length(t), 1);
P_diff = zeros(length(t), 1);
rho = zeros(length(t), 1);
mu = zeros(length(t), 1);

% Initial temperature perturbation to onset bubble
% Depends on OpenSMOKE timestep
for i = 1:length(t)
    % Find mass fractions
    m_Comp = zeros(1, n);
    for k = 1:n
        m_Comp(k) = W_Comp(i,k) / M_k(k);
    end
    m_total = sum(m_Comp);
    x_Comp = zeros(1, n);
    for k = 1:n
        x_Comp(k) = m_Comp(k) / m_total;
    end

    % Find critical mixture temperature from Li equation
    Tc_mix = dot(Tc_k .* x_Comp, 1./rhoc_k)/ dot(x_Comp,(1./rhoc_k));
    
    % Find mixture density
    rho(i) = dot(x_Comp,M_k)./dot(x_Comp,fuel.specVol(T(i)));
    
    % Find saturated vapor pressure
    % TODO : Liquid pressure set to default ambient currently
    Pv = dot(x_Comp,fuel.Psat(T(i)));
    Pl = P0; P_diff(i) = Pv - Pl;
    
    % Find Surface Tension
    sigma(i) = dot(x_Comp, fuel.sigma(T(i)));
    
    % Find viscosity
    mu(i) = dot(x_Comp, fuel.etaL(T(i)));
   
    % Energy of critical size nucleus
    dA = 16/3*pi*sigma(i)^3/(P_diff(i))^2; % Is Pe too small?
    %dA = 16/3*pi*sigma^3/(Pv-Pl)^2;

    % Collision frequency
    k = 1.3806*10^-23;
    y_m = W_Comp(i,:) * M_k.^(1/2)';
    kf = 8*Pv*sigma(i)^2/(P_diff(i))^2*(2*pi/k/T(i))^(1/2) * y_m;
    
    % Nucleation rate
    % TODO : Validate with Glassman
    gamma = 1;
    N0 = 6.023*10^23 / dot(x_Comp,fuel.specVol(T(i)));
%     J = gamma*N0*exp(-dA/k/T(i))*kf;
%     log_J = log(gamma*N0*kf)-dA/k/T(i)
%     J = exp(log_J)
    J_i(i) = gamma*N0*exp(-dA/k/T(i))*kf;
    
    % Find superheat temperature of mixture for constant pressure
    T_sl_p0(i) = (x_Comp * Tro') * Tc_mix;
    % Find single component superheat temperature
    if c ~= 0
        N = 0;
        T_sl_c(i) = Tc_k(c).*((27/32)^1/(N+1)+P0./Pc_k(c)/((N+1)*8));
    end
    % Difference in temperature between superheat and actual
    T_diff(i) = T_sl_p0(i) - T(i);
end

t_m = reshape(t, pts, []);
r_m = reshape(r, pts, []);
T_m = reshape(T, pts, []);
T_sl_m = reshape(T_sl_p0, pts, []);
T_d = reshape(T_diff, pts, []);
T_sl_cm = reshape(T_sl_c, pts, []);
J_im = reshape(J_i, pts, []);
sigma_m = reshape(sigma, pts, []);
P_diff_m = reshape(P_diff, pts, []);
rho_m = reshape(rho, pts, []);
mu_m = reshape(mu, pts, []);

% Predict onset of bubble growth (6.80 - Zeng)
t_f = 0;
J_tot = 0;
brk = false;
index = 1;
for i = 1:size(t_m,2)-1
    dt = t_m(1,i+1) - t_m(1,i);
    for j = 1:size(t_m,1)-1
        dr = r_m(j+1,i) - r_m(j,i);
        J_tot = J_tot + (4*pi*r_m(j,i)^2*dr)*J_im(j,i)*dt;
        if J_tot > 1
            t_f = t_m(j,i);
            index = i;
            brk = true;
            break
        end
    end
    if brk
        break
    end
end
t_f
J_tot
% Critical radius
Rc = 2*sigma_m(1, index)/(P_diff_m(1,index))

% TODO : Get surface viscosity from Zeng

t_indices = t_m(1,index:end);
% Values at center
t_c = t_m(1,:);
rho_c = rho_m(1,:)
P_diff_c = P_diff_m(1,:)
sigma_c = sigma_m(1,:)
mu_c = mu_m(1,:)

Rc2 = 2*sigma_m(1, index+1)/(P_diff_m(1,index+1));
% drc_dt = 1/(t_m(1, index) - t_m(1, index-1)) * ...
%     (2*sigma_m(1, index)/(P_diff_m(1,index)) - 2*sigma_m(1, index-1)/(P_diff_m(1,index-1)));
% TODO : What is this number?
drc_dt = 0;

r_init = [Rc drc_dt];

% mu = ???

function drdt = eq(t,r)
    P_diff_l = interp1(t_c, P_diff_c, t);
    rho_l = interp1(t_c, rho_c, t);
    sigma_l = interp1(t_c, sigma_c, t);
    mu_l = interp1(t_c, mu_c, t);
    % ODE from Park 2005 -> Young 1989
    d2rdt2 = 1/r(1)*(-3/2*r(2)^2+1/rho_l*(P_diff_l - 4*mu_l/r(1)*r(2)-2*sigma_l/r(1)));
    drdt = [r(2) d2rdt2]';
end

% [t2_f, r_f] = ode45(@eq, t_indices, r_init);
% t2_f
% r_f

t_test = [t_indices(1):.00001:t_indices(end)];
[t3_f, r3_f] = ode45(@eq, t_test, r_init);

% load('Data/water_droplet_data.mat')

figure
hold on
mesh(t_m, r_m, T_m);
mesh(t_m, r_m, T_sl_m);
if c ~= 0
    mesh(t_m, r_m, T_sl_cm);
end
xlabel('time [s]');
ylabel('radius [mm]');
zlabel('temperature [K]');
scatter3(t_f, 0, T_m(1,index), 'r');

% ((1-rho_amb/rho_drop)*(ri/ro)-(ri/ro)^2)*(rho_drop*ri^3/sigma)*omega^2 ...
%     + (-1+(ri/ro)^4+rho_amb/rho_drop)*We^(1/2)*sqrt(rho_drop*ri^3/sigma)*omega ...
%     + 2*(ri/ro)^2 + 2*(ri/ro)^-2 - 3*rho_bub/rho_drop== 0; 

% omega gives k_crit
% breakup when k_crit = 5?
% legend('1','2');

figure
mesh(t_m, r_m, T_d);
xlabel('time [s]');
ylabel('radius [mm]');
zlabel('temperature [K]');

% figure
% plot(t2_f, r_f(:,1));
% xlabel('time [s]');
% ylabel('radius [mm?]');
figure
plot(t3_f, r3_f(:,1));
xlabel('time [s]');
ylabel('radius [mm?]');
rout = r3_f(:,1);
MinTempDiff = min(min(T_d))

end
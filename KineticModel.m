function [MinTempDiff] = KineticModel(Comp_Prop, Antoine, W_Comp, t, r, T, P0, c)
% W_Comp = [w_IC8H18 w_TMBENZ w_NPBENZ w_NC12H26];
% Comp_Prop = [IC8H18_Prop TMBENZ_Prop NPBENZ_Prop NC12H26_Prop];
% Comp_Prop -> (1)=M, (2)=Tc, (3)=Pc, (4)=Vc 

n = size(W_Comp, 2); %  n = number of components
M_k = Comp_Prop(1,:);
Tc_k = Comp_Prop(2,:);
Pc_k = Comp_Prop(3,:);
Vc_k = Comp_Prop(4,:);
om_k = Comp_Prop(5,:);
Tb_k = Comp_Prop(6,:);

A_k = Antoine(1,:);
B_k = Antoine(2,:);
C_k = Antoine(3,:);

M_fuel = sum(M_k);

% Tsl/Tcr ratio
% Tro = ((0.11*P0./Pc_k)+0.89);
Tro = 1-.111*(1-min(1, P0./Pc_k)).^0.858; % Law

% Find # of points per time step
for i = 1:length(t)-1
    if t(i) ~= t(i+1)
        pts = i;
        break
    end
end

% Initialize vector arrays
x_Comp = zeros(length(t), n);
T_sl_ps = zeros(length(t),1);
T_sl_li = zeros(length(t),1);
T_sl_p0 = zeros(length(t),1);
T_diff = zeros(length(t),1);
T_sl_c = zeros(length(t),1);
J_i = zeros(length(t),1);
for i = 1:length(t)
    % Find mole fractions
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
    phi_tot = x_Comp * Vc_k';
    Tc_Li = sum(Tc_k .* x_Comp .* Vc_k)/ phi_tot;
    
    % Find specific volume from Yamada-Gunn modified Rackett Equation
    v_k = Vc_k.*(0.29056-0.08775.*om_k).^(1 - T(i)./Tc_k).^(2/7);
    v_Comp = W_Comp(i,:) * v_k';
    N0 = 6.023*10^23 / v_Comp;
    
    % Find pressure from Antoine equation
    Pg_k = 10.^(A_k-B_k./(T(i)+C_k)) / 750.0616; %bar
    Pg = x_Comp * Pg_k' * 10^5;
    
    % Find surface tension from Avedisian+Glassman
    Patm = 1.01325 * 10^5;
    Pl = P0 * 10^5;
    Pc_avg = x_Comp * Pc_k' * 10^5;
    Tb_avg = x_Comp * Tb_k';
    sm = Tb_avg/(Tc_Li-Tb_avg) * log(Pc_avg/Patm);
    beta = -.4412 + .2003*sm - 0.00516*sm.^2;
    sigma = beta * Pc_avg^(2/3) * Tc_Li^(1/3) * (1-T(i)/Tc_Li)^(11/9);
    % Energy of critical size nucleus
    dA = 16/3*pi*sigma^3/(Pg-Pl)^2
    s1 = sigma^3;
    sigma;
    p1 = (Pg-Pl)^2;
    
    if T(i) > 500
          (Pg-Pl);  
    end

    % Collision frequency
    k = 1.38064852*10^-23;
    y_m = W_Comp(i,:) * M_k.^(1/2)';
    kf = 8*Pg*sigma^2/(Pg-Pl)^2*(2*pi/k/T(i))^(1/2) * y_m;
    
    % Nucleation rate
    gamma = 1;
    J = gamma*N0*exp(-dA/k/T(i))*kf;
    test5 = gamma*N0*kf
    log_J = log(gamma*N0*kf)-dA/k/T(i);
    J = exp(log_J)
%     if J > 0
%         J
%     end
    % J~10^5 nuclei/cm3-s
%     test2 = N0*kf;
    test = -dA/k/T(i)
%     test3 = dA*(k*log(N0*gamma*kf/(10^5*100^3)))^-1;
%     test4 = N0*gamma*kf/(10^3);
    % Find superheat temperature of mixture for constant pressure
    T_sl_p0(i) = (x_Comp * Tro') * Tc_Li;
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

% J_tot = 0;
% for i = 1:size(t_m,1)-1
%     dt = t_m(i+1,j) - t_m(i,j);
%     for j = 1:length(t_m,2)-1
%         dr = r_m(i,j+1) - r_m(i,j);
%         J_tot = Jtot + (4*pi*r_m(i,j)^2*dr)*J*dt;
%         if Jtot > 1
%             t_f = t_m(i,j);
%             break
%         end
%     end
% end

% figure
% hold on
% mesh(t_m, r_m, T_m);
% mesh(t_m, r_m, T_sl_m);
% if c ~= 0
%     mesh(t_m, r_m, T_sl_cm);
% end
% xlabel('time [s]');
% ylabel('radius [mm]');
% zlabel('temperature [K]');
% % legend('1','2');
% 
% figure
% mesh(t_m, r_m, T_d);
% xlabel('time [s]');
% ylabel('radius [mm]');
% zlabel('temperature [K]');

MinTempDiff = min(min(T_d))

end
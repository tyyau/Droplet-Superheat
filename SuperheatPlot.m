function [MinTempDiff] = SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, c)
% W_Comp = [w_IC8H18 w_TMBENZ w_NPBENZ w_NC12H26];
% Comp_Prop = [IC8H18_Prop TMBENZ_Prop NPBENZ_Prop NC12H26_Prop];
% Comp_Prop -> (1)=M, (2)=Tc, (3)=Pc, (4)=Vc 

n = size(W_Comp, 2); %  n = number of components
M_k = Comp_Prop(1,:);
Tc_k = Comp_Prop(2,:);
Pc_k = Comp_Prop(3,:);
Vc_k = Comp_Prop(4,:);

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
% legend('1','2');

figure
mesh(t_m, r_m, T_d);
xlabel('time [s]');
ylabel('radius [mm]');
zlabel('temperature [K]');

MinTempDiff = min(min(T_d))

end


clear

rho_al = 800;  % aluminum density kg/m3
rho = 1.225; % air density kg/m3
c = 2;        % chord length m
xf = 0.45*c;    % location of elastic axis from leading edge
thicc = 0.01;   % thickness m
L = 3;        % length m
Wnat_plun = 2*2*pi; % plunge natural frequency 2hz (rad/s)
Wnat_pith = 6*2*pi; % pitch natural frequency 6 hz (rad/s)
m = rho_al*(thicc*L)*c;  % mass kg
Kh = m*Wnat_plun^2;  % plunge stiffness
a = -0.05*c;
b = c/2;
e = xf/c - 0.25;

Ialpha = m*c^2/12 + m*a^2;
Kalpha = Ialpha*Wnat_pith^2; % pitch stiffness
S = -m*a;
% Setting up the velocity
max_velocity = 100;
vel = linspace(1,max_velocity,max_velocity);

M_pk = [m S; S Ialpha];
K_pk = [Kh 0;0 Kalpha];
DK_pk = rho*pi*b^2*[-1 a;a -a^2-((b^2)/8)];
MD_pk = M_pk - DK_pk;% 2x2
om_in = [sqrt(Kh/m) sqrt(Kalpha/Ialpha)];  % first try
epsilon = 10E-6;  % degree of error allowed
tol = 1;
freq_stored = zeros(2,max_velocity);
% conducting the newton method for iteration
for i = 1:max_velocity
    for j = 1:2
        om = om_in(j);  % first try
        while tol > epsilon
            init = om;
            k=init*b/vel(i);
            C = 1 - 0.165/(1-(0.0455/k)*1i) - 0.335/(1-(0.30/k)*1i);
            EK1 = -2*pi*rho*b*vel(i)*C;
            EK2 = -rho*pi*b^2*vel(i) + 2*pi*rho*b*vel(i)*C*(a-b/2);
            EK3 = rho*pi*b^2*vel(i) - rho*pi*b*vel(i)*(b-((2*a)+b)*C);
            EK4 = -rho*pi*b^2*vel(i)*c/4 + rho*pi*b*vel(i)*(b-((2*a)+b)*C)*(a-b/2);
            EK = [EK1 EK2;EK3 EK4];
            FK2 = -2*pi*rho*b*(vel(i)^2)*C;
            FK4 = rho*pi*b^2*(vel(i)^2) - rho*pi*b*vel(i)^2*(b-((2*a)+b)*C);
            FK = [0 FK2;0 FK4];
            FK_K = FK - K_pk; % 2x2
            A_pk = [
                zeros(2,2) eye(2)
                MD_pk\FK_K MD_pk\EK]; % 4x4
            [evece_pk, cv_pk]=eig(A_pk);
            cv_pk = diag(cv_pk);
            [sorted,index]=sort(imag(cv_pk),"descend");
            cv_pkkkk = cv_pk(index);
            cv_pk_om = imag(cv_pk(index));
            om = cv_pk_om(j);   %update the value of omega
            tol = abs(init-om);
            freq_stored(j,i) = cv_pkkkk(j);  % continuously update the variable until stop
        end  
        if tol < epsilon
            tol = 1;   % reset the value of tol once while condition is reached
        end
    end
end
freq_pk = imag(freq_stored)/2/pi;
damp_pk = zeros(2,max_velocity);
for i = 1:max_velocity
    damp_pk(1,i) = - real(freq_stored(1,i))/(freq_pk(1,i));
    damp_pk(2,i) = - real(freq_stored(2,i))/(freq_pk(2,i));
end

figure
three_pk = plot(vel,damp_pk(1,:), "Marker","o");
hold on
four_pk = plot(vel,damp_pk(2,:), "Marker","o");
two_un = plot([0 max_velocity], [0 0],'Color','black');
xlim([0 100])
ylabel('$ ~\mathrm{\zeta}$', 'interpreter','latex', "FontSize",14);
xlabel('$ ~\mathrm{U (m/s)}$', 'interpreter','latex', "FontSize",14);
legend([three_pk four_pk], {'Pitch','Plunge'},...
    'interpreter','latex', "Location","northwest");
title("p-k Method Plot of Damping Ratios against Velocity (Unsteady Case)",'Interpreter',"latex","FontSize",14)
annotation('textarrow',[0.7897 0.8516],[0.2304 0.2601],'String','PK Flutter Speed: 61.66 m/s')
hold off

% Plot of Frequency Ratios against Velocity
figure
one_pk = plot(vel,freq_pk(1,:), "Marker","o");
hold on
two_pk = plot(vel,freq_pk(2,:),"Marker","o");
flut_un_pk = plot([61.66 61.66], [1.5 6.5], 'LineStyle',"-.",'Color','black');
ylabel('$ ~\mathrm{\omega_n (Hz)}$', 'interpreter','latex', "FontSize",14);
xlabel('$ ~\mathrm{U (m/s)}$', 'interpreter','latex', "FontSize",14);
legend([one_pk two_pk], {'Pitch','Plunge'},...
    'interpreter','latex', "Location","west");
title("p-k Method Plot of Natural Frequencies against Velocity (Unsteady Case)",'Interpreter',"latex","FontSize",14)
xlim([0 100])
annotation('textarrow',[0.7881 0.8373],[0.3325 0.3283],'String','2.682')
annotation('textarrow',[0.804 0.8437],[0.3918 0.3854],'String','3.14')
hold off

% Plot Root Locus
figure
grid on
hold on
for i = 1:length(vel)
    % Extract the real and imaginary parts of the eigenvalues for each velocity
    real_part = real(freq(:, i));
    imag_part = imag(freq(:, i));
    
    % Plot the real vs. imaginary parts for each eigenvalue at each velocity
    plot(real_part(1), imag_part(1), 'bo'); % First eigenvalue
    plot(real_part(2), imag_part(2), 'ro'); % Second eigenvalue
end

% Labeling the plot
xlabel('$ ~\mathrm{Re(\omega)}$', 'interpreter','latex', "FontSize",11);
ylabel('$ ~\mathrm{Im(\omega)}$', 'interpreter','latex', "FontSize",11);
title("Root Locus Plot", 'Interpreter', "latex", "FontSize", 11)
legend({'Largest Eigenvalue', 'Second Largest Eigenvalue'}, 'Location', 'northwest', 'interpreter', 'latex')
hold off

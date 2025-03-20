clear

% Define Key Terms
rho = 1.225; % Air Density

mu = 800; % Density of Airfoil

c = 2;
x_f = 0.45*c; % Location of Elastic Axis from LE
a = -0.05*c;
b = c/2;
thickness = 0.01; % thickness m
L = 3; % length m
mass = mu*thickness*c; % mass kg

omega_h = 2*2*pi; % plunge (rad/s)
omega_alpha = 6*2*pi; % pitch (rad/s)
K_h = mass*((omega_h)^2); 
I_alpha = (mass*(c^2))/12 + mass*(a^2);
K_alpha = (I_alpha)*(omega_alpha^2); % pitch stiffness
S = -mass*a ;

% Define Matrix
M = [mass S; S I_alpha];
K = [K_h 0;0 K_alpha];
D = rho*pi*b^2*[-1 a;a -(a^2)-((b^2)/8)];
m_d = M - D;

first_try = [sqrt(K_h/mass) ; sqrt(K_alpha/I_alpha)]; % Initial guess for omega and k
epsilon = 10E-5 ; % degree of error allowed


% Setting up the velocity
vel = linspace(1,100,100);
freq = zeros(2,100);


for i = 1:100 % Iterate airspeed from 1 to 100
    for j = 1:2 % Degree of Freedom, j = 1 is plunge, j = 2 is pitch
        omega = first_try(j); % first try
        error = 1;

        while error > epsilon
             init_guess = omega;
             k=init_guess * b/ vel(i);
             C = 1 - 0.165/(1-(0.0455/k)*1i) - 0.335/(1-(0.30/k)*1i);
               
             % Define E
             E_1 = -2*pi*rho*b*vel(i)*C;
             E_2 = -rho*pi*b^2*vel(i) + 2*pi*rho*b*vel(i)*C*(a-b/2);
             E_3 = rho*pi*b^2*vel(i) - rho*pi*b*vel(i)*(b-((2*a)+b)*C);
             E_4 = -rho*pi*b^2*vel(i)*c/4 + rho*pi*b*vel(i)*(b-((2*a)+b)*C)*(a-b/2);

             E = [E_1 E_2;E_3 E_4];

             % Define F
             F_2 = -2*pi*rho*b*(vel(i)^2)*C;
             F_4 = rho*pi*b^2*(vel(i)^2) - rho*pi*b*vel(i)^2*(b-((2*a)+b)*C);
             F = [0 F_2;0 F_4];

             F_K = F - K;
             A = [zeros(2,2) eye(2); m_d\F_K m_d\E]; % 4x4
             
             [~, A2]=eig(A); 
             A2 = diag(A2);
             [sorted,index]=sort(imag(A2),"descend");
             omega = sorted(j);
             
             error = abs(init_guess-omega);
             freq(j,i) = A2(index(j)); % continuously updates Frequency
         end 
         
         if error < epsilon
            error = 1; % reset the value of tol once while condition is reached

         end
     end
end

freq_pk = abs(freq)/2/pi;
damp_pk = -real(freq) ./ freq_pk ;

% Plot Damping Ratios
figure
grid on
hold on
plot_1 = plot(vel,damp_pk(1,:));
plot_2 = plot(vel,damp_pk(2,:));
Vel_plot = plot([0 100], [0 0],'k');
xlim([0 100])
xlabel('$ ~\mathrm{U (m/s)}$', 'interpreter','latex', "FontSize",11);
ylabel('$ ~\mathrm{\zeta}$', 'interpreter','latex', "FontSize",11);
legend([plot_1 plot_2], {'Pitch','Plunge'},'interpreter','latex', "Location","northwest");
title("Damping Ratios against Velocity",'Interpreter',"latex","FontSize",11)
annotation('textarrow',[0.736904761904762 0.627380952380952],[0.421222222222222 0.295238095238095],'String',{'P-k Flutter Speed = 64.285 m/s'},'FontSize', 10);
hold off

% Plot Frequency Ratios
figure
grid on
hold on
plot_3 = plot(vel,freq_pk(1,:));
plot_4 = plot(vel,freq_pk(2,:));
flut_speed = plot([64.285 64.285], [0 10], 'LineStyle',"-.",'Color','black');
xlim([0 100])
xlabel('$ ~\mathrm{U (m/s)}$', 'interpreter','latex', "FontSize",11);
ylabel('$ ~\mathrm{\omega_n (Hz)}$', 'interpreter','latex', "FontSize",11);
legend([plot_3 plot_4], {'Pitch','Plunge'},'interpreter','latex', "Location","west");
title("Natural Frequencies against Velocity",'Interpreter',"latex","FontSize",11)
% Create textarrow
annotation('textarrow',[0.70859375 0.630859375],...
    [0.364045806906272 0.307963354474982],'String',{'3.622 Hz'});
annotation('textarrow',[0.76796875 0.63125],...
    [0.504285412262157 0.540521494009866],'String',{'7.959 Hz'});

hold off

% Plot Root Locus
figure
grid on
hold on
plot_5 = plot(real(freq(1, :)),imag(freq(1, :) ), 'b', 'DisplayName', 'Pitch (Positive)');
plot_6 = plot(real(freq(2, :)), imag(freq(2, :)), 'r', 'DisplayName', 'Plunge (Positive)');
plot_7 = plot(real(freq(1, :)),-1*imag(freq(1, :) ), 'g', 'DisplayName', 'Pitch (Negative)');
plot_8 = plot(real(freq(2, :)),-1*imag(freq(2, :)), 'k', 'DisplayName', 'Plunge (Negative)');
xlabel('Re(z)', 'FontSize', 11);
ylabel('Im(z)', 'FontSize', 11);
legend('Location', 'northwest');
title('Root Locus Plot', 'Interpreter', 'latex', 'FontSize', 11);
hold off


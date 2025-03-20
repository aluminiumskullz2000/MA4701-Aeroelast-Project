clear all
close all

% Constants
L = 30;
rho = 1600;
E = 2*10^11;

n = 10;

syms y

% Define A and Second moment of Area
A = -60.78*(y/L)^3 + 121.16*(y/L)^2-79.624*(y/L) + 19.344 ;
I = 2.7462*exp(-7.958*(y/L));

% Define Assumed Functions
func1 = sym(zeros(1,n));
func_d = sym(zeros(1,n));
func_dd = sym(zeros(1,n));

for i = 1 :n
    func1(i) = (((y/L-1)/(i^2*pi^2))*cos(i*pi*y/L)-(2/(i^3*pi^3))*sin(i*pi*y/L)+y/(i^2*pi^2*L)+1/(i^2*pi^2));
    func_d(i) = diff (func1(i) , y);
    func_dd(i) = diff (func_d(i) , y);
end 

% Define Mass and Stiffness Matrix
M1 = zeros(n,n);
M2 = zeros(n,n);
M = zeros(n,n);
K = zeros(n,n);

for i = 1 : n
    for j= 1 : n
        % Calculate M1
        m1 = rho .* A .* func1(i) .* func1(j);
        M1(i,j) = int(m1, y, [0 L]);
        
        % Calculate M2
        M2(i,j) = 6120 * (((6.845/L-1)/(i^2*pi^2))*cos(i*pi*6.845 /L)-(2/(i^3*pi^3))*sin(i*pi*6.845 /L)+6.845 /(i^2*pi^2*L)+1/(i^2*pi^2)) * (((6.845/L-1)/(j^2*pi^2))*cos(j*pi*6.845 /L)-(2/(j^3*pi^3))*sin(j*pi*6.845 /L)+6.845 /(j^2*pi^2*L)+1/(j^2*pi^2)) ;
        
        % Total mass matrix M
        M(i,j) = double(M1(i,j) + M2(i,j));
        
        % Calculate Stiffness Matrix
        k1= E * I * func_dd(i) * func_dd(j);
        K(i,j) =  double(int(k1, y, [0 L]));
    end
end

[V, d] = eig(K, M);
        
% Sort eigenvalues and eigenvectors
[d_sorted, Index] = sort(diag(d), 'ascend');
for j = 1:n
    V_sort(:,j)=V(:,Index(j));
end 

% Obtain natural frequencies
nat_freq = zeros(1,n);
for k = 1:n
    nat_freq(1,k) = real(sqrt(d_sorted(k))/2/pi);
end
nat_freq;
fprintf(['First Four Natural Frequencies (Hz), ...' ...
    'overall: \n %9.2E \n %9.2E \n %9.2E \n %9.2E\n'],...
    nat_freq(1),nat_freq(2),nat_freq(3), nat_freq(4));


% Calculate and plot each mode shape
for j = 1:4
    mode_shape = 0;
    for i = 1:n
        mode_shape = mode_shape + (((y/L-1)/(i^2*pi^2))*cos(i*pi*y/L)-(2/(i^3*pi^3))*sin(i*pi*y/L)+y/(i^2*pi^2*L)+1/(i^2*pi^2)) .*V_sort(i,j);
    end
    
    figure(1)
    hold on
    mode_shape = mode_shape ./ subs(mode_shape,L);
    fplot(mode_shape,[0,L])
end

title('Mode Shapes')
xlabel('Position along the length of Wing (m)')
ylabel('$ ~\mathrm{\phi}$', 'interpreter','latex')
legend('Mode 1', 'Mode 2', 'Mode 3', 'Mode 4')
hold off;
clear all
close all

% Constants
L = 30;
rho = 1600;
E = 2*10^11;
y_engine = 6.845;

n = 7 ;

syms y

% Define A and Second moment of Area
A = -60.78*(y/L)^3 + 121.16*(y/L)^2-79.624*(y/L) + 19.344 ;
I = 2.7462*exp(-7.958*(y/L));

% Define Assumed Functions
func1 = sym(zeros(1,n));
func_d = sym(zeros(1,n));
func_dd = sym(zeros(1,n));

for i = 1 :n
    func1(i) = ((y ./ L).^(i + 1)).*(2 + i - i.*(y/L)) .* (1 / i / (i+1) / (i+2));
    func_d(i) = diff (func1(i) , y);
    func_dd(i) = diff (func_d(i) , y);
end 

% Define Mass and Stiffness Matrix
M1 = zeros(n,n);
M2 = zeros(n,n);
M = zeros(n,n);
K = zeros(n,n);
m = zeros (n,n);


for i = 1 : n
    for j= 1 : n
        % Calculate M1
        m1 = rho .* A .* func1(i) .* func1(j);
        M1(i,j) = int(m1, y, [0 L]);
       
        % Calculate M2
        M2(i,j) = 6120 * ((6.845/L).^(i + 1)).*(2+i-i.*(6.845/L)) .* (1 / i / (i+1) / (i+2)).*((6.845/L).^(j+1)).*(2+j-j.*(6.845/L)).* (1 / j / (j+1) / (j+2)) ;
        
        % Total mass matrix M
        M(i,j) = double(M1(i,j) + M2(i,j));
        m(i,j) = double (M1(i,j)) ;
        
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

[V2, d2] = eig(K, m);
        
% Sort eigenvalues and eigenvectors
[d_sorted2, Index2] = sort(diag(d2), 'ascend');
for j = 1:n
    V_sort2(:,j)=V2(:,Index2(j));
end 

difference = V_sort - V_sort2

% Obtain natural frequencies
nat_freq = zeros(1,n);
for k = 1:n
    nat_freq(1,k) = real(sqrt(d_sorted(k))/2/pi);
end
nat_freq;
fprintf(['First Four Natural Frequencies (Hz) with engine, ...' ...
    'overall: \n %9.2E \n %9.2E \n %9.2E \n %9.2E\n'],...
    nat_freq(1),nat_freq(2),nat_freq(3), nat_freq(4));

% Obtain natural frequencies
nat_freq2 = zeros(1,n);
for k = 1:n
    nat_freq2(1,k) = real(sqrt(d_sorted2(k))/2/pi);
end
nat_freq;
fprintf(['First Four Natural Frequencies (Hz) without engine, ...' ...
    'overall: \n %9.2E \n %9.2E \n %9.2E \n %9.2E\n'],...
    nat_freq2(1),nat_freq2(2),nat_freq2(3), nat_freq2(4));

% Calculate and plot each mode shape
for j = 1:4
    mode_shape = 0;
    mode_shape2 = 0;
    for i = 1:n
        mode_shape = mode_shape + ((y/L)^(i+1))*(2+i-i*(y/L)) / (i*(i+1)*(i+2)) .*V_sort(i,j) ;
        mode_shape2 = mode_shape2 + ((y/L)^(i+1))*(2+i-i*(y/L)) / (i*(i+1)*(i+2)) .*V_sort2(i,j); 
    end
    
    % Normalize mode shapes
    mode_shape = mode_shape ./ subs(mode_shape,L);
    mode_shape2 = mode_shape2 ./ subs(mode_shape2,L);
    
    % Create a new figure for each mode
    figure;
    
    % Plot mode shape with engine
    fplot(matlabFunction(mode_shape), [0, L],'-b','LineWidth', 1); % Blue for with engine
    hold on;
    
    % Plot mode shape without engine
    fplot(matlabFunction(mode_shape2), [0, L], "--r", 'LineWidth', 1); % Red for without engine
    
    % Formatting the plot
    title(['Mode Shape ' num2str(j)]);
    xlabel('Position along the length of Wing (m)');
    ylabel('$ ~\mathrm{\phi}$', 'interpreter','latex');
    legend('With Engine', 'Without Engine');
    grid on;
    
    % Release the hold for the next figure
    hold off;
end



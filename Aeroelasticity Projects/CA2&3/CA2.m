clear all
close all

l = 30;
rho = 1600;
E = 2*10^11;

n = 5;

for i = 1 : n
    for j= 1 : n
        syms y;
        A = -60.78*(y/l)^3 + 121.16*(y/l)^2-79.624*(y/l)+19.344;
        I = 2.7462*exp(-7.958*(y/l));

        phi_i = ((y/l)^(i+1))*(2+i-i*(y/l))/(i*(i+1)*(i+2));
        phi_j = ((y/l)^(j+1))*(2+j-j*(y/l))/(j*(j+1)*(j+2));
        
        % Calculate M1
        m1 = rho * A * phi_i * phi_j;
        M1(i,j) = eval(int(m1, 0, l));
        
        % Calculate M2
        M2(i,j) = 6120 * ((6.845/l)^(i+1)) * (2 + i - i*(6.845/l)) / (i*(i + 1*(i + 2)))^2;
        
        % Total mass matrix M
        M(i,j) = M1(i,j) + M2(i,j);
        
        k1= E * I * diff(phi_i,2) * diff(phi_j,2);
        K(i,j) =  eval(int (k1,0,l));

        [V, d] = eig(M\K);
        
        % Sort eigenvalues and eigenvectors

    end
end

[d, Index] = sort(diag(d), 'ascend');
V_sort(:,1:n) = V(:,Index);

% Plot the mode shapes
figure(1)
hold on

% Numeric y for plotting
y_numeric = linspace(0, l, 100);

% Multiply Assumed Shapes with EigenVector
for j = 1:n
    Bigphi = 0;

    for i = 1 : 5
        syms y
        phi_i = ((y/l)^(i+1))*(2+i-i*(y/l))/(i*(i+1)*(i+2));
        phii= phi_i*V_sort(i,j);
        Bigphi = Bigphi + phii;
    end

    figure(1)
    hold on
    Bigphi = Bigphi./subs(Bigphi,l)
    fplot(Bigphi,[0,l],'LineWidth', 1.5)
    
    % Convert Bigphi from symbolic to numeric for plotting
    Bigphi_numeric = double(subs(Bigphi, y, y_numeric));
    
    % Normalize and plot the mode shape
    Bigphi_numeric = Bigphi_numeric / max(abs(Bigphi_numeric));  % Normalize to maximum value
    plot(y_numeric, Bigphi_numeric, 'LineWidth', 1.5)
end 

title('Mode Shapes')
xlabel('Position along the length')
ylabel('Amplitude')
legend('Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5')
hold off
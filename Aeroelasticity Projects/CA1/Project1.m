clear all
close all

% Defining Variables
k_1 = 750;
k_2 = 1500;
c = 1;
a = 0.1;
m = 100;
m_p2 = 10;
S = -m * a;
I_a = 1/12*m*c^2 + m*a^2; 

i = 1; 
for p = -c/2 : c/100 : c/2 %iterating through values of p with intervals of 0.01
       j = 1;
       for m_p1 = 0 : 1 : 35 %iterating though mass of m_p1 with intervals of 1

           K = [k_1, 0; 0, k_2]; %stiffness matrix

           A = m + m_p2 + m_p1;
           B = S - m_p2*(a-c/2) - m_p1*(a-p);
           C = S - m_p2*(a-c/2) - m_p1*(a-p);
           D = I_a + m_p2*((a-c/2)^2) + m_p1*((a-p)^2);

           M = [A B; C D];%mass matrix

           omega_sq = eig(K,M);%solve for eigenvalue


           omega = sqrt(omega_sq)/2/pi ;
           
           freq_heave(i,j) = omega(1);
           freq_pitch(i,j) = omega(2);
           p_store(i,j) = p;
           m_store(i,j) = m_p1;
           
           j=j+1;

       end
       i=i+1;
end 

%Plot Heave Frequency
figure
mesh(p_store,m_store,freq_heave)
xlabel('p')
ylabel('mass')
zlabel('freq(Hz)')
title('Frequency(Heave)')

%Plot Pitch Frequency
figure
mesh(p_store,m_store,freq_pitch)
xlabel('p')
ylabel('mass')
zlabel('freq(Hz)')
title('Frequency(Pitch)')



        

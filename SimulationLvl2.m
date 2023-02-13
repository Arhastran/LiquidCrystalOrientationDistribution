%Simulation with electric field 13.02.2023

n_0 = 1.5;
n_e = 1.65;

m = 20;
n = 10;

boundaries = [pi/2, pi/2, 0, 0];

e_0 = 8.8542*10^-12;
e_o = 5.1;
e_e = 19.6;
delta_e = e_e-e_o;
k = 20*10^-12; 
E_x1 = 10^9;e
E_x2 = 10^9*0.10;
E_x3 = 10^9*0.26;

cons = (e_0*delta_e/(4*k))*abs(E_x1)^2;
onssss = (e_0*delta_e/(4*k));

x = SOR_LC_Lvl2(m,n,boundaries,e_0,delta_e,k,E_x2);
x = x.*57.324840764331;

neff = @(phi) n_0*n_e*(1/sqrt(n_0^2*sin(phi)^2+n_e^2*cos(phi)^2));
N = arrayfun(neff, x);

figure(Color='w');
% imagesc(linspace(0,n,n*2),linspace(0,m,m*2),x)
% xlabel([num2str(n) ' [\mu' 'm]']); ylabel([num2str(m) ' [\mu' 'm]']);
% set(gca,'YDir','normal')
% c = colorbar;
% c.Label.String = 'Angles distribution [\circ]';
subplot 121; imagesc(linspace(0,n,n*2),linspace(0,m,m*2),x); title("E ="+string(E_x2))
xlabel([num2str(n) ' [\mu' 'm]']); ylabel([num2str(m) ' [\mu' 'm]']);
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = 'Angles distribution [\circ]';
y = SOR_LC_Lvl2(m,n,boundaries,e_0,delta_e,k,E_x3);
y = y.*57.324840764331;
subplot 122; imagesc(linspace(0,n,n*2),linspace(0,m,m*2),y); title("E = "+string(E_x3))
xlabel([num2str(n) ' [\mu' 'm]']); ylabel([num2str(m) ' [\mu' 'm]']);
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = 'Angles distribution [\circ]';

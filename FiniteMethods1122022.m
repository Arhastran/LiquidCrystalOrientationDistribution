n_0 = 1.5;
n_e = 1.65;

m = 10;
n = 5;

boundaries = [pi/2, pi/2, 0, 0];


x = SOR_LC(m,n,boundaries);
x = x.*57.324840764331;

neff = @(phi) n_0*n_e*(1/sqrt(n_0^2*sin(phi)^2+n_e^2*cos(phi)^2));
N = arrayfun(neff, x);

imagesc(linspace(0,n,n*2),linspace(0,m,m*2),x)
xlabel([num2str(n) ' [\mu' 'm]']); ylabel([num2str(m) ' [\mu' 'm]']);
set(gca,'YDir','normal')
c = colorbar;

c.Label.String = 'Angles distribution [\circ]';



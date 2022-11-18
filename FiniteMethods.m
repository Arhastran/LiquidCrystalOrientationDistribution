%Finite methods for laplacian phi == 0

n_0 = 1.5;
n_e = 1.65;

n = 100;
e = ones(n,1);
delta_x = 0.01;
A = spdiags([-1*e 0*e e], [-1 0 1], n, n);
A(1,1)==0;
A(end,end)==0;
initial_guess = ones(n);
%initial_guess(1) == pi;
%initial_guess(end) == pi;
sol = zeros(n);
omega = 0.5;
iter = 1000;
L = (1/(2*delta_x)).*(A*A);

% We are looking for x in A*x=sol, and initial_guess is first guess for x

solution = SOR(L,initial_guess,sol,omega,iter,0.001).*pi;

figure(Color='w')
a1 = subplot(1,2,1); imagesc(e,e,solution); title("Angle distribution");
colorbar;
colormap(a1, "winter")


neff = @(phi) n_0*n_e*(1/sqrt(n_0^2*sin(phi)^2+n_e^2*cos(phi)^2));
N = arrayfun(neff, solution);
a2 = subplot(1,2,2); imagesc(e,e,N); title("RI distribution");
colormap(a2, "bone")
colorbar;

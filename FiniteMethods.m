%Finite methods for laplacian phi == 0


n = 100;
e = ones(n,1);
delta_x = 0.01
A = spdiags([-1*e*(1/(2*delta_x)) 0*e e*(1/(2*delta_x))], [-1 0 1], n, n);
A(1,1)==0
A(end,end)==0
initial_guess = ones(n);
sol = zeros(n);
omega = 0.5;
iter = 1000;

% We are looking for x in A*x=sol, and initial_guess is first guess for x

solution = SOR(A,initial_guess,sol,omega,iter,1)



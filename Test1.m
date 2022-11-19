%% Second dim tests

n_0 = 1.5;
n_e = 1.65;

n = 100;
m = 100;
mn = m*n;
ex = ones(n,1);
ey = ones(m,1);
A = spdiags([ex -2*ex ex], [-1 0 1], n, n);
B = spdiags([ey -2*ey ey], [-1 0 1], m, m);
%L = kron(A, speye(n)) + kron(speye(m), B); 

s = size(A);
initial_guess = ones(s(1),1);
initial_guess(1) = pi/2;
initial_guess(end) = pi/2;
%initial_guess(1,:) = pi/2;
%initial_guess(end,:) = pi/2;
sol = zeros(s(1),1);
omega = 0.1;
iter = 5000;

solution = SOR(A,initial_guess,sol,omega,iter,0.001);
%imagesc(solution);
initial_guess2 = [];
for i=1:m
    initial_guess2 = [initial_guess2, solution];
end
initial_guess2(:,1) = 0;
initial_guess2(:,end) = 0;
initial_guess2(1,:) = pi/2;
initial_guess2(end,:) = pi/2;
initial_guess2(1,1) = 0;
initial_guess2(1,end) = 0;
initial_guess2(end,1) = 0;
initial_guess2(end,end) = 0;
sol2 = zeros(size(initial_guess2));

solution2 = SOR(B,initial_guess2,sol2,omega,iter,0.001);
imagesc(solution2);
colormap("jet")
colorbar

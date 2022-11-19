%Finite methods for laplacian phi == 0

n_0 = 1.5;
n_e = 1.65;

n = 500;
e = ones(n,1);
delta_x = 0.001;
A = spdiags([e -2*e e], [-1 0 1], n, n);
%A(1,1)=0;
%A(1,2)=0
%A(end,end)=0;
%A(end,end-1)=0;
initial_guess = ones(n,1);
initial_guess(1) = pi/2;
initial_guess(end) = pi/2;
sol = zeros(n,1);
omega = 0.5;
iter = 5000;
%L = (A*A);

%L(1,1)==0;
%L(end,end)==0;

% We are looking for x in A*x=sol, and initial_guess is first guess for x

solution = SOR(A,initial_guess,sol,omega,iter,0.00001).*pi/2;

figure(Color='w')
colormap jet;
imagesc(e,e,solution); title("Angle distribution");
colorbar;
%colormap(a1, "winter");

figure(Color='w')
neff = @(phi) n_0*n_e*(1/sqrt(n_0^2*sin(phi)^2+n_e^2*cos(phi)^2));
N = arrayfun(neff, solution);
imagesc(e,e,N); title("RI distribution");
colormap bone;
colorbar;
% 
 Secondim = repmat(A,50,50);
 ss = size(Secondim)
 DD = spdiags(Secondim,[-1,0,1],ss(1),ss(1));
 ACell = repmat({A}, 1, 10);
 BigA = blkdiag(ACell{:});
 s = size(BigA);
 sol2 = zeros(s(1),1);
 initial_guess2 = repmat(solution,10,1); %ones(s(1));
 %initial_guess2 = ones(s(1),1);
 %initial_guess2(1) == 0;
 %initial_guess2(end) == 0;
 solution2dim = SOR(BigA,initial_guess2,sol2,omega,1000,0.001).*pi;
% 
 figure(Color='w');
 b1 = imagesc(solution2dim);
 colormap bone;
 colorbar;
 %colormap(b1, "bone");










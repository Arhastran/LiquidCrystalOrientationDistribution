function [x] = SOR_LC_Lvl2(m,n,boundaries,e_0,delta_e,k,E)
% SOR method for single mxn matrix, dedicated to solving distribution of angle in LC cell with mxn size in um
% Boundaries should come in form [a, b, c, d], and they represent edges of cell
% This function takes into account electric field 

con = 5;
ex = ones(m+con,1);
ey = ones(n+con,1);
A = spdiags([ex -2*ex ex], [-1 0 1], m+con, n+con);
B = spdiags([ey -2*ey ey], [-1 0 1], m+con, n+con);
L = kron(A, speye(m+con,n+con)) + kron(speye(m+con,n+con), B); 

constant = (e_0*delta_e/4*k)*abs(E)^2;

initial_guess = ones(size(L));
initial_guess(:,1) = boundaries(1);
initial_guess(:,end) = boundaries(2);
initial_guess(1,:) = boundaries(3);
initial_guess(end,:) = boundaries(4);

S = size(L) ;
sol = zeros(S(1),1);
omega = 0.1;
iter = 5000;
tol = 1E-5;
x = initial_guess;
x_init = initial_guess;
counter = 0;
while(counter<iter)
    for i=2:1:S(1)-1
             for j=2:1:S(2)-1
                 x(i,j)= (1/4)*(x(i,j+1)+x(i-1,j)+x(i,j-1)+x(i+1,j))+constant*sin(2*x(i,j));
             end
         end
    if norm(x-x_init)<tol
        break;
    end
    x_init=x;
    counter = counter+1;
end

end
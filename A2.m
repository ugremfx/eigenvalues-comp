

h = 0.01; % Distance between 2 grid points

% Define V(z)
z = -3:h:3; % The 1 dimensional grid with step size h
V = z.^2/2; % Potential. Outside of the interval is considered Infinite. This condition is eforced through Psi=0 outside of the box. For a non-infinite potential outside of the interval one must resolve the Schrödinger equation outside as well..
Vbump = V + 2*exp(-10*z.^2); % Potential with bump

% Plot the potential
figure
plot(z,V)
hold on
plot(z,Vbump)
xlim([min(z),max(z)]);
xlabel('z')
ylabel('Potential V(z)')
legend('Potential without bump','Potential with bump')
title('Potential energy')

%% Define the banded H matrix, i.e. construct the Hamiltonian matrix for the
% one dimensional Harmonic oscillator

no_gridpt = length(z); % determine the grid size
H = zeros(no_gridpt); % pre allocation
Hbump = H; % The hamiltonian for the bump potential


c0h =  -1/(12*h^2)*30; % Coefficients for 5 point second order derivative
ch =  1/(12*h^2)*16;
c2h =  -1/(12*h^2)*1;
cmh =  1/(12*h^2)*16; %c_-h , coefficient from the five point derivative approximation
cm2h =  -1/(12*h^2)*1;

for ii = 1:(no_gridpt-2)
 H(ii,ii) = -1/2*c0h + V(ii); % Main diagonal
 H(ii,ii+1) = -1/2*ch; % First super diagonal
 H(ii,ii+2) = -1/2*c2h; % Second super diagonal
 
 Hbump(ii,ii) = -1/2*c0h + Vbump(ii); % Main diagonal
 Hbump(ii,ii+1) = -1/2*ch; % First super diagonal
 Hbump(ii,ii+2) = -1/2*c2h; % Second super diagonal
end

H(end-1,end-1) = -1/2*c0h + V(end-1);
H(end-1,end  ) = -1/2*ch;
H(end,end) = -1/2*c0h + V(end);

H = tril(H.') + triu(H,1);  % Takes upper half of H to make H symmetric

Hbump(end-1,end-1) = -1/2*c0h + Vbump(end-1);
Hbump(end-1,end  ) = -1/2*ch;
Hbump(end,end) = -1/2*c0h + Vbump(end);

Hbump = tril(Hbump.') + triu(Hbump,1);  % Takes upper half of H to make H symmetric


%% Find eigenvectors and eigenvalues 
tic;
[X,D] = eig(H); % returns diagonal matrix D of eigenvalues and matrix X whose columns are the corresponding right eigenvectors, so that H*X = X*D.
toc;
fprintf('Time to find the eigenvectors with Matlab is %5.3f seconds \n',toc)

tic;
[X2,D2] = ipiws(H,0);
[X3,D3] = ipiws(H,1);
[X4,D4] = ipiws(H,2);
%[X5,D5] = ipiws(H,300);
%[X6,D6] = ipiws(H,500);
toc;
fprintf('Time to find the first 3 eigenvectors with Inverse power method with shift is %5.3f seconds \n',toc)


[Xbump,Dbump] = eig(Hbump);

%%
% Check if the first 2 wave functions are orthogonal for the potential
% without bump
orth = (X(:,1)')*X(:,2);
disp(orth) 
fprintf('The first two wave functions are orthogonal because psi_i^T*psi_j = %20.16f is below the machine precision. \n',orth)
% They are also normalized because the norm of the wavevector is 1.


%%
figure
plot(z,X(:,1),'Linewidth',1)
hold on
plot(z,-X(:,2),'Linewidth',1)
plot(z,-Xbump(:,1),'--','Linewidth',1)
plot(z,-Xbump(:,2),'--','Linewidth',1)
%plot(z,-sqrt(2)*sin(pi*z)/norm(sqrt(2)*sin(pi*z)),'--g') % Analytical
%solution for square infinite potential
plot(z,exp(-1/2*z.^2)*max(X(:,1)),':g')
plot(z,-2*z.*exp(-1/2*z.^2)/max(2*z.*exp(-1/2*z.^2))*max(X(:,2)),':k')
xlabel('z')
ylabel('\Psi(z)')
legend(['Numerical Lowest Energy level without bump $\frac{E}{\hbar\omega}$ = ',num2str(D(1,1),'%.3f')],['Numerical Second Energy level without bump $\frac{E}{\hbar\omega}$ = ',num2str(D(2,2),'%.3f')],['Numerical Lowest Energy level with bump $\frac{E}{\hbar\omega}$ = ',num2str(Dbump(1,1),'%.3f')],['Numerical Second Energy level with bump $\frac{E}{\hbar\omega}$ = ',num2str(Dbump(2,2),'%.3f')],'Analytical Ground state solution for V(z) = $z^2/2$','Analytical Second wavefunction solution for $V(z) = z^2/2$','Interpreter','latex','location','Best')
title('Wavefunctions numerical vs. analytic solutions')

%%
fprintf('The overlap integral for the ground states between the potential with bump and without is %8.3f . \n',(X(:,1)')*Xbump(:,1))

%%
%Energy levels
Enlevel = diag(D);
disp(Enlevel(1:10))


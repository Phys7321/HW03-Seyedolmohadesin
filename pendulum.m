function [period,sol] = pendulum(R,theta0,thetad0,grph) 
%Modify: Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0
    error('Must input length and initial conditions')
end
if nargin==1
   theta0 = pi/2;
   thetad0=0;
   grph=0;
end
if nargin==2
    thetad0 = 0;
    grph=1;
end
if nargin==3
    grph=1;
end
g=9.81;
omega = sqrt(g/R);
T= 2*pi/omega;
% number of oscillations to graph
N = 10;


tspan = [0 N*T];
%opts = odeset('events',@events,'refine',6); %Here for future event finder
opts = odeset('refine',6);
r0 = [theta0 thetad0];
[t,w] = ode45(@proj,tspan,r0,opts,g,R);
sol = [t,w];
ind= find(w(:,2).*circshift(w(:,2), [-1 0]) <= 0);
ind = chop(ind,4);
period= 2*mean(diff(t(ind)));
 
E0 = 0.5 .* R.^2 .* thetad0 .^2 + g .* R .* (1-cos(theta0));  % Initial energy of the pendulum

KE = 0.5 .* R.^2 .* w(ind(1): ind(3),2).^2 ; % Kinetic Energy
U = g .* R .* (1-cos(w(ind(1): ind(3),1))) ; % Potential Energy
E=KE+U; % Total Energy in one cycle

deltaE= (E(:) - E0) ./ E0; % function delta 

% Small-angle approximation solution
delta = atan(theta0/(omega*thetad0));
y = theta0*sin(omega*t+delta);

if grph % Plot Solutions of exact and small angle
    f1=figure;
    
    subplot(2,1,1)
    plot(t,w(:,1),'k-',t,y,'b--')
    legend('Exact','Small Angle')
    title('Exact vs Approximate Solutions')
    xlabel('t')
    ylabel('\phi')
    subplot(2,1,2)
    plot(t,w(:,1)-y,'k-')
    title('Difference between Exact and Approximate')
    xlabel('t')
    ylabel( '\Delta\phi')
    
    f2=figure ('Name', 'energy');
    
    subplot(2,1,1)
    plot(t(ind(1):ind(3)),deltaE,'m:')
    title('\Delta')
    xlabel('t')
    ylabel( '\Delta')
    
     subplot(2,1,2)
    plot(t(ind(1):ind(3)),KE, 'c-', t(ind(1):ind(3)), U, 'm-')
    legend('Kinetic Energy','Potential Energy')
    xlabel('t')
    
    f3=figure ('Name', 'Velocity and Position');
    
    subplot(2,1,1)
    plot(t, w(:,2),'c-')
    title('$\dot{\theta (t)}, \theta _0=0.1$ ', 'Interpreter','latex')
    xlabel('t')
    ylabel('$\dot{\theta}$', 'Interpreter','latex')
    
    subplot(2,1,2)
    plot(t, w(:,1),'m-')
    title('\theta (t) ')
    xlabel('t')
    ylabel( '\theta')
    
    f4 =figure('Name', 'Phase Space');
    
    plot(w(:,1),w(:,2));
    title('Phase Space ')
    ylabel('$\dot{\theta}$', 'Interpreter','latex')
    xlabel('\theta')
end
end
%-------------------------------------------
%
function rdot = proj(t,r,g,R)
    rdot = [r(2); -g/R*sin(r(1))];
end

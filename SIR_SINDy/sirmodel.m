clc; clear;
load enso

% lambda is the sparsity knob that we can turn up or down to strive for
% parsimony in the "discovered" equations
% lambdavec = 0.0000001:0.00999:0.05;
% lambdavec = 0.002;
lambdavec = linspace(0.0000001,0.3,5);
for ll = 1:length(lambdavec)
%     lambda = 0.002;      % lambda is our sparsification knob.
    lambda = lambdavec(ll);
%     noiseKnob =  .5;

%     transmissibility = 1/30; %(i.e., probability of infection given contact between a susceptible and infected individual), 
%     contact_rate_between_SI = 1/40; 
%     duration_infectiousness = 7;
%     Ninfected = 1;
    
%     b = transmissibility.*contact_rate_between_SI;
%     g = 1/duration_infectiousness;
    b = 1.5;
    g = 1/4;
    Np = 500;
    timestep = 0.1;
    end_time = 40;
    
%     fprintf('typical time between contacts:')
%     disp(1/b)
    
%     fprintf('typical time until removal:')
%     disp(1/g)
    
    fprintf('R0:')
    disp(b/g)
    
%     fprintf('people in population:')
%     disp(Np)
    
%     fprintf('initial N infected:')
%     disp(Ninfected)
    
    n = 2;
    pars = [b g];
    
    rhsode = @(t,x,pars) [-pars(1).*x(1).*x(2); pars(1).*x(1).*x(2)-pars(2).*x(2); pars(2).*x(2)];
%     x0 = [Np-Ninfected; Ninfected; 0];           % ICs
    x0 = [1; 9*10^(-6); 0];           % ICs
    
    tspan = 0:timestep:end_time;             
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,x] = ode15s(@(t,x)rhsode(t,x,pars),tspan,x0,options);
    %add in the noise
%     x = abs(x + noiseKnob*randn(size(x,1),3));
    %smooth that noise
    f1 = fit(t,x(:,1),'smoothingspline','SmoothingParam',0.9);
    f2 = fit(t,x(:,2),'smoothingspline','SmoothingParam',0.9);
    f3 = fit(t,x(:,3),'smoothingspline','SmoothingParam',0.9);
    
    % hold on
    % plot(t,x(:,1), '-g', 'LineWidth',1.7)
    % plot(t,f1(t), '-r', 'LineWidth',1.7)
    % legend('noisy', 'smooth', 'Location', 'East')
    % hold off
    
%     x(:,1) = f1(t);
%     x(:,2) = f2(t);
%     x(:,3) = f3(t);

%     x = x./Np;
    
    % Using the equation to get the derivative
%     for i=1:length(t)
%         dx(i,:) = rhsode(t(i),x(i,:)',pars);
%     end
    
    % Calculating the derivative from the raw data
    for i=2:(length(t))
        dx(i-1,:) =  (x(i,:)-x(i-1,:))./(t(i,:)-t(i-1,:));
    end
    x = x(2:end,:);
    t = t(2:end,:);
    
    %% Build library and compute sparse regression
    polyorder = 3;
    Theta = poolData(x,n,polyorder);  % up to third order polynomials
    Xi = sparsifyDynamics(Theta,dx,lambda,n);
    [terms,termspars] = poolDataLIST({'x','y','z'},Xi,n,polyorder);
    
    termsparsSINDY = cell2mat(termspars(2:end,2:end));
    rhsodeSINDY = @(t,x,pars) [
        termsparsSINDY(1,1)*1 + termsparsSINDY(2,1)*x(1) + termsparsSINDY(3,1)*x(2) + termsparsSINDY(4,1)*x(1)*x(1) + termsparsSINDY(5,1)*x(1)*x(2) + termsparsSINDY(6,1)*x(2)*x(2) + termsparsSINDY(7,1)*x(1)*x(1)*x(1) + termsparsSINDY(8,1)*x(1)*x(1)*x(2) + termsparsSINDY(9,1)*x(1)*x(2)*x(2) + termsparsSINDY(10,1)*x(2)*x(2)*x(2);
        termsparsSINDY(1,2)*1 + termsparsSINDY(2,2)*x(1) + termsparsSINDY(3,2)*x(2) + termsparsSINDY(4,2)*x(1)*x(1) + termsparsSINDY(5,2)*x(1)*x(2) + termsparsSINDY(6,2)*x(2)*x(2) + termsparsSINDY(7,2)*x(1)*x(1)*x(1) + termsparsSINDY(8,2)*x(1)*x(1)*x(2) + termsparsSINDY(9,2)*x(1)*x(2)*x(2) + termsparsSINDY(10,2)*x(2)*x(2)*x(2);
        termsparsSINDY(1,3)*1 + termsparsSINDY(2,3)*x(1) + termsparsSINDY(3,3)*x(2) + termsparsSINDY(4,3)*x(1)*x(1) + termsparsSINDY(5,3)*x(1)*x(2) + termsparsSINDY(6,3)*x(2)*x(2) + termsparsSINDY(7,3)*x(1)*x(1)*x(1) + termsparsSINDY(8,3)*x(1)*x(1)*x(2) + termsparsSINDY(9,3)*x(1)*x(2)*x(2) + termsparsSINDY(10,3)*x(2)*x(2)*x(2);
    ];
    
    [tSINDY,xSINDY] = ode15s(@(tSINDY,xSINDY)rhsodeSINDY(tSINDY,xSINDY,pars),tspan,x0,options);
    
    figure
    hold on
    plot(t,x(:,1), '-g', 'LineWidth',1.7)
    plot(t,x(:,2), '-r', 'LineWidth',1.7)
    plot(t,x(:,3), '-c', 'LineWidth',1.7)
    plot(tSINDY, xSINDY(:,1), ':', 'color', '#057B10', 'LineWidth',2)
    plot(tSINDY, xSINDY(:,2), ':', 'color', '#9D2208', 'LineWidth',2)
    % plot(tSINDY, xSINDY(:,3), ':c', 'LineWidth',2)
    plot(tSINDY, 1 - xSINDY(:,1) - xSINDY(:,2), ':', 'color', '#0076a8', 'LineWidth',2)
    legend('S', 'I', 'R', 'S_{SINDy}', 'I_{SINDy}', 'R_{SINDy}', 'Location', 'East')
    title(sprintf('Raw Data (solid) vs SINDY Model (dotted); lambda = %f', lambda))
    xlim([0 end_time])
    ylim([0 1])
    hold off
end
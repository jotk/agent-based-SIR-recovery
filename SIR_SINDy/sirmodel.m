clc; clear;
load enso
calculate_RMSE=@(a,b) sqrt(mean((a(:)-b(:)).^2));

parms = readtable("D:/Google Drive/1Chris/CSCI5423/beta_gamma_output.csv",'Delimiter', ',', "NumHeaderLines", 0);
st = (find(table2array(parms(:,1))==1));
parmsdf = table2array(parms(st:size(parms,1),[1,2,6]));
parmsdf = horzcat(parmsdf, (parmsdf(:,2))./(parmsdf(:,3)));
% mean(parmsdf(1:32,4))


M = readtable("D:/Google Drive/1Chris/CSCI5423/sir_output.csv",'Delimiter', ',', "NumHeaderLines", 0);
st = (find(table2array(M(:,1))==1)-1);
rr = st:size(M,1);
cc = [1,2,6,10];
%columns are: time, I, R, S
df = table2array(M(rr,cc));
%columns are: time, S, I, R
df = df(:,[1,4,2,3]);
%check that the total N is the same for each timestep
sum(df(:,2:4), 2);
end_time = df(size(df,1),1);
Np = sum(df(1,2:4));
% beta = table2array(M(5,2))./100;
% gamma = (table2array(M(5,1))./100)*(str2double(table2array(M(5,4)))./100);

% areaInfect = (3/61)*(3/61);
% peoplePerArea = (Np - 1)/(61*61);
% contactpertime = Np*areaInfect*peoplePerArea %people divided by number of patches times 9 - (individual of interest)
% R0 = (table2array(M(5,2))./100)*(contactpertime)*(table2array(M(5,1))) %(infection/contact)*(contact/time)*(time/infection)
R0 = mean(parmsdf(1:10,4))
beta = mean(parmsdf(1:10,2))
gamma = mean(parmsdf(1:10,3))

bestError = 1000000;
% tlen = 200;
tlenvec = linspace(50, 200, 3);
for tleni = 1:length(tlenvec)
    tlen = tlenvec(tleni);
    smoothingParametervec = linspace(.2, .99, 3);
    % smoothingParameter = 0.70;
    for smoothingParameteri = 1:length(smoothingParametervec)
        smoothingParameter = smoothingParametervec(smoothingParameteri);
    
        % lambda is the sparsity knob that we can turn up or down to strive for
        % parsimony in the "discovered" equations
        % lambdavec = 0.0000001:0.00999:0.05;
        % lambdavec = 0.01;
        lambdavec = linspace(0.0005, .1, 20);
        for ll = 1:length(lambdavec)
        %     lambda = 0.002;      % lambda is our sparsification knob.
            lambda = lambdavec(ll);
            lambda
        %     noiseKnob =  .5;
        
        %     transmissibility = 1/30; %(i.e., probability of infection given contact between a susceptible and infected individual), 
        %     contact_rate_between_SI = 1/40; 
        %     duration_infectiousness = 7;
        %     Ninfected = 1;
            
        %     b = transmissibility.*contact_rate_between_SI;
        %     g = 1/duration_infectiousness;
        %     b = 1.5;
        %     g = 1/4;
        %     Np = 500;
        %     timestep = 0.1;
        %     end_time = 40;
            
        %     fprintf('typical time between contacts:')
        %     disp(1/b)
            
        %     fprintf('typical time until removal:')
        %     disp(1/g)
            
        %     fprintf('R0:')
        %     disp(b/g)
            
        %     fprintf('people in population:')
        %     disp(Np)
            
        %     fprintf('initial N infected:')
        %     disp(Ninfected)
            
            n = 2;
        %     pars = [beta gamma];
            
        %     rhsode = @(t,x,pars) [-pars(1).*x(1).*x(2); pars(1).*x(1).*x(2)-pars(2).*x(2); pars(2).*x(2)];
        %     x0 = [Np-Ninfected; Ninfected; 0];           % ICs
        %     x0 = [1; 9*10^(-6); 0];           % ICs
        %     x0 = [df(1,2); df(1,3); df(1,4)];           % ICs
            
        %     tspan = 0:timestep:end_time;             
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        %     [t,x] = ode15s(@(t,x)rhsode(t,x,pars),tspan,x0,options);
            % add in the noise
            % x = abs(x + noiseKnob*randn(size(x,1),3));
        %     tspan = df(:,1);
            tspan = linspace(0,df(end,1),tlen)';
            t = tspan;
        %     x = df(:,2:4);
        
            %smooth that noise
            f1 = fit(df(:,1),df(:,2),'smoothingspline','SmoothingParam',smoothingParameter);
            f2 = fit(df(:,1),df(:,3),'smoothingspline','SmoothingParam',smoothingParameter);
            f3 = fit(df(:,1),df(:,4),'smoothingspline','SmoothingParam',smoothingParameter);
            
        %     hold on
        %     plot(t,x(:,1), '-g', 'LineWidth',1.7)
        %     plot(t,f1(t), '-r', 'LineWidth',1.7)
        %     legend('noisy', 'smooth', 'Location', 'East')
        %     hold off
            
        %     x(:,1) = f1(t);
        %     x(:,2) = f2(t);
        %     x(:,3) = f3(t);
        
            x = horzcat(f1(t), f2(t), f3(t));
        
            x0 = [x(1,1); x(1,2); x(1,3)];           % ICs
        
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
            polyorder = 2;
            Theta = poolData(x,n,polyorder);  % up to third order polynomials
            Xi = sparsifyDynamics(Theta,dx,lambda,n);
            [terms,termspars] = poolDataLIST({'x','y','z'},Xi,n,polyorder);
            
            termsparsSINDY = cell2mat(termspars(2:end,2:end));
            % rhsodeSINDY = @(t,x) [
            %     termsparsSINDY(1,1)*1 + termsparsSINDY(2,1)*x(1) + termsparsSINDY(3,1)*x(2) + termsparsSINDY(4,1)*x(1)*x(1) + termsparsSINDY(5,1)*x(1)*x(2) + termsparsSINDY(6,1)*x(2)*x(2) + termsparsSINDY(7,1)*x(1)*x(1)*x(1) + termsparsSINDY(8,1)*x(1)*x(1)*x(2) + termsparsSINDY(9,1)*x(1)*x(2)*x(2) + termsparsSINDY(10,1)*x(2)*x(2)*x(2);
            %     termsparsSINDY(1,2)*1 + termsparsSINDY(2,2)*x(1) + termsparsSINDY(3,2)*x(2) + termsparsSINDY(4,2)*x(1)*x(1) + termsparsSINDY(5,2)*x(1)*x(2) + termsparsSINDY(6,2)*x(2)*x(2) + termsparsSINDY(7,2)*x(1)*x(1)*x(1) + termsparsSINDY(8,2)*x(1)*x(1)*x(2) + termsparsSINDY(9,2)*x(1)*x(2)*x(2) + termsparsSINDY(10,2)*x(2)*x(2)*x(2);
            %     termsparsSINDY(1,3)*1 + termsparsSINDY(2,3)*x(1) + termsparsSINDY(3,3)*x(2) + termsparsSINDY(4,3)*x(1)*x(1) + termsparsSINDY(5,3)*x(1)*x(2) + termsparsSINDY(6,3)*x(2)*x(2) + termsparsSINDY(7,3)*x(1)*x(1)*x(1) + termsparsSINDY(8,3)*x(1)*x(1)*x(2) + termsparsSINDY(9,3)*x(1)*x(2)*x(2) + termsparsSINDY(10,3)*x(2)*x(2)*x(2);
            % ];
            rhsodeSINDY = @(t,x) [
                termsparsSINDY(1,1)*1 + termsparsSINDY(2,1)*x(1) + termsparsSINDY(3,1)*x(2) + termsparsSINDY(4,1)*x(1)*x(1) + termsparsSINDY(5,1)*x(1)*x(2) + termsparsSINDY(6,1)*x(2)*x(2);
                termsparsSINDY(1,2)*1 + termsparsSINDY(2,2)*x(1) + termsparsSINDY(3,2)*x(2) + termsparsSINDY(4,2)*x(1)*x(1) + termsparsSINDY(5,2)*x(1)*x(2) + termsparsSINDY(6,2)*x(2)*x(2);
                termsparsSINDY(1,3)*1 + termsparsSINDY(2,3)*x(1) + termsparsSINDY(3,3)*x(2) + termsparsSINDY(4,3)*x(1)*x(1) + termsparsSINDY(5,3)*x(1)*x(2) + termsparsSINDY(6,3)*x(2)*x(2);
            ];
            
        %     termsparsTruth = [0 0 0; 
        %         0 0 0; 
        %         0 -gamma gamma; 
        %         0 0 0; 
        %         -beta beta 0; 
        %         0 0 0; 
        %         0 0 0; 
        %         0 0 0; 
        %         0 0 0; 
        %         0 0 0];
            termsparsTruth = [0 0 0; 
                0 0 0; 
                0 -gamma gamma; 
                0 0 0; 
                -beta beta 0; 
                0 0 0;];
        
            err = abs(termsparsSINDY(:,1:2) - termsparsTruth(:,1:2));
            err(3,2) = (err(3,2)+2)^3;
            err(5,1) = (err(5,1)+2)^3;
            err(5,2) = (err(5,2)+2)^3;
        %     absError = sum(sum(err));
        
            [tSINDY,xSINDY] = ode15s(@(tSINDY,xSINDY)rhsodeSINDY(tSINDY,xSINDY),t,x0,options);
        
            if length(x(:,2)) == length(xSINDY(:,2))
            %     length(x(:,2))
            %     length(xSINDY(:,2))
            
                absError = calculate_RMSE(termsparsSINDY(:,1),termsparsTruth(:,1)) + calculate_RMSE(termsparsSINDY(:,2),termsparsTruth(:,2)) + calculate_RMSE(x(:,2),xSINDY(:,2))/20 + err(3,2) + err(5,1) + err(5,2);
            %     absError = calculate_RMSE(termsparsSINDY(:,1),termsparsTruth(:,1)) + calculate_RMSE(termsparsSINDY(:,2),termsparsTruth(:,2)) + err(3,2) + err(5,1) + err(5,2);
%                 absError
                %     absError = sum([err(3,2),  err(5,1),  err(5,2)]);
%                 absError = calculate_RMSE(x(:,1),xSINDY(:,1))/2 + calculate_RMSE(x(:,2),xSINDY(:,2))/2 + err(3,2) + err(5,1) + err(5,2);
            
            
                if absError < bestError
                    bestLambda = lambda;
                    besttlen = tlen;
                    bestsmoothingParameter = smoothingParameter;
                    bestError = absError;
                end
            end
        %     figure
        %     hold on
        %     plot(t,x(:,1), '-g', 'LineWidth',1.7)
        %     plot(t,x(:,2), '-r', 'LineWidth',1.7)
        %     plot(t,x(:,3), '-c', 'LineWidth',1.7)
        %     plot(tSINDY, xSINDY(:,1), ':', 'color', '#057B10', 'LineWidth',2)
        %     plot(tSINDY, xSINDY(:,2), ':', 'color', '#9D2208', 'LineWidth',2)
        % %     plot(tSINDY, xSINDY(:,3), ':c', 'LineWidth',2)
        %     plot(tSINDY, Np - xSINDY(:,2) - xSINDY(:,1), ':c', 'LineWidth',2)
        %     plot(tSINDY, 1 - xSINDY(:,1) - xSINDY(:,2), ':', 'color', '#0076a8', 'LineWidth',2)
        %     legend('S', 'I', 'R', 'S_{SINDy}', 'I_{SINDy}', 'R_{SINDy}', 'Location', 'East')
        %     title(sprintf('Raw Data (solid) vs SINDY Model (dotted); lambda = %f', lambda))
        %     xlim([0 end_time])
        %     ylim([0 Np])
        %     hold off
        end


    end   
end    
 













%% best choices



lambda = bestLambda;
% lambda = 0.0005;
tlen = besttlen;
smoothingParameter = bestsmoothingParameter;
%     noiseKnob =  .5;

%     transmissibility = 1/30; %(i.e., probability of infection given contact between a susceptible and infected individual), 
%     contact_rate_between_SI = 1/40; 
%     duration_infectiousness = 7;
%     Ninfected = 1;

%     b = transmissibility.*contact_rate_between_SI;
%     g = 1/duration_infectiousness;
%     b = 1.5;
%     g = 1/4;
%     Np = 500;
%     timestep = 0.1;
%     end_time = 40;

%     fprintf('typical time between contacts:')
%     disp(1/b)

%     fprintf('typical time until removal:')
%     disp(1/g)

%     fprintf('R0:')
%     disp(b/g)

%     fprintf('people in population:')
%     disp(Np)

%     fprintf('initial N infected:')
%     disp(Ninfected)

n = 2;
%     pars = [beta gamma];

%     rhsode = @(t,x,pars) [-pars(1).*x(1).*x(2); pars(1).*x(1).*x(2)-pars(2).*x(2); pars(2).*x(2)];
%     x0 = [Np-Ninfected; Ninfected; 0];           % ICs
%     x0 = [1; 9*10^(-6); 0];           % ICs
%     x0 = [df(1,2); df(1,3); df(1,4)];           % ICs

%     tspan = 0:timestep:end_time;             
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
%     [t,x] = ode15s(@(t,x)rhsode(t,x,pars),tspan,x0,options);
% add in the noise
% x = abs(x + noiseKnob*randn(size(x,1),3));
%     tspan = df(:,1);
tspan = linspace(0,df(end,1),tlen)';
t = tspan;
%     x = df(:,2:4);

%smooth that noise
f1 = fit(df(:,1),df(:,2),'smoothingspline','SmoothingParam',smoothingParameter);
f2 = fit(df(:,1),df(:,3),'smoothingspline','SmoothingParam',smoothingParameter);
f3 = fit(df(:,1),df(:,4),'smoothingspline','SmoothingParam',smoothingParameter);

%     hold on
%     plot(t,x(:,1), '-g', 'LineWidth',1.7)
%     plot(t,f1(t), '-r', 'LineWidth',1.7)
%     legend('noisy', 'smooth', 'Location', 'East')
%     hold off

%     x(:,1) = f1(t);
%     x(:,2) = f2(t);
%     x(:,3) = f3(t);

x = horzcat(f1(t), f2(t), f3(t));

x0 = [x(1,1); x(1,2); x(1,3)];           % ICs

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
polyorder = 2;
Theta = poolData(x,n,polyorder);  % up to third order polynomials
Xi = sparsifyDynamics(Theta,dx,lambda,n);
[terms,termspars] = poolDataLIST({'x','y','z'},Xi,n,polyorder);

termsparsSINDY = cell2mat(termspars(2:end,2:end));
% rhsodeSINDY = @(t,x) [
%     termsparsSINDY(1,1)*1 + termsparsSINDY(2,1)*x(1) + termsparsSINDY(3,1)*x(2) + termsparsSINDY(4,1)*x(1)*x(1) + termsparsSINDY(5,1)*x(1)*x(2) + termsparsSINDY(6,1)*x(2)*x(2) + termsparsSINDY(7,1)*x(1)*x(1)*x(1) + termsparsSINDY(8,1)*x(1)*x(1)*x(2) + termsparsSINDY(9,1)*x(1)*x(2)*x(2) + termsparsSINDY(10,1)*x(2)*x(2)*x(2);
%     termsparsSINDY(1,2)*1 + termsparsSINDY(2,2)*x(1) + termsparsSINDY(3,2)*x(2) + termsparsSINDY(4,2)*x(1)*x(1) + termsparsSINDY(5,2)*x(1)*x(2) + termsparsSINDY(6,2)*x(2)*x(2) + termsparsSINDY(7,2)*x(1)*x(1)*x(1) + termsparsSINDY(8,2)*x(1)*x(1)*x(2) + termsparsSINDY(9,2)*x(1)*x(2)*x(2) + termsparsSINDY(10,2)*x(2)*x(2)*x(2);
%     termsparsSINDY(1,3)*1 + termsparsSINDY(2,3)*x(1) + termsparsSINDY(3,3)*x(2) + termsparsSINDY(4,3)*x(1)*x(1) + termsparsSINDY(5,3)*x(1)*x(2) + termsparsSINDY(6,3)*x(2)*x(2) + termsparsSINDY(7,3)*x(1)*x(1)*x(1) + termsparsSINDY(8,3)*x(1)*x(1)*x(2) + termsparsSINDY(9,3)*x(1)*x(2)*x(2) + termsparsSINDY(10,3)*x(2)*x(2)*x(2);
% ];
rhsodeSINDY = @(t,x) [
    termsparsSINDY(1,1)*1 + termsparsSINDY(2,1)*x(1) + termsparsSINDY(3,1)*x(2) + termsparsSINDY(4,1)*x(1)*x(1) + termsparsSINDY(5,1)*x(1)*x(2) + termsparsSINDY(6,1)*x(2)*x(2);
    termsparsSINDY(1,2)*1 + termsparsSINDY(2,2)*x(1) + termsparsSINDY(3,2)*x(2) + termsparsSINDY(4,2)*x(1)*x(1) + termsparsSINDY(5,2)*x(1)*x(2) + termsparsSINDY(6,2)*x(2)*x(2);
    termsparsSINDY(1,3)*1 + termsparsSINDY(2,3)*x(1) + termsparsSINDY(3,3)*x(2) + termsparsSINDY(4,3)*x(1)*x(1) + termsparsSINDY(5,3)*x(1)*x(2) + termsparsSINDY(6,3)*x(2)*x(2);
];

%     termsparsTruth = [0 0 0; 
%         0 0 0; 
%         0 -gamma gamma; 
%         0 0 0; 
%         -beta beta 0; 
%         0 0 0; 
%         0 0 0; 
%         0 0 0; 
%         0 0 0; 
%         0 0 0];
termsparsTruth = [0 0 0; 
    0 0 0; 
    0 -gamma gamma; 
    0 0 0; 
    -beta beta 0; 
    0 0 0;];

err = abs(termsparsSINDY(:,1:2) - termsparsTruth(:,1:2));
err(3,2) = (err(3,2)+1)^3;
err(5,1) = (err(5,1)+1)^3;
err(5,2) = (err(5,2)+1)^3;
%     absError = sum(sum(err));

[tSINDY,xSINDY] = ode15s(@(tSINDY,xSINDY)rhsodeSINDY(tSINDY,xSINDY),t,x0,options);

if length(x(:,2)) == length(xSINDY(:,2))
%     length(x(:,2))
%     length(xSINDY(:,2))

    absError = calculate_RMSE(termsparsSINDY(:,1),termsparsTruth(:,1)) + calculate_RMSE(termsparsSINDY(:,2),termsparsTruth(:,2)) + calculate_RMSE(x(:,2),xSINDY(:,2))/20 + err(3,2) + err(5,1) + err(5,2);
%     absError = calculate_RMSE(termsparsSINDY(:,1),termsparsTruth(:,1)) + calculate_RMSE(termsparsSINDY(:,2),termsparsTruth(:,2)) + err(3,2) + err(5,1) + err(5,2);
%     absError
    %     absError = sum([err(3,2),  err(5,1),  err(5,2)]);
    absError = calculate_RMSE(x(:,1),xSINDY(:,1))/2 + calculate_RMSE(x(:,2),xSINDY(:,2))/2 + err(3,2) + err(5,1) + err(5,2);


    if absError < bestError
        bestLambda = lambda;
        besttlen = tlen;
        bestsmoothingParameter = smoothingParameter;
        bestError = absError;
    end
end
figure
hold on
plot(t,x(:,1), '-g', 'LineWidth',1.7)
plot(t,x(:,2), '-r', 'LineWidth',1.7)
plot(t,x(:,3), '-c', 'LineWidth',1.7)
plot(tSINDY, xSINDY(:,1), ':', 'color', '#057B10', 'LineWidth',2)
plot(tSINDY, xSINDY(:,2), ':', 'color', '#9D2208', 'LineWidth',2)
%             plot(tSINDY, xSINDY(:,3), ':c', 'LineWidth',2)
plot(tSINDY, Np - xSINDY(:,2) - xSINDY(:,1), ':', 'color', '#0076a8', 'LineWidth',2)
%             plot(tSINDY, 1 - xSINDY(:,1) - xSINDY(:,2), ':', 'color', '#0076a8', 'LineWidth',2)
legend('S', 'I', 'R', 'S_{SINDy}', 'I_{SINDy}', 'R_{SINDy}', 'Location', 'East')
title(sprintf('Raw Data (solid) vs SINDY Model (dotted); lambda = %f', lambda))
xlim([0 end_time])
ylim([0 Np])
hold off

bestLambda
besttlen
bestsmoothingParameter
horzcat(termsparsTruth(:,1), termsparsSINDY(:,1))
horzcat(termsparsTruth(:,2), termsparsSINDY(:,2))
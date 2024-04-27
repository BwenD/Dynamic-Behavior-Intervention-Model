clear;
clc;
% Input parameters
% s(1) to s(4) are the time lags for SN, WSA, SE, and PBC respectively [1,14];
% s(5) to s(9) are the reciprocals of the inventory time coefficients for SN, WSA, SE, PBC, and WSI respectively [0,1];
% s(10) to s(13) are the inputs for SN, WSA, SE, and PBC in the first 0-14 days [6.3,7], [6.64,7], [5.98,7], [6.19,7];
s=[1.0 11.0 11.0 11.0 0.8737936771949261 0.9488618287169458 0.09173100567391806 0.09149206684462165 0.08401530555962253 6.759295751133038 6.576326656163232 6.472543905871165 6.593459695798504]; % Result: Si=-0.8893632278085368, generation=326;
data=[0 0.002106873 0.002632234 0.002225361 0.003013402 0.001904467 0.000457595]; % Actual measured water use B3_Fb
input=[6.14 6.432 5.974 6.066 5.599238]; % y5(0) is calculated from the model
mubiao_x=[0 7 14 21 28 35 42]; % Target x values corresponding to specific time points
mubiao_y=data; % Target y values which are the actual water use measurements
lags=[s(1) s(2) s(3) s(4)]; % Input various time lags
t_range=[0 42]; % Time range for solving the differential equations

% Solve the differential equation
f=@(t,y,Z)[s(5).*(B3_Fb_a1(t,s(10))-y(1));
    s(6).*(B3_Fb_a2(t,s(11))-y(2));
    s(7).*(B3_Fb_a3(t,s(12))-y(3));
    s(8).*(B3_Fb_a4(t,s(13))-y(4));
    s(9).*(0.195.*Z(1,1)+0.24.*Z(2,1)+0.251.*Z(3,1)+0.224.*Z(4,1)-y(5))];
sol=dde23(f,lags,input,t_range); % Using MATLAB's DDE solver

% Plot the solution
plot(sol.x,sol.y,'-o'); % Plotting the solution
xlabel('Time t'); % Label for the x-axis representing time
ylabel('Solution y'); % Label for the y-axis representing the solution values
legend('SN','WSA','SE','PBC','WSI','Location','NorthWest'); % Adding a legend to the plot
output_x=sol.x; % Storing the time values for potential further use
output_y=sol.y; % Storing the solution values for potential further use

% Calculate the objective function
scx=sol.x; % Corresponding x values from the DDEs solution, which is the time axis
scy=sol.y(5,:); % Corresponding y(5) values from the DDEs solution, representing the value of WSI inventory
scy_w=[];
for i=1:length(mubiao_x)
    [~, I(1,i)]=min(abs(scx(:)-mubiao_x(1,i))); % Find the closest time step to the actual time step
    scy_w(1,i)=scy(1,I(1,i)); % The value of y(5) corresponding to the closest time step
end

Mu=corrcoef(mubiao_y,scy_w); % Objective function as the negative of the correlation coefficient (Pearson)
W=-Mu(1,2); % The objective function value
% W=sqrt(sum(W1)); % Alternative objective function (Euclidean distance), square root of the sum of squared differences
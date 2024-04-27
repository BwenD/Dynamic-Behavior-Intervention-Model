function W=B3_Inf_Mubiaohanshu(s)
% s(1) to s(4) are the time lags for SN, WSA, SE, and PBC respectively, within the range [1,14];
% s(5) to s(9) are the reciprocals of the inventory time coefficients for SN, WSA, SE, PBC, and WSI respectively, within the range [0,1];
% s(10) to s(13) are the inputs for SN, WSA, SE, and PBC in the first 0-14 days, with ranges [6.3,7], [6.64,7], [5.98,7], [6.19,7];
data=[0 0.00122268 0.003172371 0.002749691 0.002374227 0.002184742 0.001062062]; % Actual measured water use B3_Inf
mubiao_x=[0 7 14 21 28 35 42]; % Target x values corresponding to specific time points
mubiao_y=data; % Target y values which are the actual water use measurements
lags=[s(1) s(2) s(3) s(4)]; % Input various time lags
input=[5.58 6.566 5.874 6.102 5.505162]; % y5(0) is calculated from the model
t_range=[0 42]; % Time range for solving the differential equations

% Solve the differential equation
f=@(t,y,Z)[s(5).*(B3_Inf_a1(t,s(10))-y(1));
    s(6).*(B3_Inf_a2(t,s(11))-y(2));
    s(7).*(B3_Inf_a3(t,s(12))-y(3));
    s(8).*(B3_Inf_a4(t,s(13))-y(4));
    s(9).*(0.195.*Z(1,1)+0.24.*Z(2,1)+0.251.*Z(3,1)+0.224.*Z(4,1)-y(5))];
sol=dde23(f,lags,input,t_range); % Using MATLAB's DDE solver

% Plot the solution (commented out)
%plot(sol.x,sol.y,'-o'); % Plotting the solution
%xlabel('Time t');
%ylabel('Solution y');
%legend('SN','WSA','SE','PBC','WSI','Location','NorthWest');
%xlswrite('x.xlsx',sol.x); % Writing x values to an Excel file
%xlswrite('y.xlsx',sol.y); % Writing y values to an Excel file

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
end
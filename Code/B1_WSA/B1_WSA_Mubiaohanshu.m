function W=B1_WSA_Mubiaohanshu(s)
% s(1) = Time lag from WSA to WSI [1,14];
% s(2) = The reciprocal of the inventory time coefficient for WSA [0,1];
% s(3) = The reciprocal of the inventory time coefficient for WSI [0,1];
% s(4) = The input for WSA in the first 0-14 days [4.95,7];
data=[0 0.011 0.016 0.043 0.004]; % Input actual measured water use WSB for Office Building B1_WSA
mubiao_x=[0 7 14 21 28]; % Target x values corresponding to specific time points
mubiao_y=data; % Target y values which are the actual water use measurements
lags=[1 s(1) 1 1]; % Input various time lags
input=[4.3 4.95 4.4 5 4.27255];
t_range=[0 30]; % Time range for solving the differential equations

% Solve the differential equation
f=@(t,y,Z)[(B1_WSA_a1(t)-y(1));
    s(2).*(B1_WSA_a2(t,s(4))-y(2));
    (B1_WSA_a3(t)-y(3));
    (B1_WSA_a4(t)-y(4));
    s(3).*(0.227.*Z(1,1)+0.215.*Z(2,1)+0.238.*Z(3,1)+0.237.*Z(4,1)-y(5))];
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

Mu=corrcoef(mubiao_y,scy_w); % Objective function as the negative of the correlation coefficient
W=-Mu(1,2); % The objective function value
% W=sqrt(sum(W1)); % Alternative objective function (Euclidean distance), square root of the sum of squared differences
end
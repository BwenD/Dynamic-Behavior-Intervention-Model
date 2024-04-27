function W=B2_SN_Mubiaohanshu(s)
% s(1) = Time lag from SN to WSI [1,14];
% s(2) = The reciprocal of the inventory time coefficient for SN [0,1];
% s(3) = The reciprocal of the inventory time coefficient for WSI [0,1];
% s(4) = The input for SN in the first 0-14 days [4.17,7];
data=[0 0.047 0.044 0.07 0.029]; % Input actual measured water use WSB for Office Building B2_SN
mubiao_x=[0 7 14 21 28]; % Target x values corresponding to specific time points
mubiao_y=data; % Target y values which are the actual water use measurements
lags=[s(1) 1 1 1]; % Input various time lags
input=[4.17 4.89 4.65 4.57 4.3981]; % Initial conditions for the differential equation
t_range=[0 30]; % Time range for solving the differential equations

% Solve the differential equation
f=@(t,y,Z)[s(2).*(B2_SN_a1(t,s(4))-y(1));
    (B2_SN_a2(t)-y(2));
    (B2_SN_a3(t)-y(3));
    (B2_SN_a4(t)-y(4));
    s(3).*(0.214.*Z(1,1)+0.231.*Z(2,1)+0.225.*Z(3,1)+0.291.*Z(4,1)-y(5))];
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
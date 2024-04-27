function W=B1_SN_Mubiaohanshu(s)
% s(1) = Time lag from SN to WSI [1,14]; 
% s(2) = The reciprocal of the inventory time coefficient for SN [0,1];
% s(3) = The reciprocal of the inventory time coefficient for WSI [0,1]; 
% s(4) = The input for SN in the first 0-14 days [4.3,7];
data=[0.000 0.050 0.080 0.141 0.079]; % Input actual measured water use WSB, Office Building B1_SN
mubiao_x=[0 7 14 21 28]; % Target x values
mubiao_y=data; % Target y values
lags=[s(1) 1 1 1]; % Input various time lags
input=[4.3 4.95 4.4 5 4.27255];
t_range=[0 30]; % Time range for solving the differential equations

% Solve the differential equation
f=@(t,y,Z)[s(2).*(B1_SN_a1(t,s(4))-y(1)); 
    (B1_SN_a2(t)-y(2)); 
    (B1_SN_a3(t)-y(3)); 
    (B1_SN_a4(t)-y(4)); 
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
scx=sol.x; % Corresponding x values from the differential equation, which is the time axis
scy=sol.y(5,:); % Corresponding y(5) values from the differential equation, the value of WSI inventory
scy_w=[];
for i=1:length(mubiao_x)
    [~, I(1,i)]=min(abs(scx(:)-mubiao_x(1,i))); % Find the time step that is closest to the actual time step
    scy_w(1,i)=scy(1,I(1,i)); % The value of y(5) corresponding to the closest time step
end
Mu=corrcoef(mubiao_y,scy_w); % Correlation coefficient between target y values and calculated y(5) values
W=-Mu(1,2); % Objective function, negative of the correlation coefficient
% W=sqrt(sum(W1)); % Alternative objective function, square root of the sum of squared differences
end
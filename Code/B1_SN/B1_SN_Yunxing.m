clear;
clc;
A=[-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1];
b=[0;0;0;0]; % Set the optimization constraints, time lags must be greater than 0
optimtool % Call the optimization tool
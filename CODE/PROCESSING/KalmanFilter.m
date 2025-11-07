function [param, param_sigma, float] = KalmanFilter(omc, A, x_pred, Q_x_pred, Q_l)
% Calculate PPP solution with Kalman Filter following 
% [00]: p.30-32
% [01]: p.244-248 
% https://www.kalmanfilter.net/stateUpdate.html 
% https://www.kalmanfilter.net/covUpdate.html
%
% INPUT: 
%   omc             vector, observed minus computed
%   A               matrix, Design Matrix / Observation Matrix
%   x_pred          vector, predicted parameters
%   Q_x_pred        matrix, predicted covariance matrix of parameters
%   Q_l             matrix, covariance matrix of observations (from createObsCovariance.m)
% OUTPUT:
%   param           parameter vector after Kalman Filtering 
%   param_sigma     covariance matrix after Kalman Filtering
%   float           boolean, valid float solution?
%  
% Revision:
%   2025/02/18, MFWG: moved calculation of residuals outside this function
%   2025 spring, MFWG: strongly revised (e.g., new formula for covariance)
%   2025/10/14, MFWG: changed input/output from structs to variables
% 
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************


% Note: prediction (x_pred, Q_x_pred) was calculated in adjPrep_ZD/_iono/_DCM.m


%% Preparations

% Check for NaNs in observed-minus-computed, for example, missing
% observations (e.g., missing L5 observations for GPS)
idx_nan = isnan(omc);
% set rows of missing observations to zero and remove them from adjustment
omc(idx_nan) = 0;
A(idx_nan, :) = 0;

% Identity Matrix with the size of the Design Matrix / Observation Matrix
I = eye(size(A,2));


%% Filtering
% Gain computation (Kalman weight)
K = (Q_x_pred * A') / (Q_l + A*Q_x_pred*A'); 	% [01]: (7.118)

% State update / parameter update
x_ = x_pred + K*omc;           		% [01]: (7.119), omc should be right

% Covariance update (check https://www.kalmanfilter.net/covUpdate.html)
Q_x_ = (I - K*A)*Q_x_pred*(I - K*A)' + K*Q_l*K';
% Q_x_ = (I - K*A)*Q_x_pred;        % [01]: (7.120), but https://www.kalmanfilter.net/simpCovUpdate.html

% ensure that covariance matrix is symmetric (e.g., for LAMBDA method)
Q_x_ = (Q_x_ + Q_x_')/2;	



%% save results
param = x_;       	    % save state / parameter vector
param_sigma = Q_x_; 	% save covariance matrix of parameters
float = true;    	    % valid float solution


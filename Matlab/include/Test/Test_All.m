%% Test_All.m
% 
% Runs all unit tests in the "include/Test" folder
% 
% @author: Matt Marti
% @date: 2019-03-06

clear, clc, clear global
constants


%% Run tests

% Numerical
Test_skewsym
Test_trapazoidIntegration

% Attitude
Test_quat2dircos
Test_insSinc
% Test_increQuatWithAngles

% Geod
Test_ecef2latlon
Test_latlon2ecef

% InertialNavigation
Test_gravitymodel
Test_insMechanization


%% Out
fprintf('---All Passed---\n');
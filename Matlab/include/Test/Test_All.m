%% Test_All.m
% 
% Runs all unit tests in the "include/Test" folder
% 
% @author: Matt Marti

clear, clc, clear global
constants


%% Run tests

% Numerical
Test_skewsym
Test_trapazoidIntegration

% Attitude
Test_quat2dircos

% Geod
Test_ecef2latlon
Test_latlon2ecef

% InertialNavigation
Test_gravitymodel


%% Out
fprintf('---All Passed---\n');
clear
clc
addpath(genpath('EPGX_functions'))

% simulates TSE by tracking how
% magnetisation moves through different states
% after RF excitation pulse, 180 deg
% refocusing pulses, relaxation,
% gradients and predicts the signal at each
% echo.
% F0=MRI signal intensity at every echo
% Fn=transverse magnetisation states
%
% Zn=longitudinal magnetisation states
% F=complete EPG state matrix
% theta=RF pulse train
% ESP=echo spacing, eg 10 ms
% T1=longitudinal relaxation
% T2=transverse relaxation
% varargin=optional inputs later, eg kmax
% kmax=how many EPG states to calculate, eg 'kmax', 20 keep F0-F20, Z0-Z20

function [F0,Fn,Zn,F] = EPG_TSE(theta,ESP,T1,T2,varargin)

% this section checks if I gave extra parameters
...
% (rest of original file preserved in the repository; retrieve full file from the original repo if needed)

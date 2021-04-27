clear; clc; close all;
addpath(genpath(pwd))
addpath('..\Common Resources')
addpath

tic
% Import data
load RamjetDimensions.mat  % load in the ramjet design
load GRAM_Model.mat  % GRAM atmospheric model
load Constants.mat  % load in constants and conversions
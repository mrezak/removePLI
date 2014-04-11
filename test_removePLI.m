%
%   This sample script demonstrates the usage of "removePLI.m" 
%   and illustrates the performance of the power line interference
%   removal algorithm on intracortical, ECoG, EEG and ECG signals.
%
%   The algorithm is proposed in:
%   M. R. Keshtkaran and Z. Yang, "A fast, robust algorithm for power line 
%   interference cancellation in neural recording," J. Neural Eng., vol. 11,
%   no. 2, p. 026017, Apr. 2014.

%
%   This script requires the .m file "removePLI.m" and .mat files:
%   "intracortical.mat" (from Plexon), "ECoG.mat" (from Riken lab),
%   "EEG.mat" (from Physionet) and "ECG.mat" (from Physionet).
%   
%   Licence:
%   Downloaded from: https://github.com/mrezak/removePLI
%	Author: Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
%   Copyright (c) 2013, Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
%   All rights reserved.
%	This program is provided "AS IS" for non-commercial, educational 
%	and reseach purpose only. Any commercial use, of any kind, of 
%	this program is prohibited. The Copyright notice should remain intact.
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Intracortical
load 'intracortical.mat';

fs = 40000; %Sample rate (Hz)
x = data/2^15*10*1e3; % Scaling to mV

M = 6;  %number of harmonics to remove
B = [50 .2 1] ;   
P = [0.1 4 1] ;   
W = 2;

idx = fs:length(x); % portion of signal to be used for PSD estimation
figure; movegui(gcf,'northwest')
subplot(2,1,1)
pwelch(x(idx),[],[],[],fs);    %PSD of the original signal
title('Intracortical, PSD of signal BEFORE interference cancellation')
axis([0 .350 -15 55])

% Interference cancellation
s = removePLI(x, fs, M, B, P, W);
subplot(2,1,2); pwelch(s(idx),[],[],[],fs);    %PSD of the filtered signal
title('Intracortical, PSD of signal AFTER interference cancellation')
axis([0 .350 -15 55])

%% ECoG
load 'ECoG.mat';

fs = 1000; %Sample rate (Hz)
x = data; % Scaling to mV

M = 7;  %number of harmonics to remove
B = [50 .2 1] ;   
P = [0.1 4 1] ;   
W = 1;

idx = 1*fs:length(x); % portion of signal to be used for PSD estimation

figure; movegui(gcf,'north')
subplot(2,1,1)
pwelch(x(idx),[],[],[],fs);    %PSD of the original signal
title('ECoG, PSD of signal BEFORE interference cancellation')
axis([0 350 -30 45])

% Interference cancellation
s = removePLI(x, fs, M, B, P, W);
subplot(2,1,2); pwelch(s(idx),[],[],[],fs);    %PSD of the filtered signal
title('ECoG, PSD of signal AFTER interference cancellation')
axis([0 350 -30 45])

%% EEG
load 'EEG.mat';

fs = 160; %Sample rate (Hz)
x = data; % Scaling to mV

M = 1;  %number of harmonics to remove
B = [50 .05 5] ;   
P = [0.001 2 5] ;   
W = 3;

idx = fs:length(x); % portion of signal to be used for PSD estimation

figure; movegui(gcf,'west')
subplot(2,1,1)
pwelch(x(idx),[],[],[],fs);    %PSD of the original signal
title('EEG, PSD of signal BEFORE interference cancellation')
axis([0 80 -30 45])

% Interference cancellation
s = removePLI(x, fs, M, B, P, W);
subplot(2,1,2); pwelch(s(idx),[],[],[],fs);    %PSD of the filtered signal
title('EEG, PSD of signal AFTER interference cancellation')
axis([0 80 -30 45])

%% ECG
load 'ECG.mat';

fs = 250; %Sample rate (Hz)
x = data; % Scaling to mV

M = 3;  %number of harmonics to remove
B = [50 .05 1];
P = [0.01 2 1];
W = 1.8;

idx = fs:length(x); % portion of signal to be used for PSD estimation

figure; movegui(gcf,'center')
subplot(2,1,1)
pwelch(x(idx),[],[],[],fs);    %PSD of the original signal
title('ECG, PSD of signal BEFORE interference cancellation')
%axis([0 350 -15 55])

% Interference cancellation
s = removePLI(x, fs, M, B, P, W);
subplot(2,1,2); pwelch(s(idx),[],[],[],fs);    %PSD of the filtered signal
title('ECG, PSD of signal AFTER interference cancellation')
%axis([0 350 -15 55])

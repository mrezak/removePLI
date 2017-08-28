function s = removePLI(x, fs, M, B, P, W, f_ac)
%removePLI Power Line Interference Cancellation 
%   This is an implementation of the proposed algorithm in,
%   M. R. Keshtkaran and Z. Yang, "A fast, robust algorithm for power line 
%   interference cancellation in neural recording," J. Neural Eng., vol. 11,
%   no. 2, p. 026017, Apr. 2014.
%   http://iopscience.iop.org/1741-2552/11/2/026017
%   http://arxiv.org/abs/1402.6862
%
%	Usage:
%	s = removePLI(x, fs, M, B, P, W, f_ac)
%   x, input (contaminated) signal
%	s, output (clean) signal
%   fs, sample rate in Hz
%   M, number of harmonics to remove
%   B, contains three elements [B0,Binf,Bst]: 
%	- B0, Initial notch bandwidth of the frequency estimator
%	- Binf, Asymptotic notch bandwidth of the frequency estimator
%	- Bst, Rate of convergence to 95% of the asymptotic bandwidth Binf
%   P, contains three elements [P0,Pinf,Pst]: 
%	- P0, Initial settling time of the frequency estimator
%	- Pinf, Asymptotic settling time of the frequency estimator
%	- Pst, Rate of convergence to 95% of the asymptotic settling time
%	W, Settling time of the amplitude and phase estimator
% 	f_ac, Optional argument, the nominal AC frequency if known (50 Hz or 60 HZ)
%
%	EXAMPLE:
%		fs = 500;
%		n = 120*fs; %2-min sequence	
%		t = 2*pi*(1:n)/fs;
%		fline = 60 + randn; %ramdom interference frequency
%		s = filter(1,[1,-0.99],100*randn(1,n)); %1/f PSD
%		p = 80*sin(fline*t+randn) + 50*sin(2*fline*t+randn)...
%		  + 20*sin(3*fline*t+randn); % interference	
%		x = s + p;
% 		sbar = removePLI(x, fs, 3, [100,0.01,4], [0.1,2,5], 3);
% 		pwelch(s,[],[],[],fs); title('PSD of the original signal')
% 		figure; pwelch(x(fs:end),[],[],[],fs); 
%       title('PSD of the contaminated signal');
% 		figure; pwelch(sbar(fs:end),[],[],[],fs); 
%      title('PSD after interference cancellation');
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

x = x(:)' - mean(x); %removing the mean
N = length(x);
s = zeros(1,N);

% 3dB cutoff bandwidth
alpha_f = (1-atan(pi*B(1)/fs))/(1+atan(pi*B(1)/fs));	%initial, \alpha_0
alpha_inf = (1-tan(pi*B(2)/fs))/(1+tan(pi*B(2)/fs)); %asymptotic	
alpha_st = exp(log(0.05)/(B(3)*fs+1));	%rate of change

% frequency estimator's forgetting factors
lambda_f = exp(log(0.05)/(P(1)*fs+1));	%initial
lambda_inf = exp(log(0.05)/(P(2)*fs+1));	%asymptotic
lambda_st = exp(log(0.05)/(P(3)*fs+1));	%rate of change
% Smoothing parameter (cut-off freq set at 90 Hz)
gmma = (1-tan(0.5*pi*min(90,fs/2)/fs))/(1+tan(0.5*pi*min(90,fs/2)/fs));

% phase/amplitude estimator forgetting factor
lambda_a = exp(log(0.05)/(W(1)*fs+1));
if(length(lambda_a)==1), lambda_a =lambda_a*ones(1,M); end

% initializing variables
kappa_f = 0;
kappa_k = zeros(1,M);
D=10; C=5;
f_n1=0; f_n2=0;

% -- Alternative initialization:
%    kappa_f = cos(55*2*pi/fs);
% --

% initializing the first oscillator
u_kp = 1*ones(1,M); %u_k
u_k = 1*ones(1,M); %u'_k

%initializing the RLS parameters
r1 = 10*ones(1,M);
r4 = 10*ones(1,M);
a = zeros(1,M);
b = zeros(1,M);


% IIR Bandpass filtering:
if(nargin>6)      % if AC frequency is known
    if length(f_ac)==2
        Fc1 = f_ac(1);  % First Cutoff Frequency
        Fc2 = f_ac(2);  % Second Cutoff Frequency               
     else
        %Custom center frequency of pass band
        Fc1 = f_ac-2;  % First Cutoff Frequency
        Fc2 = f_ac+2;  % Second Cutoff Frequency
    end
else    %if AC frequency is not known
    %Default 40--70 Hz pass band
    Fc1 = 40;  % First Cutoff Frequency
    Fc2 = 70;  % Second Cutoff Frequency
    f_ac = [];
end

ordr   = 4;   % Order
h  = fdesign.bandpass('N,F3dB1,F3dB2', ordr, Fc1, Fc2, fs);
Hd = design(h, 'butter');
x_f = filter(Hd,x);	%Bandpass Filtering
x_f = [0 diff(x_f)];		%First Difference

%--------- Start of data processing:
for n=1:N
	% Lattice Filter
    f_n = x_f(n) + kappa_f*(1+alpha_f)*f_n1 - alpha_f*f_n2;
    
	% Frequency Estimation
    C = lambda_f*C+(1-lambda_f)*f_n1*(f_n+f_n2);
    D = lambda_f*D+(1-lambda_f)*2*f_n1^2;
    kappa_t=C/D;
    if kappa_t <-1, kappa_t=-1; end
	if kappa_t > 1, kappa_t= 1; end
    kappa_f = gmma*kappa_f + (1-gmma)*kappa_t;
    
    f_n2=f_n1; f_n1=f_n; % Updating lattice states 

    % Bandwidth and Forgetting Factor Updates
	alpha_f = alpha_st*alpha_f + (1-alpha_st)*alpha_inf;
    lambda_f = lambda_st*lambda_f + (1-lambda_st)*lambda_inf;
    
	% Discrete-Time Oscillators
    kappa_k(2) = 1; kappa_k(1) = kappa_f;

    e=x(n);
    for k=1:M %for each harmonic do:
        %calculating Cos(kw) for k=1,2...
        kappa_k(k+2) = 2*kappa_f*kappa_k(k+1) - kappa_k(k); 
        
		% Oscillator
        tmp = kappa_k(k+2)*(u_kp(k)+u_k(k));
        tmp2 =u_kp(k);
        u_kp(k) = tmp - u_k(k);
        u_k(k) = tmp + tmp2;
        
		% Gain Control
        G = 1.5 - (u_kp(k)^2 - (kappa_k(k+2)-1)/(kappa_k(k+2)+1)*u_k(k)^2);
        if G<=0, G=1;end
        u_kp(k) = G * u_kp(k);
        u_k(k) = G * u_k(k);

        % Phase/Amplitude Adaptation
        sincmp = a(k)*u_k(k) + b(k)*u_kp(k);
        e = e - sincmp;
        %--- Simplified RLS
        r1(k) = lambda_a(k)*r1(k) + u_k(k)^2;
        r4(k) = lambda_a(k)*r4(k) + u_kp(k)^2;
        a(k) = a(k) + u_k(k)*e/r1(k);
        b(k) = b(k)  + u_kp(k)*e/r4(k); 
        %------
   end
   s(n)=e;
end

end

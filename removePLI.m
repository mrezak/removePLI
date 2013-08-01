function s = removePLI(x, fs, M, B, P, W)
%removePLI Power Line Interference Cancellation 
%   This is an implementation of the proposed algorithm in,
%  M. R. Keshtkaran, Z. Yang, "Power Line Interference Cancellation 
%	in Neural Recording", Submitted to Journal of Neural Engineering, Aug-2013
%
%	Usage:
%	s = removePLI(x, fs, M, B, P, W)
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
% 		sbar = removePLI(x, fs, 3, [50,0.1,1], [0.1,4,1], 4);
% 		pwelch(s,[],[],[],fs); title('PSD of the original signal')
% 		figure; pwelch(x,[],[],[],fs); title('PSD of the contaminated signal')
% 		figure; pwelch(sbar,[],[],[],fs); title('PSD of the cleaned signal')
%
%	Author: Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
%   Copyright (c) 2013, Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
%   All rights reserved.
%	This program is provided "AS IS" for non-commercial, educational 
%	and reseach purpose only. Any commercial use, of any kind, of 
%	this program is prohibited.
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

x=x(:)'; 
N = length(x);
s = zeros(1,N);

% 3dB cutoff bandwidth
bw = (1-tan(pi*B(1)/fs))/(1+tan(pi*B(1)/fs));	%initial
bwinf = (1-tan(pi*B(2)/fs))/(1+tan(pi*B(2)/fs)); %asymptotic	
bw0 = exp(log(0.05)/(B(3)*fs+1));	%rate of change

% frequency estimator's forgetting factors
lamF = exp(log(0.05)/(P(1)*fs+1));	%initial
lamFinf = exp(log(0.05)/(P(2)*fs+1));	%asymptotic
Ftc = exp(log(0.05)/(P(3)*fs+1));	%rate of change

% phase/amplitude estimator forgetting factor
lamA = exp(log(0.05)/(W(1)*fs+1));
if(length(lamA)==1), lamA =lamA*ones(1,M); end

% initializing variables
k0 = cos(55*2*pi/fs); 
b0=0.0;b1=0.0; D=10; C=D*k0;
b00=1;
% initializing the first oscillator
xs = ones(1,M);
ys = 2*ones(1,M);

%initializing the RLS parameters
Ha = 1000*ones(1,M);
Hb = 1000*ones(1,M);
a = ones(1,M);
b = ones(1,M);
hf = zeros(1,M);

% 40-70 Hz IIR filtering:
ordr   = 4;   % Order
Fc1 = 40;  % First Cutoff Frequency
Fc2 = 70;  % Second Cutoff Frequency
h  = fdesign.bandpass('N,F3dB1,F3dB2', ordr, Fc1, Fc2, fs);
Hd = design(h, 'butter');
X = filter(Hd,x);	%Bandpass Filtering
X=[0 diff(X)];		%First Difference

%--------- Start of data processing:
for n=1:N
	% Lattice Filter
    xfilt = X(n);   
    f1 = xfilt - bw*b1;
    f0 = f1 + k0*b0;
    b1 = -k0*f0 + b0;
    
	% Frequency Estimation
    D = lamF*D+(1-lamF)*2*b0^2;
    C = lamF*C+(1-lamF)*b0*(f0+b00);
    k0=C/D;
    if k0 <-1, k0=-1; end
	if k0 > 1, k0= 1; end
    b00=b0;
    b0=f0;
    
	% Bandwidth and Forgetting Factor Updates
	bw = bw0*bw + (1-bw0)*bwinf;
    lamF = Ftc*lamF + (1-Ftc)*lamFinf;
    
	% Discrete-Time Oscillators
    hf(2) =1; hf(1) = k0;
    e=x(n);
    for k=1:M %for each harmonic do:
        hf(k+2) = 2*k0*hf(k+1) - hf(k); %calculating Cos(kw) for k=1,2...
        
		% Oscillator
        tmp = hf(k+2)*(xs(k)+ys(k));
        xsp =xs(k);
        xs(k) = tmp - ys(k);
        ys(k) = tmp + xsp;
        
		% Gain Control
        G = 1.5 - (xs(k)^2 - (hf(k+2)-1)/(hf(k+2)+1)*ys(k)^2);
        if G<0, G=1;end
        xs(k) = G * xs(k);
        ys(k) = G * ys(k);

        % Phase/Amplitude Adaptation
        sincmp = a(k)*ys(k) + b(k)*xs(k);
        e = e - sincmp;
        %--- Simplified RLS
        Ha(k) = lamA(k)*Ha(k) + ys(k)^2;
        Hb(k) = lamA(k)*Hb(k) + xs(k)^2;
        a(k) = a(k) + ys(k)*e/Ha(k);
        b(k) = b(k)  + xs(k)*e/Hb(k); 
        %------
   end
   s(n)=e;
end

end

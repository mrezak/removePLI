Power-line-interference-canceller
=================================

The source code of a proposed algorithm to remove power line interference in bio-medical signals (Intracortical, ECoG, EEG, ECG (EKG), EMG)

The provided MATLAB .m file is the implementation of the proposed algorithm in 
Mohammad Reza Keshtkaran and Zhi Yang, "Power Line Interference Cancellation in Neural Recording", Submitted to Journal of Neural Engineering, 8/2013.

You need MATLAB Software to run the program. 

For User Guide, please refer to "removePLI.m" or type "help removePLI" in the MATLAB command prompt.

## Quick Start

	Usage:
```
  $ s = removePLI(x, fs, M, B, P, W)
```
  x, input (contaminated) signal
	s, output (clean) signal
  fs, sample rate in Hz
  M, number of harmonics to remove
  B, contains three elements [B0,Binf,Bst]: 
	- B0, Initial notch bandwidth of the frequency estimator
	- Binf, Asymptotic notch bandwidth of the frequency estimator
	- Bst, Rate of convergence to 95% of the asymptotic bandwidth Binf
  P, contains three elements [P0,Pinf,Pst]: 
	- P0, Initial settling time of the frequency estimator
	- Pinf, Asymptotic settling time of the frequency estimator
	- Pst, Rate of convergence to 95% of the asymptotic settling time
	W, Settling time of the amplitude and phase estimator
	
	EXAMPLE:
		fs = 500;
		n = 120*fs; %2-min sequence	
		t = 2*pi*(1:n)/fs;
		fline = 60 + randn; %ramdom interference frequency
		s = filter(1,[1,-0.99],100*randn(1,n)); %1/f PSD
		p = 80*sin(fline*t+randn) + 50*sin(2*fline*t+randn)...
		  + 20*sin(3*fline*t+randn); % interference	
		x = s + p;
		sbar = removePLI(x, fs, 3, [50,0.1,1], [0.1,4,1], 4);
		pwelch(s,[],[],[],fs); title('PSD of the original signal')
		figure; pwelch(x,[],[],[],fs); title('PSD of the contaminated signal')
		figure; pwelch(sbar,[],[],[],fs); title('PSD of the cleaned signal')


## Author
** Mohammad Reza Keshtkaran **

## Licence
   Copyright (c) 2013, Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
   All rights reserved.
   
   "This program" refers to "removePLI.m".
   This program is provided "AS IS" for non-commercial, educational 
 	 and reseach purpose only. Any commercial use, of any kind, of
	 this program is prohibited.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

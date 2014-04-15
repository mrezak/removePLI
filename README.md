Adaptive Power Line Interference Canceller
=================================

The source code of a fast and robust algorithm for removal of power line interference from biomedical signals.
The provided MATLAB .m file is the implementation of the proposed algorithm in:
**Mohammad Reza Keshtkaran and Zhi Yang, "A Fast, Robust Algorithm for Power Line Interference Cancellation in Neural Recording", Journal of Neural Engineering 11 026017, available at http://iopscience.iop.org/1741-2552/11/2/026017**

If you find this program useful in your research, please give credit by citing the above paper. If you have any question regarding the algorithm or implementation, do not hesitate to write to the authors at one of the following addresses: keshtkaran AT ieee.org, keshtkaran AT nus.edu.sg

You need MATLAB software to use this program. 

For User Guide, please refer to "removePLI.m" or type "help removePLI" in the MATLAB command prompt.
The MATLAB script "test_removePLI.m" runs the algorithm on sample intracortical, ECoG, EEG and ECG signals.
Please refer to "removePLI_multichan.m" for an optimized implementation for multichannel data.

You can download the full package by clicking "Download Zip" on the right side of this page.

## Usage
For Graphical User Interface:
```
>> removePLI_gui
```
For running from the command line:
```
>> s = removePLI(x, fs, M, B, P, W, f_ac)
>> s = removePLI_multichan(x, fs, M, B, P, W, f_ac, freqChan)
```
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
  f_ac, Optional argument, the nominal AC frequency if known (50 Hz or 60 HZ)
  freqChan, Optional argument, The channel number to be used for frequency
     estimation. []: uses the first channel (default), 0: Performs frequency
     estimation indivisually for all the channels (similar to removePLI)
```
```
  EXAMPLE:
		fs = 500;
		n = 60*fs; %1-min sequence	
		m = 10; %number of channels
		t = 2*pi*(1:n)/fs;
		fline = 60 + randn; %ramdom interference frequency
		s = filter(1,[1,-0.99],100*randn(n,m))'; %1/f PSD
		p = bsxfun(@times,sin(fline*t+randn), (80+5*randn(m,1))) ...
          + bsxfun(@times,sin(2*fline*t+randn), (50+5*rand(m,1))) ...
		  + bsxfun(@times,sin(3*fline*t+randn), (20+5*randn(m,1))); % interference	
		x = s + p;
 		sbar = removePLI_multichan(x, fs, 3, [50,0.01,4], [0.1,2,4], 2);
 		pwelch(s(1,:),[],[],[],fs); title('PSD of the original signal')
 		figure; pwelch(x(1,fs:end),[],[],[],fs); 
		title('PSD of the contaminated signal');
 		figure; pwelch(sbar(1,fs:end),[],[],[],fs); 
		title('PSD after interference cancellation');
```

## Author
**Mohammad Reza Keshtkaran**
## Licence
  Copyright (c) 2013, Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
  
  All rights reserved.
  
  "This program" refers to the m-files in the whole package.
  This program is provided "AS IS" for non-commercial, educational 
  and reseach purpose only. Any commercial use, of any kind, of 
  this program is prohibited. The Copyright notice should remain intact.

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

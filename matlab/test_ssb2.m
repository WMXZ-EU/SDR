% Define a message signal which contains three tones at 500, 600,
% and 700 Hz with varying amplitudes.
Fs = 10e3;
dt = 1/Fs;
t = (0 : dt : .1-dt)';
m = sin(2*pi*500*t) + .5*sin(2*pi*600*t) + 2*sin(2*pi*700*t) + 1*sin(2*pi*1500*t);

fo = Fs/4; % Carrier frequency in Hz

% generate SSB signal relative to carrier frequency fo
isb=1; % -1 is LSB, +1 is USB

mc = hilbert(m);
if isb==0
    mcm = real(mc.*exp(1i*2*pi*fo*t));
else
    mcm = real(mc.*exp(1i*isb*2*pi*fo*t));
end

% demodulate SSB signal
if isb == 0
    bp=zeros(129,1); bp(65)=1;
else
    bp=fir1(128, ([400 800])/(Fs/2))';
end
bm = fft(bp,N);

if 0
    %time-based processing
    mf = hilbert(mcm);
    fi = mf.*exp(-1i*2*pi*fo*t); % baseband shift
    fi = fftfilt(bp,fi);
else
    % frequency-based processing
    N = length(m);
    fn = round(fo*dt*N);
    ish = mod((1:N)+fn+N-1,N)+1; % indices for base-band shifting
    
    % filter operation (combined hilbert transform and baseband shift)
    fs = fft(mcm);
    fs(2:N/2) = 2*fs(2:N/2); 
    fs(N/2+1:end)=0;
    fs = fs(ish);
    fs = fs .* bm;
    fi = ifft(fs);
end

f = real(fi); % is final signal 
a = abs(fi);  % is envelope of signal, which may be useful for some detection operation

figure(1)
subplot(211)
periodogram(f, [], 4096, Fs, 'power', 'centered');
ylim([-75,0]);
subplot(212)
plot(1:N,m,(1:(N-64)),f(65:N));

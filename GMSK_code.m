
clc
clear
% N la chieu dai chuoi bit DATA nguoi dung
N = 10;
% Tan so lay mau ung voi 1 hinh sin
f_sample = 1024;

InputBit = randi(2,1,N) - ones(1,N);



% Lay bit le
OddBit = zeros(1,N);
for i = 1:2:N
    OddBit(i:i+1) = [InputBit(i) InputBit(i)];
end


% Lay bit chan
EvenBit = zeros(1,N);
for i = 2:2:N
    EvenBit(i) = InputBit(i);
    if i+1<=N
        EvenBit(i+1) = InputBit(i);
    end
end



% Tim bit to hop tu bit chan va bit le (SGT thay Vu Duc Tho trang 118)
CompositeBit = zeros(size(OddBit));
for i  = 1:N
    if OddBit(i) == EvenBit(i)
        CompositeBit(i) = 0;
    else
        CompositeBit(i) = 1;
    end
end
if 1
    disp('Du lieu goc');
    disp(InputBit); % Hien thi chuoi bit DATA
    disp('Odd Bit goc');
    disp(OddBit);
    disp('Even Bit goc');
    disp(EvenBit);
    disp('Composite goc');
    disp(CompositeBit);
end

%%
% Chuyen sang tin hieu dien
f_low = 10;%(135.4 - 67.7) * 10^3;
f_high = 30;%(135.4 + 67.7) * 10^3;



% 1 bit ung voi truc thoi gian
time_bit = 0:1/f_sample:1;

% Tao tin hieu hinh sin
sin_low = sin(2*pi*f_low*time_bit);
sin_high = sin(2*pi*f_high*time_bit);




size_bit = size(sin_low,2); % chieu dai 1 bit
zero_signal = [];


% Ve tin hieu mau ung voi tan so cao va tan so thap
if 0
    plot(time_bit,sin_low,'b');hold on;
    plot(time_bit,sin_high,'r');hold off;
end


% Mot bit ung voi truc thoi gian
time_signal = 0:1/f_sample:N;


% Tin hieu dien ung voi chuoi bit



Composite_signal = zero_signal;

for i = 1:N
    size_temp = size(Composite_signal,2);
    if CompositeBit(i) == 0
        Composite_signal = [Composite_signal(1:size_temp-1) zeros(size(time_bit))];
    else
        Composite_signal = [Composite_signal(1:size_temp-1) ones(size(time_bit))];
    end
end

% Ve chuoi bit to hop tu bit chan va bit le
if 0
    plot(time_signal,Composite_signal);
    hold on;
end



% Tin hieu dieu che MSK
MSK_signal = zero_signal;

for i = 1:N
    size_temp = size(MSK_signal,2);
    if CompositeBit(i) == 0
        MSK_signal = [MSK_signal(1:size_temp-1) sin_high];
    else
        MSK_signal = [MSK_signal(1:size_temp-1) sin_low];
    end
end

% Ve tin hieu dieu che MSK
if 0
    plot(time_signal,MSK_signal);
end



% Phan tich pho MSK

fs = 5*f_sample; % Tan so lay mau
N_s = fs; % FFT N diem
transform = fft(MSK_signal,N_s); % Chuan hoa 
magTransform = abs(transform); % Lay bien do

% Ve pho MSK
if 0
    figure;
    faxis = linspace(-fs/2,fs/2,N_s); % Tao tin hieu truc tan so
    plot(faxis,fftshift(magTransform),'r'); % Ve pho FFT 
end



%% Gaussian Filter

bt = 0.3;
span = 4; % span
sps = 16;




Gauss_Filter = gaussdesign(bt,span,sps);

GMSK_signal = filter(Gauss_Filter,1,MSK_signal);

% Ve tin hieu GMSK va MSK o mien thoi gian de so sanh
% MSK : Red line
% GMSK: Blue line

if 0
    figure;
    plot(time_signal,MSK_signal,'r');
    hold on
    plot(time_signal,GMSK_signal,'b');
end

% Phan tich pho GMSK

GMSK_amp = abs(fft(GMSK_signal,N_s));


% Ve pho GMSK va MSK de so sanh
% MSK : Red line
% GMSK: Blue line
if 0
    figure;
    faxis = linspace(-fs/2,fs/2,N_s); % Tao tin hieu truc tan so
    plot(faxis,fftshift(magTransform),'r'); % Ve pho FFT
    hold on
    plot(faxis,fftshift(GMSK_amp),'b'); % Ve pho FFT
end













%% Channel






Receive_signal = awgn(GMSK_signal,1);

if 1
    figure;
    plot(time_signal,Receive_signal,'r')
    %hold on;
    figure;
    plot(time_signal,GMSK_signal,'b');
 
end


Receive_amp = abs(fft(Receive_signal,N_s));
if 1
    figure;
    faxis = linspace(-fs/2,fs/2,N_s); % Tao tin hieu truc tan so
    plot(faxis,fftshift(Receive_amp),'r'); % Ve pho FFT
    hold on
    plot(faxis,fftshift(GMSK_amp),'b'); % Ve pho FFT
end








%% Demodulation

Inph_cos = cos(2*pi*30*time_signal);
Quad_sin = sin(2*pi*30*time_signal);

Gain = 100;

[b, a] = butter(5,[1/512 90/512]);


LowPass_Inph = Gain*filter(b,a,Receive_signal.*Inph_cos);
LowPass_Quad = Gain*filter(b,a,Receive_signal.*Quad_sin);

Sum_signal = LowPass_Inph+LowPass_Quad;

Sum_amp = abs(fft(Sum_signal,N_s))/Gain;
if 0
    figure;
    faxis = linspace(-fs/2,fs/2,N_s); % Tao tin hieu truc tan so
    plot(faxis,fftshift(Receive_amp),'r'); % Ve pho FFT
    hold on
    plot(faxis,fftshift(Sum_amp),'b'); % Ve pho FFT
end


% Ve chuoi bit to hop tu bit chan va bit le
if 0
    plot(time_signal,150*Composite_signal,'r');
    hold on;
    plot(time_signal,abs(Sum_signal),'b');
end

Sum_signal_abs = abs(Sum_signal);

% Bit Reconstruction

Reconstruction_bit = [];

Threshold = sum(Sum_signal_abs)/(N-1);



for i = 1:size_bit-1:(N-1)*(size_bit-1)+1
    tmp = sum(Sum_signal_abs(i:i+size_bit-2));
    if tmp >= Threshold
        Reconstruction_bit = [Reconstruction_bit 1];
    else
        Reconstruction_bit = [Reconstruction_bit 0];
    end
end

if 0
    disp(CompositeBit);
    disp(Reconstruction_bit);
end

% Recover Odd bit and Even bit

Recover_OddBit = zeros(1,N);
Recover_EvenBit = zeros(1,N);

Recover_OddBit(1) = Reconstruction_bit(1);
Recover_OddBit(2) = Reconstruction_bit(1);

for i = 2:2:N
    Recover_EvenBit(i) = int8(xor(Recover_OddBit(i),Reconstruction_bit(i)));
    if i+1 < N 
        Recover_EvenBit(i+1) = Recover_EvenBit(i);
        Recover_OddBit(i+1) = int8(xor(Recover_EvenBit(i+1),Reconstruction_bit(i+1)));
        Recover_OddBit(i+2) = Recover_OddBit(i+1);
    end
end

OutputBit = zeros(1,N);

for i = 1:2:N
    OutputBit(i:i+1) = [Recover_OddBit(i) Recover_EvenBit(i+1)] ;
end



if 1
    disp('Composite thu');
    disp(Reconstruction_bit);
    disp('Even Bit thu');
    disp(Recover_EvenBit);
    disp('Odd Bit thu');
    disp(Recover_OddBit);
    disp('Du lieu thu');
    disp(OutputBit); % Hien thi chuoi bit DATA
end


% Tim BER

ErrorBit_count = sum(abs(OutputBit-InputBit));

BER = ErrorBit_count*1.0/N;
if 1
    disp('------------------------');
    disp('BER');
    disp(BER);
end


%% Post filtering with wiener filter
% For each segment
% Calculate Px 
% Calculate Hw = 1 - Pu / Px
% Apply Hw 
% Invert and OLA

function yw = post_filter(y, Pu, L, WINDOW, NOVERLAP, NFFT, fs)

w = hamming(NFFT);

index = 1;
yw(1 : length(y) + L) = 0;
ylen = length(y);

while index + NFFT - 1 <= ylen     
    win_frame = w .* y(index : index+NFFT-1);
    Px = pwelch(win_frame, WINDOW, NOVERLAP, NFFT, fs, 'twosided');
    Hw = 1 - Pu./Px;
    Xw = fft(win_frame, NFFT);
    Yw = Hw .* Xw;
    yw_frame = ifft(Yw);
    yw(index : index + NFFT - 1)= yw(index : index + NFFT - 1)+ yw_frame';
    index = index+ L / 2;
    
end
    yw = real(yw(1:ylen));
end


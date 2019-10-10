%Ερώτημα 1.6

function Vector = ttdecode(signIn)
digits = zeros(800,1200);       %initialize array with the right digits according to each frequency range
digits(1:550,1:970)=1;
digits(1:550,971:1060)=2;
digits(1:550,1061:1200)=3;
digits(551:630,1:970)=4;
digits(551:630,971:1060)=5;
digits(551:630,1061:1200)=6;
digits(631:700,1:970)=7;
digits(631:700,971:1060)=8;
digits(631:700,1061:1200)=9;
digits(701:800,971:1060)=0;

Vector = zeros(1,floor(length(signIn)/1000));
k = 1;
for i = 1:1100:length(signIn)
    if i+999>length(signIn)
        f = fft(signIn(i:length(signIn)));
    else
        f = fft(signIn(i:i+999));   %DFT of every tone of the signal
    end
    e = f.*conj(f);             %Energy around every frequence   
    fbetween = 0.8/(2*pi);      %0.8 rad/sec is a good threshold (according to the array given
                                %to seperate high from low frequencies
    fstop = 1.14/(2*pi);        %maximum frequency is 1.1328 Hz, we give it some range until 1.14 Hz
    start = fbetween*length(e); %find index of matrix e where higher frequencies begin
    start=floor(start);         
    stop = fstop*length(e);     %find index of matrix e where we have to stop 
                                %(after that energy is close to zero and at some point begin the negative frequencies)
    stop = floor(stop);
    flo = find(e(1:start)==max(e(1:start)));    %find low frequency with max energy
    fh = find(e(start:stop)==max(e(start:stop)));   %find high frequency with max energy
    fh = fh + start; %calculate correct frequency (find(e(start:stop)==max(e(start:stop))) starts from zero index)
    flo = 2*pi*flo/length(e);   
    fh = 2*pi*fh/length(e);
  
    indexh = floor(fh*length(e));
    indexlo = floor(flo*length(e));
    Vector(k) = digits(indexlo,indexh);
    k = k+1;
    
end
disp(Vector);
end
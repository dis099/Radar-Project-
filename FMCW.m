[Y,FS] = audioread('Umer_Range.m4a');
c = 3e8;
f_start=2.4e9;
f_stop=2.5e9;
BW=f_stop-f_start;
T_pulse = 20e-3;
N = FS*T_pulse;
range_resolution=c/(2*BW);
range_max=range_resolution*N/2;

trig=Y(:,2);
data=Y(:,1);

%f=(1:length(Y))*FS;
%plot(f,trig);

time=zeros(length(Y),1);

value=0;
tresh=0;
start=(trig>tresh);

for j = 100:(size(start,1)-N)
    if start(j) == 1 && mean(start(j-11:j-1)) == 0
        value = value + 1; 
       
    end
end

s=zeros(value,N);
count=0;
zpad=4*N;


for k = 100:(size(start,1)-N)
    if start(k) == 1 && mean(start(k-11:k-1)) == 0
        count = count + 1; 
        s(count,:) = data(k:k+N-1); 
        time(count) = k*1/FS;
    end
end

ave = mean(s,1);
for l = 1:size(s,1)
    s(l,:) = s(l,:) - ave;
end

FFT=fft(s,zpad,2);
v=20*log10(abs(FFT));
v1=v(:,1:size(v,2)/2);
v2=v1-max(max(v1));
range=linspace(0,range_max,zpad);
clims = [max(max(v2))-50 max(max(v2))-0];
imagesc(range,time(1:size(v2,1)),v2,clims);
xlabel('Range (m)');
ylabel('Time (sec)');
axis([0 100 -Inf Inf]);
colorbar;

figure()
s2 = s(2:size(s,1),:)-s(1:size(s,1)-1,:);
FFT=fft(s2,zpad,2);
v=20*log10(abs(FFT));
v1=v(:,1:size(v,2)/2);
v2=v1-max(max(v1));
range=linspace(0,range_max,zpad);
clims = [max(max(v2))-50 max(max(v2))-0];
imagesc(range,time(1:size(v2,1)),v2,clims);
xlabel('Range (m)');
ylabel('Time (sec)');
axis([0 100 -Inf Inf]);
colorbar;
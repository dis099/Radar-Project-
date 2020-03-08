[Y,FS] = audioread('Umer_Running.m4a');
c = 3e8;
f_0 = 2.424e9;
T_pulse = 100e-3;
N = FS*T_pulse;
number_of_row=floor(length(Y)/N)+1;
Array=zeros(number_of_row,N);

Mean_Y=mean(mean(Y));

for i = 1:number_of_row
    for j = 1:N
        Array(i,j)=Y((i-1)*N+j)-Mean_Y;
    end
end

FFT=fft(Array,4*N,2);

v=20*log10(abs(FFT));

v1=v(:,1:size(v,2)/2);

v2=v1-max(max(v1));
delta_f=linspace(0,FS/2,size(v2,2));
velocity=(delta_f*c)/(2*f_0);
time = linspace(1, T_pulse*size(v,1),size(v,1));
imagesc(velocity,time,v2,[(max(max(v2))-35) (max(max(v2))-0)]);

colormap('default'); 
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
axis([0 40 -Inf Inf]);
colorbar
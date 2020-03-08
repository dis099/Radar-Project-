%FMCW
c = 3e8;
f_0 = 2.424e9;
T_pulse = 100e-3;
FS=44100;
N = FS*T_pulse;

%SamplesPerFrame=1024;
%Microphone = dsp.AudioFileReader('Nouvel enregistrement 4.m4a');

SamplesPerFrame=FS;
Microphone = dsp.AudioRecorder('SampleRate',FS,'NumChannels',2,...
    'OutputNumOverrunSamples',true,'SamplesPerFrame',SamplesPerFrame);
test_condition=1;

figure();
clims=[-30 0];
time_data_selection=1;
seconde_display=15;
zpad=4*N;

Z=step(Microphone);
Y=zeros(time_data_selection*FS,2);
number_of_row=floor(length(Y)/N)+1;
Array=zeros(number_of_row,N);
seconde=0;

while(test_condition)
    i=0;
    while(i<(time_data_selection*FS)/SamplesPerFrame)
        Y(i*SamplesPerFrame+1:(i+1)*SamplesPerFrame,:)=Z;
        Z=step(Microphone);
        i=i+1;
    end
    
    Mean_Y=mean(mean(Y));
    
    for i = 1:number_of_row
        for j = 1:N
            Array(i,j)=Y((i-1)*N+j)-Mean_Y;
        end
    end
    
    FFT=fft(Array,zpad,2);
    v=20*log10(abs(FFT));
    v1=v(:,1:size(v,2)/2);
    v2=v1-max(max(v1));
    delta_f=linspace(0,FS/2,size(v2,2));
    velocity=(delta_f*c)/(2*f_0);
    time = time_data_selection*seconde+linspace(0, T_pulse*size(v,1),size(v,1));
    seconde=seconde+1;
    imagesc(velocity,time,v2,clims);
    hold on
    xlabel('Velocity (m/sec)');
    ylabel('Time (sec)');
    axis([0 40 time_data_selection*(seconde-seconde_display/time_data_selection) Inf]);
    colorbar;
    drawnow;
    Y=zeros(time_data_selection*FS,2);
end

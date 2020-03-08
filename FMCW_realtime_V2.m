%FMCW
c = 3e8;
f_start=2.4e9;
f_stop=2.5e9;
BW=f_stop-f_start;
T_pulse = 20e-3;
FS=44100;
N = FS*T_pulse;
range_resolution=c/(2*BW);
range_max=range_resolution*N/2;

%SamplesPerFrame=1024;
%Microphone = dsp.AudioFileReader('Mesure/Nouvel enregistrement 7 copie.m4a');

SamplesPerFrame=1024;
Microphone = dsp.AudioRecorder('SampleRate',FS,'NumChannels',2,...
    'OutputNumOverrunSamples',true,'SamplesPerFrame',SamplesPerFrame);
test_condition=1;

figure();
clims=[-30 0];
time_data_selection=1;
seconde_display=15;

zpad=4*N;

Y=step(Microphone);
Z=zeros(time_data_selection*FS,2);
seconde=0;

while(test_condition)
    i=0;
    value=0;
    while(i<(time_data_selection*FS)/SamplesPerFrame)
        Z(i*SamplesPerFrame+1:(i+1)*SamplesPerFrame,:)=Y;
        Y=step(Microphone);
        i=i+1;
    end
    
    trig=Z(:,2);
    data=Z(:,1);
    time=zeros(length(Z),1);
    tresh=0;
    start=(trig>tresh);
    for j = 100:(size(start,1)-N)
        if start(j) == 1 && mean(start(j-11:j-1)) == 0
            value = value + 1; 
        end
    end

    s=zeros(value,N);
    count=0;


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
    s2 = s(2:size(s,1),:)-s(1:size(s,1)-1,:);
    FFT=fft(s2,zpad,2);
    v=20*log10(abs(FFT));
    v1=v(:,1:size(v,2)/2);
    MAX=max(max(v1));
    v2=v1-MAX;
    range=linspace(0,range_max,zpad);
    %clims = [max(max(v2))-35 max(max(v2))-0];
    time1=time_data_selection*seconde+time(1:size(v2,1));
    seconde=seconde+1;
    imagesc(range, time1,v2,clims);
    hold on
    xlabel('Range (m)');
    ylabel('Time (sec)');
    %axis([0 30 -Inf Inf]);
    axis([0 40 time_data_selection*(seconde-seconde_display/time_data_selection) Inf]);
    colorbar;
    drawnow;
    Z=zeros(time_data_selection*FS,2);

end

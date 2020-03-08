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
t=zeros(value-1,N);
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

s2 = s(2:size(s,1),:)-s(1:size(s,1)-1,:);
FFT=fft(s2,zpad,2);
v=20*log10(abs(FFT));
v1=v(:,1:size(v,2)/2);

v3=sgolayfilt(v1,3,11);
v4=v3-max(max(v3));
[max2,indice2]=max(v4,[],2);
vitessefilt=zeros(size(v4,1)-1,1);
vitessefilt1=zeros(size(v4,1)-1,1);

v2=v1-max(max(v1));

% frequency=linspace(0,FS/2,size(v2,2));
% plot(frequency,v1(50,:));
% hold on;
% plot(frequency,v2(50,:));
% hold on;
% plot(frequency,v4(50,:));

[max1,indice1]=max(v2,[],2);
%vitesse=zeros(size(v2,1)-1,1);
time = linspace(0,length(trig)/FS,size(v2,1));
range=linspace(0,range_max,zpad);
frequence=linspace(0,range_max,zpad);
tabrange = zeros(size(v2,1),1);
tabrange1 = zeros(size(v2,1),1);
for o = 1 : size(v2,1)-1
    tabrange(o)=range(indice1(o));
    tabrange1(o)=range(indice2(o));
end
tabrange2=sgolayfilt(tabrange,3,11);

for o = 1 : size(v2,1)-1
%     vitesse(o,1)=(tabrange2(o)-tabrange2(o+1));
    vitessefilt(o,1)=(range(indice2(o))-range(indice2(o+1)))*FS;
    vitessefilt1(o,1)=(tabrange1(o)-tabrange1(o+1))*FS ;
end

vitesse1=zeros(size(v2,1)-2,1);
for o = 2 : size(v2,1)-1
     vitesse1(o,1)=(tabrange2(o-1)-tabrange2(o+1))*FS/2;
end
% clims = [max(max(v2))-50 max(max(v2))-0];
% imagesc(range,time(1:size(v2,1)),v2,clims);
% xlabel('Range (m)');
% ylabel('Time (sec)');
% axis([0 100 -Inf Inf]);
% colorbar; 
% figure();
% plot(vitesse,time(1:size(v2)-1));
% vitess=sgolayfilt(vitesse,15,101);
% figure();
% plot(vitess,time(1:size(v2)-1));
vitf=sgolayfilt(vitessefilt,15,101);
% figure()
% plot(vitess,time(1:size(v2)-1));
% figure()
% plot(vitessefilt,time(1:size(v2)-1));
figure()
plot(vitf,time(1:size(v2)-1));
figure();
vitf1=sgolayfilt(vitessefilt1,15,101);
plot(vitf1,time(1:size(v2)-1));

figure();
vitdouble=sgolayfilt(vitesse1,15,101);
plot(vitdouble,time(1:size(v2)-1));

% for y = 1 :size(vitf,1)
%     for z = 1 : size(vitf,2)
%         if vitf(y,z)>3000
%             vitf(y,z)=0;
%         end
%         if vitf(y,z)<-3000
%             vitf(y,z)=0;
%         end
%     end
% end
% figure();
% plot(vitf,time(1:size(v2)-1));


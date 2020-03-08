[Y,FS] = audioread('Malvinas.m4a');
c = 3e8;
f=50e6;
f_c=2.45e9;
f_start=f_c-f;
f_stop=f_c+f;
T_pulse = 20e-3;
T_rp=250e-3;
R_s=0;

N=FS*T_pulse;
Nrp=FS*T_rp;
BW=f_stop-f_start;
lambda=c/f_c;
cr=BW/T_pulse;
delta_x=lambda/2;


%separate trigger and data
trig=Y(:,2);
data=Y(:,1);

rp_start=abs(trig)>mean(abs(trig));
value=0;
for j = Nrp+1:(size(rp_start,1)-Nrp)
    if rp_start(j) == 1 && sum(rp_start(j-Nrp:j-1)) == 0
        value = value + 1;
    end
end

RP=zeros(value,Nrp);
RP_trig=zeros(value,Nrp);
count=0;

for k = Nrp+1:(size(rp_start,1)-Nrp)
    if rp_start(k) == 1 && sum(rp_start(k-Nrp:k-1)) == 0
        count = count + 1;
        RP(count,:) = data(k:k+Nrp-1); 
        RP_trig(count,:) = trig(k:k+Nrp-1);
    end
end

thresh=0.08;
sif=zeros(size(RP,1),N/2);
for l= 1:size(RP,1)
    clear SIF;
    SIF = zeros(N,1); 
    start=(RP_trig(l,:)>thresh);
    count=0;
    for m = 12:(size(start,2)-2*N)
        [Y, I] = max(start(1,m:m+2*N));
        if mean(start(m-10:m-2)) == 0 && I == 1
            count = count + 1;
            SIF = RP(l,m:m+N-1)+SIF; 
        end
    end
    SI = SIF/count;
    FF = ifft(SI);
    clear SI;
    sif(l,:) = fft(FF(size(FF,1)/2+1:size(FF,1)));
end

for n=1:size(sif,1)
    for o=1:size(sif,2)
        if isnan(sif(n,o)) ==1 
            sif(n,o)=1e-30;
        end
    end
end
L = delta_x*(size(sif,1));
Xa=linspace(-L/2, L/2, (L/delta_x));
time = linspace(1, T_pulse, size(sif,2));
Kr = linspace(((4*pi/c)*(f_c - BW/2)), ((4*pi/c)*(f_c + BW/2)), (size(time,2)));
ave = mean(sif,1); 
for p = 1:size(sif,1)
    sif(p,:) = sif(p,:) - ave; 
end

clear N;
N = size(sif,2);
H=zeros(1,N);
for q = 1:N
    H(q) = 0.5 + 0.5*cos(2*pi*(q-N/2)/N); 
end

sif_h=zeros(size(RP,1),N);
for r = 1:size(sif,1)
    sif_h(r,:) = sif(r,:).*H;
end
sif = sif_h;

%figure 1
S_image = angle(sif);
imagesc(Kr, Xa, S_image); 
colormap('default'); 
xlabel('K_r (rad/m)'); 
ylabel('SAR Position, Xa (m)'); 
colorbar;

% zpad = 2048;
% s_zeros = zeros(zpad,size(sif,2)); 
% for i = 1:size(sif,2)
%     index = round((zpad - size(sif,1))/2);
%     s_zeros(index+1:(index + size(sif,1)),i) = sif(:,i); 
% end
% sif = s_zeros;
% S = fftshift(fft(sif,[],1),1);
% Kx = linspace((-pi/delta_x), (pi/delta_x), (size(S,1)));
% 
% %figure 2
% figure()
% S_image = 20*log10(abs(S));
% imagesc(Kr,Kx,S_image,[max(max(S_image))-40,max(max(S_image))]);
% colormap('default'); 
% xlabel('K_r (rad/m)'); 
% ylabel('K_x (rad/m)'); 
% colorbar;
% 
% %figure 3
% figure()
% S_image = angle(S); 
% imagesc(Kr, Kx, S_image); 
% colormap('default'); 
% xlabel('K_r (rad/m)'); 
% ylabel('K_x (rad/m)'); 
% colorbar;
% 
% S_matched = S;
% kstart=floor(min(Kr)-20);
% kstop=ceil(max(Kr));
% Ky_e=linspace(kstart,kstop,1024);
% count = 0;
% Ky=zeros(zpad,N);
% S_st=zeros(zpad,1024);
% for ii = 1:zpad
%     count = count + 1;
%     Ky(count,:) = sqrt(Kr.^2 - Kx(ii)^2);
%     S_st(count,:) = (interp1(Ky(count,:), S_matched(ii,:), Ky_e));
% end
% 
% for x=1:size(S_st,1)
%     for y=1:size(S_st,2)
%         if isnan(S_st(x,y))==1
%             S_st(x,y)=1e-30;
%         end
%     end
% end
% 
% %figure 4
% figure()
% S_image = angle(S_st); 
% imagesc(Ky_e, Kx, S_image); 
% xlabel('K_y (rad/m)'); 
% ylabel('K_x (rad/m)'); 
% colorbar;
% 
% v = ifft2(S_st,(size(S_st,1)*4),(size(S_st,2)*4));
% bw = c*(kstop-kstart)/(4*pi);
% max_range = (c*size(S_st,2)/(2*bw));
% figure(fig_count);
% S_image = v;
% S_image = fliplr(rot90(S_image));
% cr1 = -20; % depends on the Kx of the Stolt Interpolation 
% cr2 = 20; % depends on the Kx of the Stolt Interpolation 
% dr1 = 1;
% dr2=100;
% dr_index1 = round((dr1/max_range)*size(S_image,1)); 
% dr_index2 = round((dr2/max_range)*size(S_image,1)); 
% cr_index1 = round(((cr1+zpad*delta_x/(2*1))/(zpad*delta_x/1))*size(S_image,2));
% cr_index2 = round(((cr2+zpad*delta_x/(2*1))/(zpad*delta_x/1))*size(S_image,2));
% trunc_image = S_image(dr_index1:dr_index2,cr_index1:cr_index2); 
% downrange = linspace(-1*dr1,-1*dr2, size(trunc_image,1)); 
% crossrange = linspace(cr1, cr2, size(trunc_image, 2));
% for jj = 1:size(trunc_image,2)
%     trunc_image(:,jj) = (trunc_image(:,jj)').*(abs(downrange*1)).^(3/2); 
% end
% trunc_image = 20 * log10(abs(trunc_image));
% clims=[(max(max(trunc_image))-30) (max(max(trunc_image))-0)];
% 
% %figure 5
% figure()
% imagesc(crossrange, downrange, trunc_image,clims);
% colormap('default'); 
% axis equal;
% xlabel('Crossrange(m)');
% ylabel('Downrange (m)'); 
% colorbar;

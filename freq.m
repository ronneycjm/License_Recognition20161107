close all;
clear all;
%ģ���Ƕ�1~90  >90�޷�ʶ��  atan��-90~90��
%ģ������>5���� <5�޷�ʶ��
%�����ʮ�����⣬��edge canny����
%% ���벢��ʾͼ��
filename = 'test_xie.jpg';
I = imread(filename);
subplot(231)
imshow(uint8(I));
title('ԭͼ');

%% �����˶�ģ��ͼ��
PSF = fspecial('motion',100, 130);
g = imfilter(I, PSF, 'circular');
subplot(232)
imshow(uint8(g));
title('�˶�ģ��ͼ');

%% ���˶�ģ��ͼ����лҶȻ��������ж�ά���ٸ���Ҷ�任��������Ƶ��ͼ
gb = rgb2gray(g);
figure
imshow(uint8(gb));
PQ = paddedsize(size(gb));
F = fft2(gb);
subplot(233)
imshow(uint8(F));
title('Ƶ��ͼ');

%% ��Ƶ��ѹ��������
H = log(1+abs(F));
Hc = fftshift(H);
subplot(234)
imshow(uint8(Hc));
title('ѹ������');
%% ��canny���ӽ�ѹ�����к��Ƶ��ͼ���б�Ե��⣬��ֵ��
T = graythresh(Hc);
bw=edge(Hc, 'canny', T);
subplot(235)
imshow(bw);
title('canny���� ��ֵ��');

%% �Զ�ֵ�����Ƶ��ͼ����radon�任
thetastep = 0:0.1:180;
[R,xp] = radon(bw, thetastep);
subplot(236)
imshow(R);
title('radon�任');
figure 
imagesc(thetastep,xp,R);colormap(hot);colorbar;
%% �����ͨ��radon�任�����ģ���Ƕ�
MAX = max(max(R));
[m, n] = find(R == MAX);
[M,N] = size(Hc);
beita = (n-1)*0.1;%atan(tan(n*pi/180)*N/M)*180/pi;  %tan�ڻ���  atan�����ǻ��ȣ����ɽǶ�
alpha = beita + 90 ;
theta = atan ( tan(alpha*pi/180 - pi/2)*M/N ) *180/pi;
if theta < 0
    theta = theta +180;
end;
disp(theta);

%% �������
% figure;
% imshow(gb);
% gb_crop =gb(50:160,40:150);
% figure;
% imshow(gb_crop);

gb_rotate = imrotate(gb, -theta);%, 'crop');% 
figure;
imshow(gb_rotate);

[r,c] = size(gb_rotate);
ave = mean(mean(gb_rotate));
for ii = 1:r
    for jj = 1:c
        if gb_rotate(ii,jj) == 0
            gb_rotate(ii,jj)=150; %avg
        end
    end
end
figure;
imshow(gb_rotate);
dif = conv2(double(gb_rotate),[0.5,-0.5]);       %��ˮƽ�᷽���ϵ�һ��΢��ͼ��
dif(:,1) = 0;                               %����һ�к����һ����0����ֹ�߽�Ӱ��
dif(:,size(dif,2)) = 0;
for i = 1:size(dif,1)
    s(i,:) = smooth(xcorr(dif(i,:)));             %��difÿ�н������������
end
sum_s = sum(s);
 figure;
 plot(sum_s);
[n1,n2] = Find2Min(sum_s);    %�ҵ���С����ֵ���±�
dy= abs((n2 - n1) / 2);       %ģ���߶ȵ��ڹ�����ط�������2
disp(dy);

% %% Ƶ�׳��ȼ��
% %HcΪѹ��Ƶ��ͼ
% im_rotate_1 = wiener2(Hc,[2,2]);
% balance_im = histeq(uint8(im_rotate_1));
% level = graythresh(balance_im);  %T=level
% bw_im_rotate_1=im2bw(balance_im,level);
% subplot(232)
% imshow(im_rotate_1,[]);
% subplot(233)
% imshow(bw_im_rotate_1,[]);
% im_rotate = imrotate(bw_im_rotate_1,90-alpha);
% subplot(234)
% imshow(im_rotate,[]);
% s = sum(im_rotate);
% max_s = max(s);
% for num = 1:length(s)
%     if s(num)>(max_s/2)
%         s(num) = max_s;
%     end
% end
% mid =  (1+length(s))/2 ;
% for l1=floor(mid):-1:1
%     if s(l1+1) > s(l1) && s(sl1)<s(l1-1)
%         resleft = l1;
%         break;
%     end
% end
%     
% subplot(235)
% plot(s);

%% ά���˲�



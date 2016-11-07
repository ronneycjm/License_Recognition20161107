close all;
clear all;
%模糊角度1~90  >90无法识别  atan是-90~90！
%模糊长度>5均可 <5无法识别
%解决亮十字问题，用edge canny算子
%% 读入并显示图像
filename = 'test_xie.jpg';
I = imread(filename);
subplot(231)
imshow(uint8(I));
title('原图');

%% 生成运动模糊图像
PSF = fspecial('motion',100, 130);
g = imfilter(I, PSF, 'circular');
subplot(232)
imshow(uint8(g));
title('运动模糊图');

%% 对运动模糊图像进行灰度化，并进行二维快速傅里叶变换，生成其频谱图
gb = rgb2gray(g);
figure
imshow(uint8(gb));
PQ = paddedsize(size(gb));
F = fft2(gb);
subplot(233)
imshow(uint8(F));
title('频谱图');

%% 将频谱压缩，居中
H = log(1+abs(F));
Hc = fftshift(H);
subplot(234)
imshow(uint8(Hc));
title('压缩居中');
%% 用canny算子将压缩居中后的频谱图进行边缘检测，二值化
T = graythresh(Hc);
bw=edge(Hc, 'canny', T);
subplot(235)
imshow(bw);
title('canny算子 二值化');

%% 对二值化后的频谱图进行radon变换
thetastep = 0:0.1:180;
[R,xp] = radon(bw, thetastep);
subplot(236)
imshow(R);
title('radon变换');
figure 
imagesc(thetastep,xp,R);colormap(hot);colorbar;
%% 计算出通过radon变换求出的模糊角度
MAX = max(max(R));
[m, n] = find(R == MAX);
[M,N] = size(Hc);
beita = (n-1)*0.1;%atan(tan(n*pi/180)*N/M)*180/pi;  %tan内弧度  atan出来是弧度，换成角度
alpha = beita + 90 ;
theta = atan ( tan(alpha*pi/180 - pi/2)*M/N ) *180/pi;
if theta < 0
    theta = theta +180;
end;
disp(theta);

%% 求出长度
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
dif = conv2(double(gb_rotate),[0.5,-0.5]);       %求水平轴方向上的一阶微分图像
dif(:,1) = 0;                               %将第一列和最后一列置0，防止边界影响
dif(:,size(dif,2)) = 0;
for i = 1:size(dif,1)
    s(i,:) = smooth(xcorr(dif(i,:)));             %对dif每行进行自相关运算
end
sum_s = sum(s);
 figure;
 plot(sum_s);
[n1,n2] = Find2Min(sum_s);    %找到最小两个值的下标
dy= abs((n2 - n1) / 2);       %模糊尺度等于共轭相关峰距离除以2
disp(dy);

% %% 频谱长度检测
% %Hc为压缩频谱图
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

%% 维纳滤波



% 把模板高度变成车牌高度 162*620  即  115->162    115*60
clc;
clear;
%%read the image and get the height and width
I=imread('test_blur.jpg');
Ib = rgb2gray(I);
% I0=imread('0.jpg');
%imshow('im1.jpg');
[license_height,license_width]=size(Ib);

% I1 = imresize(I0,[162 63]);
% imwrite(I1,'10.jpg','jpg');
Ib_double=double(Ib);  %将uint8型转换为double型，否则不能计算统计量

thresh=graythresh(Ib);%二值化
Ib_d=im2bw(Ib,thresh);
%figure,imshow(Ib_d)
%%
%%problem 1   split not correct , try to get the exact height and width
for i=1:license_width
    sum_bw(i)= sum(Ib_d(:,i));
end
for i=1:license_height
    sum_bw_hei(i)=sum(Ib_d(i,:));
end

for i=1:length(sum_bw_hei)
    if sum_bw_hei(i)>400
        sum_bw_hei(i)=0;
    end
end
      
for i=1:length(sum_bw_hei)
    if(sum_bw_hei(i)>2)
        head=i-1;
        break;
    end
end

for i=1:length(sum_bw_hei)
    if(sum_bw_hei(length(sum_bw_hei)-i)>2)
        tail=length(sum_bw_hei)-i;
        break;
    end
end
%useless 



%hard code
head=10;
tail=55;


line_list=find(sum_bw>2);
line_list_rec = zeros(1,length(line_list));

for i = 2:length(line_list)
    if line_list(i)-line_list(i-1)==1
        line_list(i-1) = 0;      %line_list 不为0的点位跳变的前边缘_index 
    else
        line_list_rec(i) = line_list(i);  %line_list_rec 不为0的点对应跳变的后边缘_index
    end
end

head_index = find(line_list~=0); %找到跳变的前边缘
right_index = line_list(head_index);     %对应列号
last_index = find(line_list_rec~=0); %找到跳变的后边缘
myi = line_list_rec(last_index);
left_idnex = [1 myi];

for k = 1:length(right_index)
    %figure(k+1);
    %imshow(Ib_d(head:tail,left_idnex(k):right_index(k)));
end

%left(k)字符左边缘   right_index(k)字符右边缘

%head=23;
%tail=168;

figure,imshow(Ib_d(head:tail,:));
%figure,plot(sum_bw_hei)

%useless

%license_height2 = round(license_height*90/120);
license_height2 = tail-head+1;
digit_begin_index = round(license_width*139.5/440)-1;

%license_width2 = round(license_width*45/440);
license_width2 = round(license_height2 / 2 + 3);
    

%%
%求图片的数据信息，根据论文上的公式进行模糊处理
mean = mean2(Ib_double);
var_deviation=var(Ib_double(:));
std_deviation = std(Ib_double(:))   ;

%定义
template_slide_window=zeros(license_height2,license_width2);
template_distance=zeros(10,license_width - license_width2);

ti_l=zeros(1,10);
tcor_l=zeros(1,10);

for k= 0:35% 35
    template_name = strcat('temp\',int2str(k),'.jpg');
    template = imread(template_name);
    PSF = fspecial('motion',20, 30);
    template = imfilter(template, PSF, 'circular','conv');
    %imcomplement 取反
    template = imcomplement(template);
    template_resized = imresize(template,[license_height2 license_width2]);
    template_double = double(template_resized);
    mean_template = mean2(template_double);
    template_var_deivation = var(template_double(:));
    template_std_deviation = std(template_double(:));
    mean_template_mat=mean_template*ones(license_height2,license_width2);
    for i =1 : license_height2
        for j = 1 : license_width2
            template_norm(i,j) = ( template_double(i,j) - mean_template )* std_deviation / template_std_deviation + mean;
        end
    end
%    resized_name=strcat('template_norm1\',int2str(k),'.jpg');
%    imwrite(uint8(template_norm),resized_name,'jpg') ;

    if k == 15
        figure;
        imshow(template_norm,[]);
    end
    %todo 论文公式2
    
    ti=100;
    tcor=0;
    for i = 1 : license_width - license_width2 + 1
        %Tc
        template_slide_window = template_norm;

        temp1 = round((license_height-license_height2)/2)+1;
        temp2 = round((license_height+license_height2)/2);
        target_slide_window = Ib_double(head:tail,i:i+license_width2-1);
        %It 
         
        
        % corr
            %cor_res=ones(license_height2,license_width2);
            %for i =1 : license_height2
                %cor_res(i,:)=xcorr(template_double(i,:),mean_template_mat(i,:));
            %end

            %figure,plot(template_norm);

            cor_res(:,:,k+1)=corrcoef(template_slide_window,target_slide_window);
            %cor_res(:,:,k+1)
            if ti>sum(sum(abs(template_slide_window - target_slide_window)))/(license_width2*license_height2)
                ti=sum(sum(abs(template_slide_window - target_slide_window)))/(license_width2*license_height2);
            end
            tmp=abs(sum(sum(cor_res(:,:,k+1)))-2);
            if tmp>tcor
               tcor=tmp;
            end
            
        % corr
        
        
        
        
        template_distance(k+1,i)=10*(3-tmp);
        %template_distance(k+1,i) = sum(sum(abs(template_slide_window - target_slide_window)))/(license_width2*license_height2);

    end    
    ti_l(k+1)=ti;
    tcor_l(k+1)=3-tcor;
    
end

%匹配

interval=round(57*license_width/440);
d_min=zeros(1,license_width-license_width2+1);
d_min_tag=ones(1,license_width-license_width2+interval)*15;

for i=1:license_width-license_width2+1
     %公式三
    tmp_template=template_distance(:,i); 
    %tmp_template
    [d_min(i),d_min_tag(i)] = min(tmp_template,[],1); 
    %problem2 threshold of d_min
    %if(d_min(i) > 18)
    %    d_min_tag(i) = 15;
   % end
    
    
    
end    

figure,plot(d_min);
%%
%find min around license_width2
flag = 0;
cnt = 1;
xcand=[];
for index = 1 : license_width - license_width2 + 1
    flag = 0;
    j = index - round(license_width2/2) + 1 ;
    if j <=0 
        j = 1;
    end
    while( abs(index-j)<license_width2/2)
        if d_min(index)>d_min(j)
            flag = 1;
        end
        j = j + 1;
        if j > license_width - license_width2 + 1
            break;
        end
    end
    if flag == 0 
        xcand(cnt) = index;
        cnt = cnt + 1;
    end
end
xcand
%%
%figure,plot(ti_l);
%figure,plot(tcor_l);
%plot(d_min)
%figure,plot(template_distance(1,:));
%figure,plot(template_distance(2,:));

%d_min 为 fig3(a)的d_min曲线

figure,plot(d_min_tag);

%interval
answer=zeros(1,7);
for i=1:length(xcand)
   answer(i) = d_min_tag(xcand(i)) - 1 ;
end

answer
%figure,plot(answer);

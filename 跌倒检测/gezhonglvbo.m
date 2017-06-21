clear all 
close all
clc
 
%source = mmreader('test.avi');
 source = VideoReader('0fall.avi');%读入视频文件 
%for x = 15:17
%for x = 1:source.NumberOfFrames
for x=1:length(source)%获得总帧数，并依次读取 
        b = read(source,x);
        %b = source(x).cdata;  
        
        %P=imnoise(b,'gaussian',0.02');
        %subplot(2,2,2)
        %imshow(P)
        %bb=rgb2gray(b);
          
        %d=size(b);
        P=double(b); 
        [m,n]=size(P); 
        
      %{ 
       % fprintf('m=%d,n=%d\n', m,n);
        
        k=floor([(m+1)/2 (n+1)/2]);%计算图象中心 
        sigma=20; %滤波器标志值，默认为0.5
        bb=zeros(m,n); 
        for i=1:m      
            for j=1:n  
                bb(i,j) =exp(-((i-k(1))^2+(j-k(2))^2)/(4*sigma))/(4*pi*sigma);  %       b(i,j)=exp(-((i^2+j^2)/(2*sigma^2)));     
            end
        end
        Img1=conv2(P,bb,'same'); %用生成的高斯序列卷积运算，进行高斯滤波 
        d=uint8(Img1); 
        subplot(2,3,2);
        imshow(d);
        title('自编高斯滤波去噪结果');       
     
        A2=fspecial('gaussian',20,20);      %生成高斯序列  
        Y3=filter2(A2,b)/255;              %用生成的高斯序列进行滤波  
        subplot(2,3,3);
        imshow(Y3);
        title('高斯滤波去噪结果');   
        %}
        
        %根据需要调节[n*n]的大小！！
        G=medfilt2(b,[10,10]);
        
        %{
        figure(1),subplot(2,2,1),imshow(uint8(b)),title('original')
        subplot(2,2,2),imshow(uint8(G)),title('中值滤波去噪结果');
        %}
        
       %{
        %自编的中值滤波函数。x是需要滤波的图像,nn是模板大小(即n×n)   
        P1=P;  
        nn=21;
        for i=1:m-nn+1  
            for j=1:n-nn+1  
                c=P(i:i+(nn-1),j:j+(nn-1)); %取出x1中从(i,j)开始的n行n列元素,即模板(n×n的)  
                e=c(1,:);      %是c矩阵的第一行  
                for u=2:nn  
                    e=[e,c(u,:)];     %将c矩阵变为一个行矩阵      
                end  
            mm=median(e);      %mm是中值  
            P1(i+(nn-1)/2,j+(nn-1)/2)=mm;   %将模板各元素的中值赋给模板中心位置的元素  
            end  
        end   
        %未被赋值的元素取原值  
        d2=uint8(P1);
        subplot(2,2,3);
        imshow(d2);
        title('自编中值滤波去噪结果');
        %}
        
        %{
        A=fspecial('average',20); %生成系统预定义的3X3滤波器  
        Y=filter2(A,b)/233; %用生成的滤波器进行滤波,并归一化
        subplot(2,3,6);
        imshow(Y);
        title('均值滤波去噪结果'); 
        %}
       
        %i=imread('tuxiang.jpg');
        %i1=rgb2gray(b);%i1灰度图像
        %level = graythresh(G)%使用最大类间方差法寻找一个合适的阈值
        t=im2bw(G);% t是二值图像，不需要求阈值 
        
        %{
        figure(1),subplot(2,1,1),imshow(uint8(b)),title('original')
        subplot(2,1,2),imshow(t),title('二值化');
        %}

       
        se = strel('rectangle',[5,5]);
        e1 = imerode(t,se);
%h = imshow(e);
%figure(3),imshow(e),title('腐蚀')

        e2 = imdilate(e1,strel('rectangle',[10,10]));
%h = imshow(e);
%figure(4),imshow(e),title('膨胀')

%{


se = strel('rectangle',[7,1]);
e = imerode(e,se);
%h = imshow(e);
figure(5),imshow(e),title('腐蚀')

e = imdilate(e,strel('rectangle',[7,1]));
%h = imshow(e);
figure(6),imshow(e),title('膨胀')

        se = strel('rectangle',[8,1]);
        b2 = imerode(b,se);
%}
        
        figure(1),subplot(2,2,1),imshow(b),title('original')
        subplot(2,2,2),imshow(uint8(G)),title('中值滤波');
        subplot(2,2,3),imshow(e1),title('腐蚀');
        subplot(2,2,4),imshow(e2),title('膨胀');
        
        
        %Mov1(n)  = im2frame(uint8(G),gray);              % put frames into movie        
                
end 

%movie2avi(Mov1,'zhongzhilvbo_output','fps',5);               % save movie as avi 
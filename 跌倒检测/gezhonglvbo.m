clear all 
close all
clc
 
%source = mmreader('test.avi');
 source = VideoReader('0fall.avi');%������Ƶ�ļ� 
%for x = 15:17
%for x = 1:source.NumberOfFrames
for x=1:length(source)%�����֡���������ζ�ȡ 
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
        
        k=floor([(m+1)/2 (n+1)/2]);%����ͼ������ 
        sigma=20; %�˲�����־ֵ��Ĭ��Ϊ0.5
        bb=zeros(m,n); 
        for i=1:m      
            for j=1:n  
                bb(i,j) =exp(-((i-k(1))^2+(j-k(2))^2)/(4*sigma))/(4*pi*sigma);  %       b(i,j)=exp(-((i^2+j^2)/(2*sigma^2)));     
            end
        end
        Img1=conv2(P,bb,'same'); %�����ɵĸ�˹���о�����㣬���и�˹�˲� 
        d=uint8(Img1); 
        subplot(2,3,2);
        imshow(d);
        title('�Ա��˹�˲�ȥ����');       
     
        A2=fspecial('gaussian',20,20);      %���ɸ�˹����  
        Y3=filter2(A2,b)/255;              %�����ɵĸ�˹���н����˲�  
        subplot(2,3,3);
        imshow(Y3);
        title('��˹�˲�ȥ����');   
        %}
        
        %������Ҫ����[n*n]�Ĵ�С����
        G=medfilt2(b,[10,10]);
        
        %{
        figure(1),subplot(2,2,1),imshow(uint8(b)),title('original')
        subplot(2,2,2),imshow(uint8(G)),title('��ֵ�˲�ȥ����');
        %}
        
       %{
        %�Ա����ֵ�˲�������x����Ҫ�˲���ͼ��,nn��ģ���С(��n��n)   
        P1=P;  
        nn=21;
        for i=1:m-nn+1  
            for j=1:n-nn+1  
                c=P(i:i+(nn-1),j:j+(nn-1)); %ȡ��x1�д�(i,j)��ʼ��n��n��Ԫ��,��ģ��(n��n��)  
                e=c(1,:);      %��c����ĵ�һ��  
                for u=2:nn  
                    e=[e,c(u,:)];     %��c�����Ϊһ���о���      
                end  
            mm=median(e);      %mm����ֵ  
            P1(i+(nn-1)/2,j+(nn-1)/2)=mm;   %��ģ���Ԫ�ص���ֵ����ģ������λ�õ�Ԫ��  
            end  
        end   
        %δ����ֵ��Ԫ��ȡԭֵ  
        d2=uint8(P1);
        subplot(2,2,3);
        imshow(d2);
        title('�Ա���ֵ�˲�ȥ����');
        %}
        
        %{
        A=fspecial('average',20); %����ϵͳԤ�����3X3�˲���  
        Y=filter2(A,b)/233; %�����ɵ��˲��������˲�,����һ��
        subplot(2,3,6);
        imshow(Y);
        title('��ֵ�˲�ȥ����'); 
        %}
       
        %i=imread('tuxiang.jpg');
        %i1=rgb2gray(b);%i1�Ҷ�ͼ��
        %level = graythresh(G)%ʹ�������䷽�Ѱ��һ�����ʵ���ֵ
        t=im2bw(G);% t�Ƕ�ֵͼ�񣬲���Ҫ����ֵ 
        
        %{
        figure(1),subplot(2,1,1),imshow(uint8(b)),title('original')
        subplot(2,1,2),imshow(t),title('��ֵ��');
        %}

       
        se = strel('rectangle',[5,5]);
        e1 = imerode(t,se);
%h = imshow(e);
%figure(3),imshow(e),title('��ʴ')

        e2 = imdilate(e1,strel('rectangle',[10,10]));
%h = imshow(e);
%figure(4),imshow(e),title('����')

%{


se = strel('rectangle',[7,1]);
e = imerode(e,se);
%h = imshow(e);
figure(5),imshow(e),title('��ʴ')

e = imdilate(e,strel('rectangle',[7,1]));
%h = imshow(e);
figure(6),imshow(e),title('����')

        se = strel('rectangle',[8,1]);
        b2 = imerode(b,se);
%}
        
        figure(1),subplot(2,2,1),imshow(b),title('original')
        subplot(2,2,2),imshow(uint8(G)),title('��ֵ�˲�');
        subplot(2,2,3),imshow(e1),title('��ʴ');
        subplot(2,2,4),imshow(e2),title('����');
        
        
        %Mov1(n)  = im2frame(uint8(G),gray);              % put frames into movie        
                
end 

%movie2avi(Mov1,'zhongzhilvbo_output','fps',5);               % save movie as avi 
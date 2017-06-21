
function varargout = untitled(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 12-Apr-2016 17:53:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled is made visible.

function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA

[filename, pathname] = uigetfile( '*.avi', '开始');
source= VideoReader([pathname,filename]);
vidFrames = read(source);
numFrames = source.NumberOfFrames;
 

% -----------------------  定义一些与视频帧有关的变量 ----------------------
 fr= read(source,1)                     % 先读入第一帧
 if ndims(fr)==3
    fr_bw = rgb2gray(fr);                     % 将第一帧灰度化
 end
fr_size = size(fr);                              % fr_size：视频画面的尺寸大小
width = fr_size(2);
height = fr_size(1);
fg = zeros(height, width);            % 将矩阵fg零初始,之后作为前景目标的存储矩阵 
bg_bw = zeros(height, width);   % 将矩阵bg_bw零初始化
  sigma=30; %滤波器标志值，默认为0.5

% ---------------------定义一些混合高斯算法中需要的变量 ------------------

C = 3;                                                  % C：单高斯模型的数量（通常为3到5）
M = 3;                                                 % M：背景的数量
D = 2.5;                                               % D：偏差阈值
alpha = 0.01;                                     % alpha：学习率更新的速度（0到1之间）
thresh = 0.45;                                    % thresh：前景阈值（0.25到0.75之间）
sd_init = 6;                                         % sd_init：初始化标准方差
w = zeros(height,width,C);            % w：初始化权值数组
mean = zeros(height,width,C);    % mean：像素均值
sd = zeros(height,width,C);          % sd：像素标准方差
u_diff = zeros(height,width,C);    % u_diff：与像素均值的差
p = alpha/(1/C);                              % p：参数学习率
rank = zeros(1,C);                            % rank：w/sd

% --------------------- 一些模型变量的初始赋值 ---------------------------

pixel_depth = 8;                                            % pixel_depth：像素深度为8位
pixel_range = 2^pixel_depth -1;              % pixel range：像素范围,0-255

for i=1:height
    for j=1:width
        for k=1:C
            
            mean(i,j,k) = rand*pixel_range;     % 随机获取均值mean
            w(i,j,k) = 1/C;                                      % 统一赋值权值w
            sd(i,j,k) = sd_init;                                % 初始化方差sd  
            
        end
    end
end

%------------------------------ 处理每一帧 -----------------------------------
for n = 1:100 
   fr = read(source,n);                                                                          % 向fr中读入每一帧
   if ndims(fr)==3
    fr_bw = rgb2gray(fr); 
   end
    P=double(fr_bw); 
    [mlen,nlen]=size(P);
    k=floor([(mlen+1)/2 (nlen+1)/2]);%计算图象中心 
        bb=zeros(mlen,nlen); 
        for i=1:mlen      
            for j=1:nlen  
                bb(i,j) =exp(-((i-k(1))^2+(j-k(2))^2)/(4*sigma))/(4*pi*sigma);  %       b(i,j)=exp(-((i^2+j^2)/(2*sigma^2)));     
            end
        end
        Img1=conv2(P,bb,'same'); %用生成的高斯序列卷积运算，进行高斯滤波  
        fr_bw=uint8(Img1); 
    for m=1:C
        u_diff(:,:,m) = abs(double(fr_bw) - double(mean(:,:,m)));   % 计算像素差值
    end
     
    % 更新每个像素的背景模型
    for i=1:height
        for j=1:width
            
            match = 0;                                                                                   % 像素是否匹配模型的标记变量
            
            for k=1:C                       
                
                if (abs(u_diff(i,j,k)) <= D*sd(i,j,k))                                       % 判断此像素是否匹配模型                    
                    match = 1;                                                                           % 标记匹配                    
                    % 更新权值,均值,方差
                    w(i,j,k) = (1-alpha)*w(i,j,k) + alpha;
                    p = alpha/w(i,j,k);                  
                    mean(i,j,k) = (1-p)*mean(i,j,k) + p*double(fr_bw(i,j));
                    sd(i,j,k) =   sqrt((1-p)*(sd(i,j,k)^2) + p*((double(fr_bw(i,j)) - mean(i,j,k)))^2);
                else                                    
                    w(i,j,k) = (1-alpha)*w(i,j,k);                                               % 与几个模型均不匹配,则权值减小
                    
                end
            end
            
            w(i,j,:) = w(i,j,:)./sum(w(i,j,:));                                                    % 归一化            
            bg_bw(i,j)=0;
            
            for k=1:C
                bg_bw(i,j) = bg_bw(i,j)+ mean(i,j,k)*w(i,j,k);                   % 更新背景
            end
            
            % 若无匹配模型,创建一个新模型
            if (match == 0)
                [min_w, min_w_index] = min(w(i,j,:));  
                mean(i,j,min_w_index) = double(fr_bw(i,j));
                sd(i,j,min_w_index) = sd_init;              
            end

            rank = w(i,j,:)./sd(i,j,:);                                                              % 计算模型优先级
            rank_ind = [1:1:C];

            % 优先级算法
            for k=2:C               
                for m=1:(k-1)                 
                    if (rank(:,:,k) > rank(:,:,m))                     
                        rank_temp = rank(:,:,m);  
                        rank(:,:,m) = rank(:,:,k);
                        rank(:,:,k) = rank_temp;                        

                        rank_ind_temp = rank_ind(m);  
                        rank_ind(m) = rank_ind(k);
                        rank_ind(k) = rank_ind_temp;    
                    end                    
                end
            end
            
            % 计算前景
            match = 1;
            k = 1;            
            fg(i,j) = 0;
            fg2(i,j) = 0;
            
            while ((match == 1)&&(k<=M))    % 判断是否是前景
                if (w(i,j,rank_ind(k)) >= thresh)                                               
                    if (abs(u_diff(i,j,rank_ind(k))) <= D*sd(i,j,rank_ind(k)))  
                        fg(i,j) = 0;                                 % 非前景目标的话,将此点黑化
                        fg2(i,j) = 0;                               % 同上
                        match = 0;                               % 非前景标记
                    else
                       %fg(i,j) = fr_bw(i,j);                 % 是前景的话,将此点保留
                       %fg(i,j) = fr(i,j); 
                       fg(i,j) = 1;                               
                       fg2(i,j)=255;                             % 是前景的话，将此点变白
                    end
                end
                k = k+1;
            end
        end
    end
    
    se=strel('rectangle',[2 2]);  % 构建结构元素,2*2矩阵
    fg2=imfill(fg2,'holes');        % 填充区域
    fg2=imdilate(fg2,se);          % 膨胀
    fg2=bwareaopen(fg2,1)
    fg2 = bwareaopen(fg2,500)
    %figure(1);imshow(fg2);   
    skin=fg;
    fg(:,:,1)= double(fr(:,:,1)).*skin(:,:,1); 
    fg(:,:,2)= double(fr(:,:,2)).*skin(:,:,1);  
    fg(:,:,3)= double(fr(:,:,3)).*skin(:,:,1); 
    %figure(2),imshow(fr)
    %figure(3),imshow(uint8(bg_bw))
    %figure(4),imshow(uint8(fg))
    
  %imshow(uint8(bg_bw))%显示背景图像三行一列定位第2个 subplot(3,1,2),
    imshow(uint8(fg2)) %显示前景图像  subplot(3,1,3),
    axes(handles.axes7);
     %将图像矩阵转换为视频帧
    Mov1(n)  = im2frame(uint8(fg),gray);% put frames into movie前景帧
    % 显示结果图像
    %figure(1),subplot(2,2,1),imshow(fr)                 % 原始视频帧
    %subplot(2,2,2),imshow(uint8(bg_bw))            % 背景帧
    %subplot(2,2,3),imshow(uint8(fg))                     % 前景帧
    %subplot(2,2,4),imshow(uint8(fg2))                   % 白色前景帧
    % 保存成视频形式
   %Mov1(n)  = im2frame(uint8(fg),gray);            
    %Mov2(n)  = im2frame(uint8(bg_bw),gray);           
     %Mov3(n)  = im2frame(uint8(fg2),gray); 
    imshow(fr)%显示输入图像三行一列定位第一个figure(1),subplot(3,1,1),
    axes(handles.axes6);
I = uint8(fg2);
 
threshold = graythresh(I);
 
bw = im2bw(I,threshold);
 
imshow(bw)
 
% remove all object containing fewer than 30 pixels
 
bw = bwareaopen(bw,30);
 

% fill a gap in the pen's cap
 
se = strel('disk',2);
 
bw = imclose(bw,se);
 

% fill any holes, so that regionprops can be used to estimate
 
% the area enclosed by each of the boundaries
 
bw = imfill(bw,'holes');
 
%imshow(bw)
 
ed=edge(bw);
 
%imshow(ed)
 

 L=bwlabel(bw);

L1 = bwlabel(ed);
 
Pwl=zeros(100,max(L1(:)));
 
     for i=1:max(L(:))
        
 
[y,x]=find(L==i);
 
x0=min(x(:));
 
x1=max(x(:));
 
y0=min(y(:));
 
y1=max(y(:));
if x0 > 0 && y0 > 0 
    rectangle('Position',[x0,y0,x1-x0,y1-y0],'edgeColor','g','LineWidth',1)
end
 

Pwl(i)=(x1-x0)/(y1-y0);
if(Pwl(i)>2)
     set(handles.edit2,'string','异常');
 else
   set(handles.edit2,'string','正常'); 
 end
     end   
end
 
% 保存成avi视频
%movie2avi(Mov1,'mixture_of_gaussians_output','fps',15);              
%movie2avi(Mov2,'mixture_of_gaussians_background','fps',15);           
%movie2avi(Mov3,'mixture_of_gaussians_out','fps',15);
 




 %movie2avi(Mov1,'mixture_of_gaussians_output','fps',30);           % save movie as avi 
%movie2avi(Mov2,'mixture_of_gaussians_background','fps',30);           % save movie as avi 


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%clear;
%figure,imshow(uint8(fg));




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

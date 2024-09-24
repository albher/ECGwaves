function varargout = ECGwave(varargin)
% ECGWAVE M-file for ECGwave.fig
%      ECGWAVE, by itself, creates a new ECGWAVE or raises the existing
%      singleton*.
%
%      H = ECGWAVE returns the handle to a new ECGWAVE or the handle to
%      the existing singleton*.
%
%      ECGWAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECGWAVE.M with the given input arguments.
%
%      ECGWAVE('Property','Value',...) creates a new ECGWAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ECGwave_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ECGwave_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% e_dim the above text to modify the response to help ECGwave

% Last Modified by GUIDE v2.5 11-May-2010 12:16:09

% Begin initialization code - DO NOT E_DIM
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ECGwave_OpeningFcn, ...
                   'gui_OutputFcn',  @ECGwave_OutputFcn, ...
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
% End initialization code - DO NOT E_DIM


% --- Executes just before ECGwave is made visible.
function ECGwave_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ECGwave (see VARARGIN)

% Choose default command line output for ECGwave

handles.ref= [0,0];
if ~isempty(varargin),
   handles.sig= varargin{1}; 
else
   handles.sig=[];
end
ReadSig(varargin);
handles.pref= [];
handles.w=[];
handles.Kwave=[];
handles.clus=[];
% Initial values

if size(varargin,1)>0,
   hold off;
   plot(handles.sig); 
   grid on;
 %  xlim([0,size(handles.sig,1)]);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ECGwave wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ECGwave_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles; %.output;


% --- Executes on button press in B_Onset.
function B_Onset_Callback(hObject, eventdata, handles)
% hObject    handle to B_Onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zoom off
handles.figure1;
hold on
if handles.ref(1)~=0,
  x= handles.ref(1);    
  plot(x, handles.sig(x,:),'>w');    
end
[x,y,b]= ginput(1);
x= floor(x);
handles.ref(1)= x;
plot(x,handles.sig(x,:),'>');

guidata(hObject, handles);

% --- Executes on button press in B_End.
function B_End_Callback(hObject, eventdata, handles)
% hObject    handle to B_End (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zoom off
handles.figure1;
hold on
if handles.ref(2)~=0,
  x= handles.ref(2);    
  plot(x, handles.sig(x,:),'<w');    
end
[x,y,b]= ginput(1);
x= floor(x);
handles.ref(2)= x;
handles.figure1;
hold on
plot(x,handles.sig(x,:),'<');

guidata(hObject, handles);

% --- Executes on button press in B_Find.
function B_Find_Callback(hObject, eventdata, handles)
% hObject    handle to B_Find (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N= floor(get(handles.S_Dim,'Value'));
Thr= get(handles.S_Thr,'Value');
ref= handles.ref;
w= ref(2)-ref(1);
d= phasedist(handles.sig, [w, N], ref(1));

valQRS= get(handles.RB_QRS,'Value');
valP= get(handles.RB_P,'Value');
valT= get(handles.RB_T,'Value');

set(handles.E_Mes, 'String', 'Busy');
pause(0.2);
if valQRS==1,
    py= mean(d)+ Thr*std(d);
    point= maxfilt(d, py, (ref(2)-ref(1)));
    point(1)=[]; point(end)=[];
    handles.pref= point;
    handles.w(1)= w;
elseif valP==1,
    if isempty(handles.pref),
       set(handles.E_Mes, 'String', 'Error: First QRS complex');
       return 
    end
    point= maxpoint(d, handles.pref(:,1), -2*w);
    handles.pref(:,2)= point;
    handles.w(2)= w;
elseif valT==1,
    if isempty(handles.pref),
       set(handles.E_Mes, 'String', 'Error: First QRS complex');
       return 
    end
    point= maxpoint(d, handles.pref(:,1), 2*w);       
    handles.pref(:,3)= point;
    handles.w(3)= w;
else
    set(handles.E_Mes, 'String', 'Error: Mark one wave');
    return
end
set(handles.E_Mes, 'String', 'Done');
figure(handles.figure1);
hold off
plot(handles.sig);
grid on
% xlim([0,size(handles.sig,1)]);
hold on;
plot(point, handles.sig(point,:), '>k');
plot(point+w, handles.sig(point+w,:), '<k');

guidata(hObject, handles);


% --- Executes on button press in RB_P.
function RB_P_Callback(hObject, eventdata, handles)
% hObject    handle to RB_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RB_P
val= get(handles.RB_P,'Value');
if val==1,
   set(handles.RB_QRS,'Value',0);
   set(handles.RB_T,'Value',0);
end
handles.ref= [0,0];
guidata(hObject, handles);


% --- Executes on button press in RB_T.
function RB_T_Callback(hObject, eventdata, handles)
% hObject    handle to RB_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RB_T

val= get(handles.RB_T,'Value');
if val==1,
   set(handles.RB_QRS,'Value',0);
   set(handles.RB_P,'Value',0);
end
handles.ref= [0,0];
guidata(hObject, handles);

% --- Executes on button press in RB_QRS.
function RB_QRS_Callback(hObject, eventdata, handles)
% hObject    handle to RB_QRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RB_QRS

val= get(handles.RB_QRS,'Value');
if val==1,
   set(handles.RB_P,'Value',0);
   set(handles.RB_T,'Value',0);
end
handles.ref= [0,0];
guidata(hObject, handles);

function E_Dim_Callback(hObject, eventdata, handles)
% hObject    handle to E_Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Dim as text
%        str2double(get(hObject,'String')) returns contents of E_Dim as a double


% --- Executes during object creation, after setting all properties.
function E_Dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: E_Dim controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', '5'); 

% --- Executes on slider movement.
function S_Thr_Callback(hObject, eventdata, handles)
% hObject    handle to S_Thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Thr= get(hObject, 'Value');
set(handles.E_Thr, 'String', num2str(Thr));


% --- Executes during object creation, after setting all properties.
function S_Thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S_Thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Value', 2);
set(hObject,'Min', 0)
set(hObject,'Max', 10)



function E_Thr_Callback(hObject, eventdata, handles)
% hObject    handle to E_Thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Thr as text
%        str2double(get(hObject,'String')) returns contents of E_Thr as a double

 

% --- Executes during object creation, after setting all properties.
function E_Thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: E_Dim controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', '2');

% --- Executes on slider movement.
function S_Dim_Callback(hObject, eventdata, handles)
% hObject    handle to S_Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
N= get(hObject, 'Value');
set(handles.E_Dim, 'String', num2str(floor(N)));


% --- Executes during object creation, after setting all properties.
function S_Dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S_Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Value', 5);
set(hObject,'Min', 0);
set(hObject,'Max', 10);

% --- Executes on button press in B_Wav.
function B_Wav_Callback(hObject, eventdata, handles)
% hObject    handle to B_Wav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


nw= size(handles.w);
valQRS= get(handles.RB_QRS,'Value');
valP= get(handles.RB_P,'Value');
valT= get(handles.RB_T,'Value');
Nclus= str2num(get(handles.E_Clu,'String'));
N= str2num(get(handles.E_Dim,'String'));
set(handles.E_Mes, 'String', 'Busy');
pause(0.1);
if valQRS==1,
   if nw<1,
      set(handles.E_Mes, 'String', 'No limits QRS marked');
      return
   else
      [wave, idx]= getwave(handles.sig, ...
          [handles.pref(:,1),handles.pref(:,1)+handles.w(1)]);
      tit= 'QRS cluster waves';
   end
elseif valP==1,
   if nw<2,
      set(handles.E_Mes, 'String', 'No limits P marked');
      return
   else
      [wave, idx]= getwave(handles.sig, ...
         [handles.pref(:,2),handles.pref(:,2)+handles.w(2)]);
      tit= 'P cluster waves';
   end
elseif valT==1,
   if nw<3,
      set(handles.E_Mes, 'String', 'No limits T marked');
      return
   else
      [wave, idx]= getwave(handles.sig, ...
          [handles.pref(:,3),handles.pref(:,3)+handles.w(3)]);
      tit= 'T cluster waves';
   end
else
   set(handles.E_Mes, 'String', 'Error: Mark one wave');
   return
end

[clus, Kwave]= Kclus(wave, idx, N*2, Nclus);
set(handles.E_Mes, 'String', 'Done');
handles.figure1;
hold off
plot(wave);
grid on;
hold on
st= {'*b','*g','*r','*c','dm','db','dg','dr','dc','dm'};

for i= 1:Nclus,
   I= find(clus==i);
   plot(idx(I), wave(idx(I),:), st{i}, 'linew', 5);
end
title(tit, 'FontWeight', 'bold', 'FontSize', 12);

handles.Kwave= Kwave;
handles.clus= clus;
guidata(hObject, handles);


function E_Mes_Callback(hObject, eventdata, handles)
% hObject    handle to E_Mes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Mes as text
%        str2double(get(hObject,'String')) returns contents of E_Mes as a double


% --- Executes during object creation, after setting all properties.
function E_Mes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Mes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: E_Dim controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', 'Choise wave onset and end');




% --- Executes on button press in B_RR.
function B_RR_Callback(hObject, eventdata, handles)
% hObject    handle to B_RR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[N, dim]= size(handles.pref);
if dim==0,
   set(handles.E_Mes, 'String', 'Error: Mark QRS complex.');
   return
end

delay= 2;
RR= handles.pref(2:end,1)- handles.pref(1:end-1,1);
figure;
subplot(2,1,1);
plot(RR, '*', 'linew', 1);
hold on
plot(RR,'k');
grid on; box on;
xlim([0,length(RR)]);
title('RR distance', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Cycle(T)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('RR(T)', 'FontWeight', 'bold', 'FontSize', 10);

subplot(2,1,2);
plot(RR(delay:end), RR(1:end-delay+1), '*', 'linew', 1);
grid on; box on;
tit= ['RR Poicare Map, (\Delta= ', num2str(delay), '). '...
    sprintf('RR= %.1f +/- %.1f.', ...
    mean(RR), std(RR))];
title(tit, 'FontWeight', 'bold', 'FontSize', 12);
xlabel('RR(T)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('RR(T-\Delta)', 'FontWeight', 'bold', 'FontSize', 10);


% --- Executes on button press in B_PQ.
function B_PQ_Callback(hObject, eventdata, handles)
% hObject    handle to B_PQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[N, dim]= size(handles.pref);
if dim<2,
   set(handles.E_Mes, 'String', 'Error: Mark P wave.');
   return
end


delay= 2;
PQ= handles.pref(:,1)- handles.pref(:,2);
figure;
subplot(2,1,1);
plot(PQ, '*', 'linew', 1);
hold on
plot(PQ,'k');
grid on; box on;
xlim([0,length(PQ)]);
title('PQ distance', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Cycle(T)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('PQ(T)', 'FontWeight', 'bold', 'FontSize', 10);

subplot(2,1,2);
plot(PQ(delay:end), PQ(1:end-delay+1), '*', 'linew', 1);
grid on; box on;
tit= ['PQ Poicare Map, (\Delta= ', num2str(delay), '). '...
    sprintf('RR= %.1f +/- %.1f.', ...
    mean(PQ), std(PQ))];
title(tit, 'FontWeight', 'bold', 'FontSize', 12);
xlabel('PQ(T)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('PQ(T-\Delta)', 'FontWeight', 'bold', 'FontSize', 10);



% --- Executes on button press in B_QT.
function B_QT_Callback(hObject, eventdata, handles)
% hObject    handle to B_QT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[N, dim]= size(handles.pref);
if dim<3,
   set(handles.E_Mes, 'String', 'Error: Mark T wave.');
   return
end

delay= 2;
QT= handles.pref(:,3)- handles.pref(:,1);
figure;
subplot(2,1,1);
plot(QT, '*', 'linew', 1);
hold on
plot(QT,'k');
grid on; box on;
xlim([0,length(QT)]);
title('QT distance', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Cycle(T)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('QT(T)', 'FontWeight', 'bold', 'FontSize', 10);

subplot(2,1,2);
plot(QT(delay:end), QT(1:end-delay+1), '*', 'linew', 1);
grid on; box on;
tit= ['QT Poicare Map, (\Delta= ', num2str(delay), '). '...
    sprintf('QT= %.1f +/- %.1f.', ...
    mean(QT), std(QT))];
title(tit, 'FontWeight', 'bold', 'FontSize', 12);
xlabel('QT(T)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('QT(T-\Delta)', 'FontWeight', 'bold', 'FontSize', 10);



function E_Clu_Callback(hObject, eventdata, handles)
% hObject    handle to E_Clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Clu as text
%        str2double(get(hObject,'String')) returns contents of E_Clu as a double



% --- Executes during object creation, after setting all properties.
function E_Clu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%set(handles.E_Clu,'String','2');



% --- Executes on button press in B_Clu.
function B_Clu_Callback(hObject, eventdata, handles)
% hObject    handle to B_Clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N= size(handles.sig,2);
Kwave= handles.Kwave;
clus= handles.clus;
valQRS= get(handles.RB_QRS,'Value');
valP= get(handles.RB_P,'Value');
valT= get(handles.RB_T,'Value');

if isempty(Kwave),
   set(handles.E_Mes, 'String', 'Error: Make wave');
   return   
end

[Nclus, w]= size(Kwave);
w= w/N;

if valQRS==1,
   lab= 'QRS ';  
elseif valP==1,
   lab= 'P ';
elseif valT==1,
   lab= 'T ';
else
   set(handles.E_Mes, 'String', 'Error: Mark one wave');
   return
end

stl= {'b','g','r','c','m','b--','g--','r--','c--','m--'};
stp= {'*b','*g','*r','*c','dm','db','dg','dr','dc','dm'};


k=0;
while N>0,
   figure;
   Ni= min(3,N);
   for j= 1:Ni,
      k= k+1;
      idx= (k-1)*w+1:k*w;
      subplot(Ni+1,1,j);
      for i= 1:Nclus,
         box on; grid on; hold on
         plot(Kwave(i, idx), stl{i}, 'linew', 2);
      end
      axis([1, w, min(min(Kwave(:,idx))), max(max(Kwave(:,idx)))]);
      ylabel(['L-',num2str(k)],'FontWeight', 'bold', 'FontSize', 8);
   end
   xlabel('sample','FontWeight', 'bold', 'FontSize', 8); 
   subplot(Ni+1,1,Ni+1);
   for i= 1:Nclus,
      I= find(clus==i);
      plot(I,clus(I), stp{i}, 'linew', 2);
      hold on; box on;
      xlim([0,length(clus)]);
      xlabel('ECG cycle', 'FontWeight', 'bold', 'FontSize', 8);
      ylabel('Cluster', 'FontWeight', 'bold', 'FontSize', 8);
   end
   plot(clus,'k');
   subplot(Ni+1, 1, 1);
   tit= [lab, ' Clusters. Nº Clus=', num2str(Nclus), ...
       ' . Signal leads and cycles.'];
   title(tit, 'FontWeight', 'bold', 'FontSize', 12);
   N= N-3;
end


% --- Executes on button press in B_Fil.
function B_Fil_Callback(hObject, eventdata, handles)
% hObject    handle to B_Fil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Filter dimension
N=5;

handles.figure1;
hold off
%handles.sig= notch(handles.sig, 0.05, 0.05);
handles.sig= filter(1/N*ones(1,N),1,handles.sig);
handles.sig(1:N,:)=[];
plot(handles.sig);
grid on

guidata(hObject, handles);


% --- Executes on button press in B_Ret.
function B_Ret_Callback(hObject, eventdata, handles)
% hObject    handle to B_Ret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.figure1;
hold off
plot(handles.sig);
grid on
% xlim([0,size(handles.sig,1)]);


% --- Executes on button press in B_Exi.
function B_Exi_Callback(hObject, eventdata, handles)
% hObject    handle to B_Exi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pat.ref= handles.pref;
pat.Kwave= handles.Kwave';
pat.clus= handles.clus;
save data pat

delete(gcf);


function sig= ReadSig(name, N)


return

function [tmax, dmax]= maxfilt(d, ref, int)

t= find(d>ref);
t= t(find((t(2:end)-t(1:end-1))>int));
t= [1; t; length(d)];

nt= length(t)-1;
tmax= zeros(nt,1);
dmax= zeros(nt,1);

for ii= 1:nt,
   [dmax(ii), I]= max(d(t(ii)+1:t(ii+1)));
   tmax(ii)= t(ii)+I;
end

return

function point= maxpoint(dist, ref, Tc)

nref= size(ref(:),1);
ndist= size(dist(:),1);

for ii=1: nref,
   if Tc>0,
      [val,jj]= max(dist(ref(ii):min(ndist,ref(ii)+Tc)));
      point(ii)= ref(ii)+jj-1;
   else
      [val,jj]= max(dist(max(1,ref(ii)+Tc):ref(ii)));
      point(ii)= max(1,ref(ii)+Tc+jj-1);
   end
end
point= point';

return


function [dist]= phasedist(y, win, ref)

% [dist]= phasedist(y, win, ref)
% Input: y signal multilead
%        win(1): window weigth
%        win(2): dimension phase space
%        win(3): reference position, default 0.
%        ref: reference point to measure the distance.
% Output: dist: distance with respect to the reference point.

if nargin<3,
   error('[dist]= phase(y, win, ref)');
end

nw= length(win);
if nw==1,
   win= [win, win, 0];
elseif nw==2,
   win= [win, 0];
end
   
delay= floor(win(1)/win(2));
idx= cumsum([1, delay*ones(1, win(2))]);
idx(end)= idx(end)+ rem(win(1), win(2));
[ny, my]= size(y);

y= [repmat(y(1,:), [win(3), 1]); ...
    y; ...
    repmat(y(end,:), [win(1), 1])];
dist= zeros(ny, 1);

for i= 1:win(2),
   k1= idx(i):idx(i)+ny-1;
   k2= idx(i+1):idx(i+1)+ny-1;
   yref= repmat(y(k2(ref),:)- y(k1(ref),:), [ny,1]);
   dist= dist+ sum((y(k2,:)-y(k1,:)- yref).^2, 2);
end

% Normal Inverse Distance
dist= sqrt(dist);
I= find(dist~=0); J= find(dist==0);
dist(I)= 1./dist(I);
dist(I)= dist(I)/max(dist(I)); dist(J)= 1;

return


function [wave, idx]= getwave(sig, lim)

ncycle= size(lim,1);
nsig= size(sig,2);

wave=[]; idx(1)=1;
for ii=1:ncycle,
   idx(ii+1)= idx(ii)+ lim(ii,2)- lim(ii,1)+1;    
   kk= [];
   for jj= 1:nsig,
      kk(:,jj)= linspace(sig(lim(ii,1),jj), sig(lim(ii,2),jj), lim(ii,2)-lim(ii,1)+1)';
   end
   wave= [wave; sig(lim(ii,1):lim(ii,2),:)-kk];
end
idx(ii+1)=[];

return


function [clus, Kwave]= Kclus(wave, idx, N, Nclus)

[n, dim]= size(wave);
nidx= length(idx);
idx(nidx+1)= n+1;
w= idx(2)-idx(1);
widx= cumsum([floor(w/N/2), floor(w/N)*ones(1, N*dim-1)]);

Twave= zeros(nidx, w*dim);
for i= 1:nidx,
   tmp= wave(idx(i):idx(i+1)-1,:);
   Twave(i,:)= tmp(:)';
end

clus= kmeans(Twave(:,widx), Nclus);

Kwave= zeros(Nclus, w*dim);
for i= 1:Nclus,
   I= find(clus==i);
   Kwave(i,:)= mean(Twave(I,:)); 
end

return


function [sig, B, A]= notch(sig, w0, BW)

if nargin<1,
   error('sig= notch(sig, [w0], [BW])');
elseif nargin==1,
   w0= 0.05; BW= 0.05;
elseif nargin==2,
   BW= w0;
end
    
a1= 2*cos(2*pi*w0)/(1+tan(pi*BW));
a2= (1-tan(pi*BW))/(1+tan(pi*BW));
A=[1,-a1,a2];  B= 1/2*[(1+a2),-2*a1,(1+a2)];
sig= filter(B,A,sig);



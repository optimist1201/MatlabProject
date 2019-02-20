
tic
clear all
% results_CA10716_SiPM_29.2V_short
% results_CA10928_SiPM_29.2V_short
% results_CA10929_SiPM_29.2V_short

% results_FS6_SiPM_29.2V_short
% results_FF01_375_110_SiPM_29.2V_short
% results_SL280_380_SiPM_29.2V_short
% results_SL360_50_SiPM_29.2V_short
% results_UFS1__SiPM_29.2V_short
% results_UFS5__SiPM_29.2V_short
% name='results_FF01_375_110_SiPM_29.2V_short';

% Example(name,1);
name='results_CA10716_SiPM_29.2V_short';
fid = fopen(['D:\Diplom\Угловые характеристики\Результаты_19_12_18\',name,'.txt'],'r');

% fid = fopen('D:\Diplom\results_SiPM_BOOM_MC_W.txt','r');
%Считываем темновой ток
for j=1:5
fgets(fid);
end
fgets(fid,8);
I_SiPM_LED_OFF=fscanf(fid,'%f');
%Считываем все данные, пропуская текстовые вставки
j=1;
tic
while ischar(fgets(fid))
f=fscanf(fid,'%f');
for i=1:length(f)   
    AllData(j)=f(i);
    j=j+1; 
end
end

%Проверяем что в файле нет ошибок записи и выдаем строчку с ошибкой
%Проверка работает только если данных нет, или записано что-то лишнее, если
%какая-то часть данных заменена на нули и тп, то проверка может не
%сработать
    i=0;
    j=0;
    while i==j&&(i*8+1)<length(AllData)
      j=AllData(i*8+1);
      i=i+1;
    end
    if i<length(AllData)/8||~(AllData((i-1)*8+1)==length(AllData)/8)
    X=['Ошибка в строчке:', num2str(i-1)];
    disp(X);
    pause
    end
%Создаем прямоугольную матрицу только со значениями тока (это удобно)
%В столбцах значения для различных вертикальных углов
RecData=reshape(AllData,8,length(AllData)/8);
RecData=RecData';
f=RecData(1,6);
j=1;
i=1;
while f~=100
    f = RecData(j,6);
    O(i)=f;
    k=1;
    while RecData(j,6)==f
        U_f(i,k)=RecData(j,8);
        T(i,k)=RecData(j,3);
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
[S,I]=sort(O);
for i=1:length(I)
   U_f(i,:)= U_f(I(i),:);
end
A=0;
P=0;
% (U_f(i,k-1)+U_f(i,k+1))/2*0.8>U_f(i,k)||
O=size(U_f);
for i=1:O(1)
    y=1;
    for k=3:O(2)-3
        if (U_f(i,k))>(U_f(i,k-1)+U_f(i,k+1))/2*1.2
            U_f(i,k)=(U_f(i,k-1)+U_f(i,k+1))/2;
            P(i,y)=k;
            y=y+1;
        end
        if (U_f(i,k))>(U_f(i,k-2)+U_f(i,k+2))/2*1.5
            U_f(i,k)=(U_f(i,k-2)+U_f(i,k+2))/2;
             P(i,y)=k;
            y=y+1;
        end
    end
end



for i=5
    l=1;
     [k,I]=max(U_f(int16(i/2),:)-I_SiPM_LED_OFF);
     clear g;
    for j=46+(-1)^i*45:135+(-1)^i*45
    if U_f(int16(i/2),j)>I_SiPM_LED_OFF*0.9
        g(l)=j;
        l=l+1;
    end
    end
    Rangel(i)=min(g);
    Ranger(i)=max(g);
%     k=-(min(g)+max(g))/2:;
    %Аппроксимируем данные с помощью синуса для всех вертикальных углов,
%вычитая значения темнового тока
% a=@(x) (x(4)*pi/90*(g-(min(g)+max(g))/2)+x(2)*pi/90);
% F0=@(x)(x(1)*cos(a(x)));
 a=@(x) (pi/90*(g-x(4)));
 a1=@(x) (pi/90*(g-x(4))).^2;
 F3=@(x) x(1)*exp(-(a(x).^2/.2./x(3)^2).^x(5))+I_SiPM_LED_OFF;
% %Эти две формулы для коэффициентов пропускания Френеля не определены для
% %значения угла а(х)=0. Работает похоже только из-за того что фаза не равна точно 0 при итерациях, что легко может случится и код не заработает. 
F  =@(x)(U_f(int16(i/2),g)) - F3(x)  ;
low=  [max(U_f(int16(i/2),:)-I_SiPM_LED_OFF)/2,   0,         0,      0,                    1   ];
d0=   [max(U_f(int16(i/2),:)-I_SiPM_LED_OFF)/1,   0.1,       1.2,   (min(g)+max(g))/2,    1.2 ];
upper=[max(U_f(int16(i/2),:)-I_SiPM_LED_OFF)*2,   4,         8,     (min(g)+max(g)),      1.5 ];
X(:,i) = lsqnonlin(F,d0,low,upper);
A=A+sum((F(X(:,i))).^2);

% if i/2==int16(i/2)
% figure(f)
% plot(-(g-(min(g)+max(g))/2)*2,(U_f(int16(i/2),g) - F3(X(:,i)))./U_f(int16(i/2),g),'r'); hold on
% xlabel('HorizontalAngle');
% ylabel('RelativeError');
% xlim([-90 90]);
% figure(m)
% plot(-a(X(:,i))*180/pi,U_f(int16(i/2),g),'r'); hold on,
% plot(-a(X(:,i))*180/pi, F3(X(:,i)),'g'); hold on
% xlabel('HorizontalAngle');
% ylabel('I_Sipm, uA');
% xlim([-90 90]);
% end
% if i/2~=int16(i/2)
% figure(n)
% plot((g-(min(g)+max(g))/2)*2,(U_f(int16(i/2),g) - F3(X(:,i)))./U_f(int16(i/2),g),'r'); hold on
% xlabel('HorizontalAngle');
% ylabel('RelativeError');
% xlim([-90 90]);
% figure(b)
% plot(a(X(:,5))*180/pi,U_f(int16(5/2),g),'r'); hold on,
% plot(a(X(:,5))*180/pi, F3(X(:,5)),'g'); hold on
% xlabel('HorizontalAngle');
% ylabel('I_Sipm, uA');
% xlim([-90 90]);
% end
% plot(a(X(:,i))*180/pi,(U_f(i,g)-I_SiPM_LED_OFF),'.g'); hold on
% plot(a(X(:,i))*180/pi, (F3(X(:,i))-I_SiPM_LED_OFF),'--b'); hold on
end
% b=figure('Name',[name,'appr_dir']);
% figure(b)
% plot(a(X(:,5))*180/pi,(U_f(5,g)-I_SiPM_LED_OFF)/max(F3(X(:,5))-I_SiPM_LED_OFF),'*r'); hold on
% plot(a(X(:,5))*180/pi, (F3(X(:,5))-I_SiPM_LED_OFF)/max(F3(X(:,5))-I_SiPM_LED_OFF),'m'); hold on

plot(a(X(:,5))*180/pi,(U_f(5,g)-I_SiPM_LED_OFF),'*y'); hold on
plot(a(X(:,5))*180/pi, (F3(X(:,5))-I_SiPM_LED_OFF),'m'); hold on
legend(name,[name , '_Approximation'])
Example()

% saveas(b,[name,'appr.png']);
% legend('only SiPM',['only SiPM' , ' Approximation'],'CA10928',['CA10928' , ' Approximation'])
%Ищем средние значения и стандартное отклонение
% AvA=mean(X(1,:));
% AvPh=mean(X(2,:));
% AvK=mean(X(3,:));
% AvC=mean(X(4,:));
% StdA=std(X(1,:));
% StdPh=std(X(2,:));
% StdK=std(X(3,:));
% StdX=std(X(4,:));
% 
% for i=1:O(1)/2
% T_s1(i)=sum(T(i,:))/181;
% A_s(i)=(X(1,2*i-1)+X(1,2*i))/2-AvA;
% end
% i=1:O(1)/2;
% T_s(i)=(T_s1(i) - sum(T_s1)/18);
% R_c = (sum(T_s(i).*A_s(i)))/sqrt(sum(T_s(i).^2)*sum((A_s(i).^2)))
% Trc = R_c*sqrt(16)/sqrt(1-R_c*R_c)
% figure
% plot(T_s,'.');hold on
% plot(A_s);
% 
%  j=figure('Name',[name,'AmplMax']);
% i=1:O(1)/2;
% plot(S(i)-90,X(1,i*2),S(i)-90,X(1,i*2-1));hold on
% xlabel('VerticalAngle');
% ylabel('AmplitudeMax');
% 
%  e=figure('Name',[name,'Sigma']);
% i=1:O(1)/2;
% plot(S(i)-90,1./X(3,i*2),S(i)-90,1./X(3,i*2-1));hold on
% xlabel('VerticalAngle');
% ylabel('Sigma');
% 
%  t=figure('Name',[name,'Power']);
% i=1:O(1)/2;
% plot(S(i)-90,1./X(5,i*2),S(i)-90,1./X(5,i*2-1));hold on
% xlabel('VerticalAngle');
% ylabel('Power');

% saveas(f,[name,'Diff_back.png']);
% saveas(n,[name,'Diff_dir.png']);
% saveas(m,[name,'appr_back.png']);

% saveas(j,[name,'AmplMax.png']);
% saveas(e,[name,'Sigma.png']);
% saveas(t,[name,'Power.png']);

toc

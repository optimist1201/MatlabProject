function Example()
tic
clear all
% results_FS6_SiPM_29.2V_short
% results_FF01_375_110_SiPM_29.2V_short
% results_SL280_380_SiPM_29.2V_short
% results_SL360_50_SiPM_29.2V_short
% results_UFS1__SiPM_29.2V_short
% results_UFS5__SiPM_29.2V_short
% results_only_SiPM_29.2V_short
name1='results_only_SiPM_29.2V_short';
fid = fopen(['D:\Diplom\Угловые характеристики\Результаты_19_12_18\',name1,'.txt'],'r');
%Считываем темновой ток
for j=1:5
fgets(fid);
end
fgets(fid,8);
I_SiPM_LED_OFF=fscanf(fid,'%f')
%Считываем все данные, пропуская текстовые вставки
j=1;
while ischar(fgets(fid))
f=fscanf(fid,'%f');
for i=1:length(f)   
    AllData(j)=f(i);
    j=j+1; 
    i=i+1;
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
   T(i,:)= T(I(i),:);
end
% f=figure('Name',[name,'diff_back']);
% n=figure('Name',[name,'diff_dir']);
% m=figure('Name',[name,'appr_back']);
% b=figure('Name',[name,'appr_dir']);

A=0;
U_f1=U_f-I_SiPM_LED_OFF;
for i=3
    l=1;
     [k,I]=max(U_f(int16(i/2),:)-I_SiPM_LED_OFF);
     clear g;
    for j=46+((-1)^i)*45:135+((-1)^i)*45
    if U_f(int16(i/2),j)>I_SiPM_LED_OFF*1.1
        g(l)=j;
        l=l+1;
    end
    end
%     Rangel(i)=min(g);
%     Ranger(i)=max(g);
%     k=min(g):max(g);
    %Аппроксимируем данные с помощью синуса для всех вертикальных углов,
%вычитая значения темнового тока
a=@(x) (pi/90*g-x(2)*pi/90);
F0=@(x)(x(1)*cos(a(x)));
% %Эти две формулы для коэффициентов пропускания Френеля не определены для
% %значения угла а(х)=0. Работает похоже только из-за того что фаза не равна точно 0 при итерациях, что легко может случится и код не заработает. 
F1 =@(x) 1-(sin(a(x)-asin(x(3)*sin(a(x))))./sin(a(x)+asin(x(3)*sin(a(x))))).^2;
F2 =@(x) 1-(tan(a(x)-asin(x(3)*sin(a(x))))./tan(a(x)+asin(x(3)*sin(a(x))))).^2;
F3 =@(x) F0(x).*(F1(x)*0.5+F2(x)*0.5)+I_SiPM_LED_OFF;
F  =@(x)(U_f(int16(i/2),g)) -   F3(x)   ;
low=  [max(U_f(int16(i/2),:)-I_SiPM_LED_OFF)/2,   0,  0.00  ];
d0=   [max(U_f(int16(i/2),:)-I_SiPM_LED_OFF)/1,   (min(g)+max(g))/2,  0.67]; %Фаза это второе значение, не ставить равным 0 
upper=[max(U_f(int16(i/2),:)-I_SiPM_LED_OFF)*2,   (min(g)+max(g)),  1.00  ];
X(:,i) = lsqnonlin(F,d0,low,upper);
% if i/2==int16(i/2)
% figure(f)
% plot(-a(X(:,i))*180/pi,(U_f(int16(i/2),g) - F3(X(:,i)))./U_f(int16(i/2),g),'r'); hold on
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
% plot(a(X(:,i))*180/pi,(U_f(int16(i/2),g) - F3(X(:,i)))./U_f(int16(i/2),g),'r'); hold on
% xlabel('HorizontalAngle');
% ylabel('RelativeError');
% xlim([-90 90]);
% figure(b)
% plot(a(X(:,1))*180/pi,U_f(int16(1/2),g),'*b'); hold on
% plot(a(X(:,1))*180/pi, F3(X(:,1)),'r'); hold on
% xlabel('HorizontalAngle');
% ylabel('I_Sipm, uA');
% xlim([-90 90]);
% end
% U_f1(int16(i/2),g)=U_f(int16(i/2),g) - F3(X(:,i));
A=A+sum((F(X(:,i))).^2);
end

xlabel('HorizontalAngle');
ylabel('I_Sipm, uA');
xlim([-90 90]);
% plot(a(X(:,3))*180/pi,(U_f(3,g)-I_SiPM_LED_OFF)/max(F3(X(:,3))-I_SiPM_LED_OFF),'*c'); hold on
% plot(a(X(:,3))*180/pi, (F3(X(:,3))-I_SiPM_LED_OFF)/max(F3(X(:,3))-I_SiPM_LED_OFF),'g'); hold on

plot(a(X(:,3))*180/pi,(U_f(3,g)-I_SiPM_LED_OFF),'*c'); hold on
plot(a(X(:,3))*180/pi, (F3(X(:,3))-I_SiPM_LED_OFF),'g'); hold on

grid on
grid minor
%Ищем средние значения и стандартное отклонение
AvA=mean(X(1,:));
AvPh=mean(X(2,:));
AvK=mean(X(3,:));
% AvC=mean(X(4,:));
StdA=std(X(1,:));
StdPh=std(X(2,:));
StdK=std(X(3,:));
% StdX=std(X(4,:));

%Коэффициент корреляции между температурой и током
% for i=1:5
% T_s1(i)=sum(T(i,:))/181;
% A_s(i)=(X(1,2*i-1)+X(1,2*i))/2-AvA;
% end
% i=1:5;
% T_s(i)=(T_s1(i) - sum(T_s1)/5);
% R_c = (sum(T_s(i).*A_s(i)))/sqrt(sum(T_s(i).^2)*sum((A_s(i).^2)))
% Trc = R_c*sqrt(16)/sqrt(1-R_c*R_c)
% figure
% plot(T_s,'.');hold on
% plot(A_s,'r*');
% grid on

% for i=1:181
% T_s1(i)=T(5,i);
% A_s(i)=U_f1(5,i);
% end
% i=1:181;
% T_s(i)=(T_s1(i) - mean(T_s1));
% A_s(i)=A_s(i)-mean(A_s);
% sqrt(sum(T_s(i).^2)*sum(A_s(i).^2))
% R_c = (sum(T_s(i).*A_s(i)))/sqrt(sum(T_s(i).^2)*sum((A_s(i).^2)))
% Trc = R_c*sqrt(16)/sqrt(1-R_c*R_c)
% figure
% plot(T_s,'.');hold on
% plot(A_s);

%Строим красивые графики

% j=figure('Name',[name,'AmplMax']);
% i=1:18;
% plot(S(i)-90,X(1,i*2),S(i)-90,X(1,i*2-1));hold on
% xlabel('VerticalAngle');
% ylabel('AmplitudeMax');
%     e=figure('Name',[name,'RefrInd']);
% i=1:18;
% plot(S(i)-90,1./X(3,i*2),S(i)-90,1./X(3,i*2-1));hold on
% xlabel('VerticalAngle');
% ylabel('RefractionIndex');

% saveas(f,[name,'Diff_back.png']);
% saveas(n,[name,'Diff_dir.png']);
% saveas(m,[name,'appr_back.png']);
% saveas(b,[name,'appr_dir.png']);
% saveas(e,[name,'RefrInd.png']);
% saveas(j,[name,'AmplMax.png']);

toc
end
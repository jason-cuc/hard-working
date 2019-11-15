clc;
clear;
syms x y r ;
f = (x+y)^2 + (x+1)^2 + (y+3)^2;
syms x_tmp y_tmp;
f_tmp = (x_tmp+y_tmp)^2 + (x_tmp+1)^2 + (y_tmp+3)^2;

% һ�׵���:
fx = diff(f,x);
fy = diff(f,y);
grad_f1 = [fx,fy]';  % �ݶ�(������)

% ��ͼ:ԭʼ3d����ͼ
x = -20:0.1:20;
y = -15:0.1:15;
[X,Y] = meshgrid(x,y); 
Z = (X+Y).^2 + (X+1).^2 + (Y+3).^2;
figure(1);
mesh(X,Y,Z);
xlabel('������x'); ylabel('������y'); zlabel('�ռ�����z');
hold on;
% ��ͼ:ԭʼ��
x0 = 10; y0 = -1.5;
z0 = (x0+y0)^2 + (x0+1)^2 + (y0+3)^2;
plot3(x0,y0,z0,'r*')
hold on

% ��ʼ��:
acc = 0.0001;     % ����
x = 10; 
y = -1.5;  % ���
k = 0;     % �½�����:�������⣡��������"�µ�"�����жϵģ�
B=eye(2);
theta=0.25;
sigma=0.05;
fprintf('�����ݶ��½���ʼ:\n');

while 1
    g_old=eval(grad_f1);
    dk=-inv(B)*eval(grad_f1);      %���·���
    xy_old=[x;y];
    m=0;
    mk=0;
    while(m<500)                   %Armijo���� 
    x_tmp=x+theta^m*dk(1);
    y_tmp=y+theta^m*dk(2);
    if(eval(f_tmp)<= eval(f)+sigma*theta^m*g_old'*dk)
        mk=m;
        break;
    end
    m=m+1;
    end
    
    x=x+theta^mk*dk(1);
    y=y+theta^mk*dk(2);
    xy_new=[x;y];
    
    acc_tmp = sqrt( (xy_new(1)-xy_old(1))^2 + (xy_new(2) - xy_old(2))^2 );
    
    if (acc_tmp<= acc  )
       S=double(eval(f));
        x=double(x);
        y=double(y);
        fprintf('��ֵ����Ϊ:(%.5f,%.5f,%0.5f)\n',x,y,S)
        fprintf('��������:%d\n',k)
        plot3(x,y,S,'*')
        break;
    end
    
     g_new=eval(grad_f1);  %����B
     Sk=xy_new-xy_old;
     Yk=g_new-g_old;
     if(Yk'*Sk <=0)
         B=B;
     else
         B=B-(B*Sk*Sk'*B)/(Sk'*B*Sk)+(Yk*Yk')/(Yk'*Sk); 
     end
     k=k+1;
end
    
    
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
H=eye(2);
fprintf('�����ݶ��½���ʼ:\n');
while 1
    g_old=eval(grad_f1);
    dk=-H*eval(grad_f1);    %���·���
    xy_old=[x;y];
    
    x_tmp=x+r*dk(1);
    y_tmp=y+r*dk(2);
    if diff(eval(f_tmp)) == 0     %��ȷ��������
        r_result=0;
    else
    r_result = solve(diff(eval(f_tmp)));
    end
    x= x+r_result*dk(1);
    y= y+r_result*dk(2);
    
    xy_new=[x;y];
    acc_tmp = sqrt( (xy_new(1)-xy_old(1))^2 + (xy_new(2) - xy_old(2))^2 );
    if acc_tmp<= acc    
       S=double(eval(f));
        x=double(x);
        y=double(y);
        fprintf('��ֵ����Ϊ:(%.5f,%.5f,%0.5f)\n',x,y,S)
        fprintf('��������:%d\n',k)
        break;
    end
    g_new=eval(grad_f1);
    Sk=xy_new-xy_old;
    Yk=g_new-g_old;
    H=H+(Sk-H*Yk)*(Sk-H*Yk)'/((Sk-H*Yk)'*Yk);   %����H �Գ���1
    k=k+1;
end
    



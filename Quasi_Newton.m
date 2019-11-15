clc;
clear;
syms x y r ;
f = (x+y)^2 + (x+1)^2 + (y+3)^2;
syms x_tmp y_tmp;
f_tmp = (x_tmp+y_tmp)^2 + (x_tmp+1)^2 + (y_tmp+3)^2;

% 一阶导数:
fx = diff(f,x);
fy = diff(f,y);
grad_f1 = [fx,fy]';  % 梯度(列向量)

% 做图:原始3d曲面图
x = -20:0.1:20;
y = -15:0.1:15;
[X,Y] = meshgrid(x,y); 
Z = (X+Y).^2 + (X+1).^2 + (Y+3).^2;
figure(1);
mesh(X,Y,Z);
xlabel('横坐标x'); ylabel('纵坐标y'); zlabel('空间坐标z');
hold on;
% 做图:原始点
x0 = 10; y0 = -1.5;
z0 = (x0+y0)^2 + (x0+1)^2 + (y0+3)^2;
plot3(x0,y0,z0,'r*')
hold on

% 初始化:
acc = 0.0001;     % 精度
x = 10; 
y = -1.5;  % 起点
k = 0;     % 下降次数:这里特殊！共轭是用"新点"来做判断的！
H=eye(2);
fprintf('共轭梯度下降开始:\n');
while 1
    g_old=eval(grad_f1);
    dk=-H*eval(grad_f1);    %更新方向
    xy_old=[x;y];
    
    x_tmp=x+r*dk(1);
    y_tmp=y+r*dk(2);
    if diff(eval(f_tmp)) == 0     %精确搜索步长
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
        fprintf('极值坐标为:(%.5f,%.5f,%0.5f)\n',x,y,S)
        fprintf('迭代次数:%d\n',k)
        break;
    end
    g_new=eval(grad_f1);
    Sk=xy_new-xy_old;
    Yk=g_new-g_old;
    H=H+(Sk-H*Yk)*(Sk-H*Yk)'/((Sk-H*Yk)'*Yk);   %更新H 对称秩1
    k=k+1;
end
    




%% No-pill trials
 % t=0:0.5:8;
 % x = [97 145 193	206	227	228	222	231	227	215	199	183	158	148	136	120 115;]
  t=0:0.5:4;
  x = [92 129	136	116	115	109	96	93	90;]
%x = [85 101 106 131 141 135 139 133 117 114 108 110 93;]
% x=[85 101 106 131 141 135 139 133 117 114 108 110 93;
%     92 105 136 184 211 187 167 153 156 134 117 109 98;
%     86 98 129 124 140 143 124 132 115 130 113 111 112;
%     93 103 107 144 155 152 152 158 162 151 145 125 127;
%     104 105 120 159 157 155 154 137 137 130 124 119 103];
% x=[85 101 106 131 141 135 139 133 117 ;
%     92 105 136 184 211 187 167 153 156 ;
%     86 98 129 124 140 143 124 132 115 ;
%     93 103 107 144 155 152 152 158 162 ;
%     104 105 120 159 157 155 154 137 137 ];
 [f,omega,beta,s_y,cv,s_yx,r2]=GlucoseModelTemp(t,x)

 

function [f,omega,beta,s_y,cv,s_yx,r2]=GlucoseModelTemp(t,x)
e=0.5*std(x()); % Lấy ra độ lệch chuẩn của 5 dữ liệu
x00 = x 
%x00=mean(x()); % Thực hiện trung bình các mảng dữ liệu
x0=x00-x00(1); % Tạo giá trị trên giá trị glucose cơ sở
a=AUC(t,x0); % Tính ra AUC 
xmax=max(x0); % Trả về giá trị Max trên đường cơ sở
tau=t(find(x0==xmax)); % Trả về thời gian có nồng độ max
omega=pi/t(end); % Tính omega (Đãng nhẽ phải là 2*pi/t)
beta=2*omega/tan(tau*omega); % Tính hằng số Beta
omega0=sqrt(omega^2+(beta/2)^2); % Tính Omega 0
f=xmax*omega0*exp(0.5*beta*tau/omega); % Tính F
p0=[f beta omega omega0];
p1=fminsearch(@p4,p0,[],a); % Bước tối ưu hóa 4
p2=fminsearch(@model4,p1,[],t,x0); % Bước tối ưu hóa 5

tt=0:.1:t(end)+3; % Chia thời gian nhỏ ra làm 10 lần (0.1)
xx=x00(1)+(p2(1)/p2(3)).*exp(-0.5*p2(2).*tt).*sin(p2(3).*tt);; % Thực hiện tính toán với từng đoạn nhỏ thời gian 0.1
%xx= tt
xa=[-1 t(end)+1]; % Giới hạn trục hoành
ya=[x00(1) x00(1)]; % Giới hạn trục tung
figure(1)
% Tạo trục hoành với data trục hoành (tt), data trục tung (xx)
   % Giới hạn trục hoành (xa), giới hạn trục tung (ya)
   % 'k--' kiểu đường
   % t, x00: dữ liệu thực nghiệm | 'ro' kiểu đường
plot(tt,xx,xa,ya,'k--',t,x00,'ro') 
title('Sự thay đổi đường huyết ở bệnh nhân đái tháo đường type 1')
xlabel('Thời gian (giờ)')
ylabel('Nồng độ glucose (mg/dL)')
%plot(t,x00,'ro') 
f=p2(1);
beta=p2(2);
omega=p2(3);
omega0=sqrt(omega^2+(beta/2)^2);

area_under_the_curve=(f/omega0^2)*(1+exp(-0.5*pi*beta/omega))
x_max=(f/omega0)*exp(-0.5*beta*atan(2*omega/beta)/omega)
t_max=atan(2*omega/beta)/omega

n=length(t);
x_mean=mean(x00); % Giá trị trung bình thực nghiệm
st=sum((x00-x_mean).^2); % Độ lệch bình phương
s_y=sqrt(st/(n-1)); % Độ lệch chuẩn
cv=s_y/x_mean; % Độ lệch chuẩn
% Tính dữ liệu với giá trị p ước tính
xp=x00(1)+(p2(1)/p2(3)).*exp(-0.5*p2(2).*t).*sin(p2(3).*t);;
sr=sum((x00-xp).^2);    % Độ lệch bình phương
s_yx=sqrt(sr/(n-3)); % Độ lệch chuẩn
r2=1-sr/st; % Độ lệch chuẩn

figure(2)
plot(tt,xx,xa,ya,'k--'), hold on
title('Sự thay đổi đường huyết ở bệnh nhân đái tháo đường type 1')
xlabel('Thời gian (giờ)')
ylabel('Nồng độ glucose (mg/dL)')
% Vẽ đường thể hiện độ lệch
% e thể hiện độ lệch đầu cuối của 5 đầu dữ liệu vào
errorbar(t,x00,e,'ko'), hold off 
end

function area=AUC(t,x)
x(2:end-1)=2*x(2:end-1);
area=((t(end)-t(1))/(length(t)-1))*sum(x)/2; % Tính diện tích đường cong
end

function y=p4(p,a) % Trả về giá trị chênh lệch giữa AUC thực nghiệm và AUC lý thuyết
a1=(p(1)/(p(4)^2))*(1+exp(-0.5*pi*p(2)/p(3))); % Công thức (5)
y=(a-a1)^2; % [AUCdata - AUC(F, β, ω)]
end

function r=model4(p,t,x0) % Trả về tổng bình phương phần dư
x1=(p(1)/p(3)).*exp(-0.5*p(2).*t).*sin(p(3)*t); % Công thức (3)
r=sum((x0-x1).^2);
end



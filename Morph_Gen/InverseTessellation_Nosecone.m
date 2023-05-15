clear all; close all;
clc
% Input parameters
eps      = 0.05; % The inner circle radius
W0       = 1;
L0       = 1;
t0       = 1; t1  = 1;
N        = 12;  % Number of slices cut

% The given shape (curvature distribution)
syms s
p1 =  26.52; p2 = -116.3; p3 =  162.3; p4 =   -132;
p5 =   63.5; p6 = -47.14; p7 =  18.64; p8 =  41.56;
% f(x) = p1*x.^7 + p2*x.^6 + p3*x.^5 + p4*x.^4 + p5*x.^3 + p6*x.^2 + p7*x + p8;
theta    =       p1*s.^7 + p2*s.^6 + p3*s.^5 + p4*s.^4 + p5*s.^3 + p6*s.^2 + p7*s + p8;
dtheta   =  diff(p1*s.^7 + p2*s.^6 + p3*s.^5 + p4*s.^4 + p5*s.^3 + p6*s.^2 + p7*s + p8,s);
sintheta = @(s) sin((p1*s.^7 + p2*s.^6 + p3*s.^5 + p4*s.^4 + p5*s.^3 + p6*s.^2 + p7*s + p8)/180*pi);
costheta = @(s) cos((p1*s.^7 + p2*s.^6 + p3*s.^5 + p4*s.^4 + p5*s.^3 + p6*s.^2 + p7*s + p8)/180*pi);

% Looking for the max point
theta1 = @(s) -(p1*s.^7 + p2*s.^6 + p3*s.^5 + p4*s.^4 + p5*s.^3 + p6*s.^2 + p7*s + p8);
s_star = fminunc(theta1,0); % theta_max can be found at s_star

xs_int0s = int(costheta,0,s);
delta    = 1 - integral(costheta,0,1);
ys_int0s = int(sintheta,0,s);
y1_int01 = integral(sintheta,0,1);

W_s      = tan(pi/N)*(L0/W0)*(1+eps-xs_int0s-delta); % 2*, half width

s0 = 0;s1 = 1;
dthetas0 = subs(dtheta,s,s0);
dthetas1 = subs(dtheta,s,s1);
W_s0     = vpa(subs(W_s,s,s0));
W_s1     = vpa(subs(W_s,s,s1));
    
% P_tilde  = (W_s0*t0^3*dthetas0 - W_s1*t1^3*dthetas1)/y1_int01;
% y_star   = W_s0*t0^3*dthetas0/P_tilde;
y_star  = integral(sintheta,0,s_star);
P_tilde = W_s0*t0^3*dthetas0/y_star;

% Determine h(s): Thickness of the elastica
T_s      = (P_tilde*(y_star - ys_int0s)/(W_s*dtheta)).^(1/3); 
% Determine E(s): Modulus of the elastica
E_s      = (P_tilde*(y_star - ys_int0s)/(W_s*dtheta)).^(3/3); 

% Plot the calculated results
% s10      = [linspace(0,0.80,50),linspace(0.85,1,2)];
s10      = [linspace(0,1,100)];
theta0   = subs(theta,s,s10);
dtheta0  = subs(dtheta,s,s10);
x0       = subs(xs_int0s,s,s10);
h0       = subs(ys_int0s,s,s10);
w0       = subs(W_s,s,s10);
t0       = subs(T_s,s,s10);
E0       = subs(E_s,s,s10);

s10 = double(s10);
x0  = double(x0);
h0  = double(h0);
w0  = double(w0);
t0  = double(t0);
E0  = double(E0); E0=E0/max(E0);

% theta0   = double(theta0);
% dtheta0  = double(dtheta0);
% 
% figure(10)
% plot(s10,theta0,'*r','linewidth',2)
% xlabel('s');ylabel('\theta(s)');
% title('Width distribution')
% 
% figure(11)
% plot(s10,dtheta0,'*r','linewidth',2)
% xlabel('s');ylabel('d\theta(s)');
% title('Width distribution')

figure(1)
plot(s10,w0,'*r','linewidth',2)
xlabel('s');ylabel('W(s)');
title('Width distribution')
% SaveString = strcat('w_vs_s_eps_',num2str(eps),'_N_',num2str(N),'.txt'); 
% save(fullfile(SaveString),'s10','w0','-ascii','-double');

figure(2)
plot(s10,t0,'b','linewidth',2)
xlabel('s');ylabel('T(s)');
title('Thickness distribution')

figure(3)
plot(s10,E0,'r','linewidth',2)
xlabel('s');ylabel('Modulus(s)');
title('Modulus distribution')
% axis([0 1 0 1])
% SaveString = strcat('E_vs_s_eps_',num2str(eps),'_N_',num2str(N),'.txt'); 
% save(fullfile(SaveString),'s10','E0','-ascii','-double');

figure(4)
plot(x0,h0,'b','linewidth',2)
xlabel('x');ylabel('h');
axis equal;

%% Generate Tesselation Plane
% Export the geometric coordinates
x1 = fliplr(s10-1-eps);  y1 = fliplr(w0);
x2 = fliplr(x1);       y2 = fliplr(-y1);
X0 = [x1,x2]; Y0 = [y1,y2];
% X0 = [x1,x1(end)-0.1,x2(1)-0.1,x2]; Y0 = [y1,y1(end),y2(1),y2];
% X0 = [x1,x2]; Y0 = [y1,y2];
figure(5)
plot(X0,Y0,'b'); axis equal;

if mod(N/2,2) == 0
Theta = atan2(X0,Y0)+0*pi/N;% N/2是偶函数
else
Theta = atan2(X0,Y0)+1*pi/N;% N/2是奇函数
end

R = sqrt(X0.^2+Y0.^2);

XR0 = []; YR0 = [];
for ii = 1:N
    angle  = ii*2*pi/N;
    Theta2 = Theta-angle;
    XR     = R.*cos(Theta2);
    YR     = R.*sin(Theta2);
    XR0    = cat(2,XR0,XR(1:end-1));
    YR0    = cat(2,YR0,YR(1:end-1));
end
figure(6)
plot(XR0,YR0,'-b','linewidth',2);axis equal;
length(XR0)
% SaveString = strcat('XR_vs_YR_eps_',num2str(eps),'_a_b_',num2str(a_b),'_N_',num2str(N),'.txt'); 
% save(fullfile(SaveString),'XR0','YR0','-ascii','-double');

%% Plot 3D shape
h    = (h0);
x    = x0(end)-x0+eps;
L    = max(x);
y    = linspace(-tan(pi/N)*L, tan(pi/N)*L, length(x));
X    = repmat(x, length(y), 1);
Y    = repmat(y', 1, length(x));
Z    = repmat(h, length(y), 1);

Z(find(abs(Y)>(tan(pi/N).*(X)))) = NaN; % this sets the width of the beam

figure(7)
subplot 121
plot(x, h)
xlabel('x')
ylabel('h')
legend('Shape fo the elastica')

subplot 122
plot(s10, w0)
hold on
plot(s10, tan(pi/N).*(s10(end:-1:1)+eps), 'k--')
xlabel('s')
ylabel('w')
legend('Shape of the strip','Shape of ''flat branch''')

% Plot one strip
figure(8)
surf(X, Y, Z)
hold on
axis off
% Plot the other strips by rotating the first one
Theta = atan2(Y,X);
R = sqrt(X.^2+Y.^2);
for ii = 1:N
    angle = ii*2*pi/N;
    Theta2 = Theta+angle;
    X2 = R.*cos(Theta2);
    Y2 = R.*sin(Theta2);
    
    surf(X2, Y2, Z)
%     light('position',[2,2,2],'style','local')
    rotate3d on;
end
title('3D shape of tesselated dome')
shading interp
axis equal

s10_new = s10';
w0_new = w0';
fit2 = fit(s10_new,w0_new,'fourier3');
coeffs2=coeffvalues(fit2)';
coeffs2(end+1) = x(1); 

str2 = append('w_newshape_nosecone');
str2 = strrep(str2,'.','');
writematrix(coeffs2, str2);

Width = 40; 
Length = 40; 
Square_Size = 0.5;  
x_vf = linspace(0,1,Length/Square_Size);
Em = 1;
Ef = 30;
Ei = 0.5*Ef;
E_interp = interp1(s10,E0,x_vf,'spline');


for i=1:length(x_vf)
    if E_interp(i) <= Ei/Ef
        vf(i) = (E_interp(i)-(Em/Ei))/(1-(Em/Ei));
    elseif E_interp(i) > Ei/Ef
        vf(i) = (E_interp(i)-(Ei/Ef))/(1-(Ei/Ef));
    end
end

% for i=1:length(x_vf)
%     
%     vf(i) = (E_interp(i)-(Em/Ef))/(1-(Em/Ef));
% end
vf = vf';

str1 = append('vf_newshape_nosecone_3mat');
str1 = strrep(str1,'.','');

writematrix(vf, str1);
figure(10)
plot(x_vf, vf)

% for i=1:length(x_vf)
%     
%     vf(i) = (E_interp(i)-(Em/Ef))/(1-(Em/Ef));
% end
% vf = vf';
% 
% str1 = append('vf_newshape_nosecone');
% str1 = strrep(str1,'.','');
% 
% writematrix(vf, str1);

% % Here is to plot the axisymmetric shape, to compare with the firt representation - it should be the same for large N
% [XX YY] = meshgrid(linspace(-1.5, 1.5, 1000));
% R = hypot(XX, YY);
% Z3 = interp1(x, h, R);
% figure(9)
% surf(XX, YY, Z3)
% shading interp; axis equal
% title('3D shape of axisymmetric dome')
%%
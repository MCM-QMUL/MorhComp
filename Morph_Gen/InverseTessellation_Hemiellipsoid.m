% Matlab code to generate discretised morphing structure
% This code was developed by Hirak Kansara 
% Corresponding author: Wei Tan (wei.tan@qmul.ac.uk), Liu Mingchao (mingchao.liu@ntu.edu.sg)

clc;
clear all;
% close all;

a = 1 ; b = 1 ; 

eps      = 0.1; % The inner circle radius
W0       = 1;
L0       = 1;
t0       = 1; t1  = 1;
N        = 8;  % Number of slices cut


[s10, x, E0, w0, h0] = inverse(a,b);


%             var1 = s10;
%             var2 = x;
%             var3 = E0;
%             var4 = w0;
% me = mean(E0);
% A= 30;
% for i = 1:length(E0)
%     if abs(E0(i)) > A
%         E0_new(i) = 1;
%     else
%         E0_new(i) = abs(E0(i)/A);
%     end
%     
% end
% E0_new = E0_new';

A = max(E0);

E0_new = E0'/A;

m = mean(E0_new);

s10_new = s10';
mean = trapz(s10_new, E0_new);
x_new = x';
w0_new = w0';
vol = trapz(s10_new*100, w0_new*100)*2*8*2;
figure(12)
plot(s10_new*100, w0_new*100);


fit1 = fit(x_new,E0_new,'gauss8'); % fourier4 gauss8
coeffs = coeffvalues(fit1)';

fit2 = fit(s10_new,w0_new,'fourier3');
coeffs2=coeffvalues(fit2)';
coeffs2(end+1) = x(1); 
coeffs2(end+1) = h0(end); 
str1 = append('E_fit_eps_',num2str(eps),'_a_',num2str(a),'_b_',num2str(b),'_N_',num2str(N));
str1 = strrep(str1,'.','');

str2 = append('w_fit_eps_',num2str(eps),'_a_',num2str(a),'_b_',num2str(b),'_N_',num2str(N)); % 'w_a_b_075' 'w_fit_eps_',num2str(eps),'_a_',num2str(a),'_b_',num2str(b),'_N_',num2str(N)
str2 = strrep(str2,'.','');

writematrix(coeffs, str1);
writematrix(coeffs2, str2);

figure(3)
plot(s10,E0_new,'r','linewidth',2)
xlabel('s');ylabel('Modulus(s)');
title('Modulus distribution')

figure(10)
plot(s10,E0,'r','linewidth',2)
xlabel('s');ylabel('Modulus(s)');
title('Modulus distribution')

Width = 40; 
Length = 40; 
Square_Size = 0.5;  
x_vf = linspace(0,1,Length/Square_Size);

Ef = 40;
Em = 5;
Ei = 0.2*Ef;

E_interp = interp1(s10_new,E0_new,x_vf,'spline');

% assuming 50% intermediate material anad 50% rigid material
% E_new = E_interp - (Ei/Ef);
% vf_r = (E_new -(Em/Ef))/(1 - (Em/Ef)); 
% 
% E_test = (Ei/Ef) + vf_r*(Ef/Ef) + (1-vf_r)*(Em/Ef);
% rigid = vf_r;
% 
% soft = (1-vf_r)*Em/Ef;
% 
% figure(11)
% hold all
% plot(x, E_interp)
% plot(x, rigid)
% plot(x, soft)
% plot(x, E_new)
% patch([x fliplr(x)], [E_new fliplr(rigid)], 'g')
% patch([x fliplr(x)], [zeros(size(x)) fliplr(rigid)], 'b')
% figure(13)
% plot(x, soft)
% 
% str1 = append('vf_hemi_3mat');
% str1 = strrep(str1,'.','');
% 
% writematrix(vf_r', str1);


% for i=1:length(x)
%     vf_r = (E_interp(i) - vf_s*Em)/(2*Ef);
%     vf_i = 0.8 - vf_r;
% end

for i=1:length(x_vf)
    
    vf(i) = (E_interp(i)-(Em/Ef))/(1-(Em/Ef));
end

% for i=1:length(x)
%     if E_interp(i) < Ei/Ef
%         vf(i) = (E_interp(i)-(Em/Ei))/(1-(Em/Ei));
%     elseif E_interp(i) > Ei/Ef
%         vf(i) = (E_interp(i)-(Ei/Ef))/(1-(Ei/Ef));
%     end
% end

vf = vf';

str1 = append('vf_a_b_075');
str1 = strrep(str1,'.','');

writematrix(vf, str1);
figure(10)
plot(x_vf, vf)

function [var1, var2, var3, var4, var5] = inverse(a,b)

    % Input parameters
    eps      = 0.1; % The inner circle radius
    W0       = 1;
    L0       = 1;
    t0       = 1; t1  = 1;
    N        = 8;  % Number of slices cut

    % The given shape (curvature distribution)
    f = ArcLengthofEllipticCurve(a,b);

    syms s
    
    a1 = f(1:1);
    a2 = f(2:2);
    a3 = f(3:3);
    a4 = f(4:4);
    a5 = f(5:5);
    
    theta    =       a1*s.^4+a2*s.^3+a3*s.^2+a4*s+a5;
    dtheta   =  diff(a1*s.^4+a2*s.^3+a3*s.^2+a4*s+a5);
    sintheta = @(s) sin((a1*s.^4+a2*s.^3+a3*s.^2+a4*s+a5)/180*pi);
    costheta = @(s) cos((a1*s.^4+a2*s.^3+a3*s.^2+a4*s+a5)/180*pi);
    
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
        
    P_tilde  = (W_s0*t0^3*dthetas0 - W_s1*t1^3*dthetas1)/y1_int01;
    y_star   = W_s0*t0^3*dthetas0/P_tilde;
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
    E0  = double(E0);
    
    x    = x0(end)-x0+eps;
    
    var1 = s10;
    var2 = x;
    var3 = E0;
    var4 = w0;
    var5 = h0;
    figure(1)
    plot(s10,w0,'*r','linewidth',2)
    xlabel('s');ylabel('W(s)');
    title('Width distribution')
    
    % figure(2)
    % plot(s10,t0,'b','linewidth',2)
    % xlabel('s');ylabel('T(s)');
    % title('Thickness distribution')
    

%     axis([0 1 0 1])
%     SaveString = strcat('E_vs_s_eps_',num2str(eps),'_a_b_',num2str(a_b),'_N_',num2str(N),'.txt'); 
%     save(fullfile(SaveString),'s10','E0','-ascii','-double');
    
    figure(4)
    plot(x0,h0,'b','linewidth',2)
    xlabel('x');ylabel('h');
%     axis equal;
%     SaveString = strcat('w_s_data_eps_',num2str(eps),'_a_b_',num2str(a_b),'_N_8.txt'); 
%     save(fullfile(SaveString),'s10','w0','-ascii','-double');
    
    
    
    %% Generate Tesselation Plane
    % Export the geometric coordinates
%     x1 = fliplr(s10-1-eps);  y1 = fliplr(w0);
%     x2 = fliplr(x1);       y2 = fliplr(-y1);
    % X0 = [x1,x2]; Y0 = [y1,y2];
    % % X0 = [x1,x1(end)-0.1,x2(1)-0.1,x2]; Y0 = [y1,y1(end),y2(1),y2];
%     X0 = [x1,x2]; Y0 = [y1,y2];
%     figure(5)
%     plot(X0,Y0,'b'); axis equal;
    
%     if mod(N/2,2) == 0
%     Theta = atan2(X0,Y0)+0*pi/N;% N/2ÊÇÅ¼º¯Êý
%     else
%     Theta = atan2(X0,Y0)+1*pi/N;% N/2ÊÇÆæº¯Êý
%     end
%     % 
%     R = sqrt(X0.^2+Y0.^2);
%     
%     XR0 = []; YR0 = [];
%     for ii = 1:N
%         angle  = ii*2*pi/N;
%         Theta2 = Theta-angle;
%         XR     = R.*cos(Theta2);
%         YR     = R.*sin(Theta2);
%         XR0    = cat(2,XR0,XR(1:end-1));
%         YR0    = cat(2,YR0,YR(1:end-1));
%     end
%     figure(6)
%     plot(XR0,YR0,'-b','linewidth',2);axis equal;
%     length(XR0)
    % SaveString = strcat('XR_vs_YR_eps_',num2str(eps),'_a_b_',num2str(a_b),'_N_',num2str(N),'.txt'); 
    % save(fullfile(SaveString),'XR0','YR0','-ascii','-double');
    
    %% Plot 3D shape
    h    = (h0);
%     
    L    = max(x);
    y    = linspace(-tan(pi/N)*L, tan(pi/N)*L, length(x));
    X    = repmat(x, length(y), 1);
    Y    = repmat(y', 1, length(x));
    Z    = repmat(h, length(y), 1);
%     
%     Z(find(abs(Y)>(tan(pi/N).*(X)))) = NaN; % this sets the width of the beam
    
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
%     hold on
    
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

end 
        
%         system(['abaqus cae noGUI=Discretised_morph.py']);
        
function f = ArcLengthofEllipticCurve(a,b)
    % sin(30/180*pi)
    k = sqrt(1-(b/a)^2);
    E_fi = @(x) sqrt(1-k^2*sin(x).^2);

    N = 91;
    theta = linspace(0,90,N);
    for i = 1:N
        theta0 = theta(i)/180*pi;
        fi = asin(b*tan(theta0)./sqrt(a^2+b^2*tan(theta0).^2));
        l1(i) = a*integral(E_fi,0,fi);
    end
    l0 = l1/max(l1); theta0 = 90-theta;
    l0_new = l0';
    theta0_new = theta0';
    fitting = fit(l0_new,theta0_new,'poly4');
    f = coeffvalues(fitting);
end    
        
        % % Here is to plot the axisymmetric shape, to compare with the firt representation - it should be the same for large N
        % [XX YY] = meshgrid(linspace(-1.5, 1.5, 1000));
        % R = hypot(XX, YY);
        % Z3 = interp1(x, h, R);
        % figure(9)
        % surf(XX, YY, Z3)
        % shading interp; axis equal
        % title('3D shape of axisymmetric dome')
        %%




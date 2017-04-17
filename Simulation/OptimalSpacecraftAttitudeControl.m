%Boris Polonsky 2017
%The simulation of spacecraft optimal control introduced by Yonmook Park
%Reference:
%Park, Y. (2015). Robust and optimal attitude control of spacecraft with disturbances. International Journal of Systems Science, 46(7), 1222-1233.
%The simulation result is a little bit different from the results in Park's
%article. Still checking possible issues in the code.
%Configuration
timeInterval=0.001;
timeTerminal=20;
%Initialization
G=diag([1,1.5,2]);
J=diag([5,10,15]);
W=eye(3);
alpha=2;
K_omega=diag([5,10,15]);
K_r=diag([5,10,15]);
fprintf('alpha==%d >=spectral norm of inv(W)*G''==%d\n',alpha,norm(inv(W)*G',2))
fprintf('K_omega-G*inv(W)*G''>0\n')
K_omega-G*inv(W)*G'
omega_init=zeros([3,1]);
r_init=ones([3,1])*pi/4;
%Plot
[omega,r,u,xi,L,omegaNorms,xiNorms]=plotAll(omega_init,r_init,J,G,K_omega,K_r,W,timeInterval,timeTerminal);
timeAxis=(0:size(omega,2)-1)*timeInterval;

figure(1)
plot(timeAxis,omega)
xlabel('time,s')
ylabel('Angular velocity,rad/s')
title('Time histories of angular velocities')
legend('\omega_{1}(t)','\omega_{2}(t)','\omega_{3}(t)')

figure(2)
plot(timeAxis,r)
xlabel('time,s')
ylabel('Euler angles,rad')
title('Time histories of Euler angles')
legend('\phi(t)','\theta(t)','\psi(t)')

figure(3)
plot(timeAxis,u)
xlabel('time,s')
ylabel('control torque,Nm')
title('Time histories of control torque')
legend('u_{1}(t)','u_{2}(t)','u_{3}(t)')

figure(4)
plot(timeAxis,xi)
xlabel('time,s')
ylabel('disturbance')
title('Time histories of disturbances')
legend('\xi_{1}(t)','\xi_{2}(t)','\xi_{3}(t)')

figure(5)
plot(timeAxis,alpha*omegaNorms)
hold on
plot(timeAxis,xiNorms)
xlabel('time,s')
ylabel('Euclidean norms')
%The syntax "\lVert\rVert didn't work..."
title(sprintf('||\\xi(t)|| and \\alpha||\\omega(t)||, (\\alpha=%.1f)',alpha))
legend('\alpha||\omega(t)||','||\xi(t)||')

figure(6)
plot(timeAxis,L)
xlabel('time,s')
title('Time histories of the right-hand side trem of the performance index')
function S=updateS(omega)
S=[0,omega(3,1),-omega(2,1);
    -omega(3,1),0,omega(1,1);
    omega(2,1),-omega(1,1),0];
end
function F=updateF(r)
F=[1,sin(r(1,1))*tan(r(2,1)),cos(r(1,1))*tan(r(2,1));
    0,cos(r(1,1)),-sin(r(1,1));
    0,sin(r(1,1))*sec(r(2,1)),cos(r(1,1))*sec(r(2,1))];
end
function omega=updateOmega(omega,u,xi,J,G,S,timeInterval)
%omega=omega+timeInterval*inv(J)*(S*J*omega+u+G*xi)
omega=omega+timeInterval*(J\(S*J*omega+u+G*xi));
end
function r=updateEulerAngles(omega,r,F,timeInterval)
r=r+timeInterval*F*omega;
end
function u=updateU(omega,r,K_omega,K_r,F)
u=-K_omega*omega-F'*K_r*r;
end
function xi=updateXi(omega,W,G)
%xi=inv(W)*G'*omega;
xi=W\(G'*omega);
end
function Q=updateQ(omega,r,u,xi,K_omega,K_r,F,G,W)
Q=diag(K_omega-G*inv(W)*G',K_r*F*inv(K_omega)*F'*K_r);
end
function N=updateN(F,K_omega,K_r)
N=[zeros(3,3),inv(K_omega)*F'*K_r];
end
function L=UpdateL(omega,r,u,xi,K_omega,K_r,W,F,G)
R=inv(K_omega);
L=omega'*(K_omega-G*inv(W)*G')*omega+r'*(K_r*F*R*F'*K_r)*r...
    +2*u'*R*F'*K_r*r...
    +u'*R*u-xi'*W*xi;
end
%The term in the intergral of performance index
function [omega,r,u,xi,L,omegaNorms,xiNorms]=plotAll(omega_init,r_init,J,G,K_omega,K_r,W,timeInterval,timeTerminal)
%caculate epoch
epoch=ceil(timeTerminal/timeInterval);
%preallocate
omega=zeros(3,epoch+1);
r=zeros(3,epoch+1);
u=zeros(3,epoch+1);
xi=zeros(3,epoch+1);
L=zeros(1,epoch+1);
omega(:,1)=omega_init;
omegaNorms=zeros(1,epoch+1);
xiNorms=zeros(1,epoch+1);
r(:,1)=r_init;
F=updateF(r(:,1));
S=updateS(omega(:,1));
u(:,1)=updateU(omega(:,1),r(:,1),K_omega,K_r,F);
xi(:,1)=updateXi(omega(:,1),W,G);
L(1,1)=UpdateL(omega(:,1),r(:,1),u(:,1),xi(:,1),K_omega,K_r,W,F,G);
omegaNorms(1,1)=norm(omega(:,1),2);
xiNorms(1,1)=norm(xi(:,1),2);
for i=2:epoch+1
    omega(:,i)=updateOmega(omega(:,i-1),u(:,i-1),xi(:,i-1),J,G,S,timeInterval);
    r(:,i)=updateEulerAngles(omega(:,i-1),r(:,i-1),F,timeInterval);
    F=updateF(omega(:,i));
    S=updateS(omega(:,i));
    u(:,i)=updateU(omega(:,i),r(:,i),K_omega,K_r,F);
    xi(:,i)=updateXi(omega(:,i),W,G);
    L(1,i)=UpdateL(omega(:,i),r(:,i),u(:,i),xi(:,i),K_omega,K_r,W,F,G);
    omegaNorms(1,i)=norm(omega(:,i),2);
    xiNorms(1,i)=norm(xi(:,i),2);
end
end
    
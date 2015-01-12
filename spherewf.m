%%
%**************************************************
%                                                 *
% Incident, reflected, and scattered wave field   *
%                                                 *
%                Miles Barnhart                   *
%                                                 *
%**************************************************
clear all;clc;format compact;close all;
% ============= MATERIAL PROPERTIES ===============
v=0.3;%poissons 
E=200E9;%Elastic modulus
rho=7.86E3;%Medium density (Stainless Steel)
nu=E/(2*(1+v));%Lame constant (Shear Modulus)
lamb=v*E/((1+v)*(1-2*v));%Lame constant
% =================================================


w=2*pi*10000;%Frequency
beta=sqrt(w^2*rho/(lamb));
alpha=sqrt((w^2*rho)/(lamb+2*nu));
lambda=2*pi/alpha;
rad=0.5*lambda;%Cavity Radius
N = 50;%Iter count

theta = 0:pi/100:2*pi;

raddef = rad:rad/N:10*rad;
[th,Rad] = meshgrid(theta,raddef);
[x2,y] = pol2cart(th,Rad);
h = waitbar(0,'Running...');
% ======= SPHERICAL BESSEL AND HANKEL EQN's ========
for n1=0:N+1
sph_besa1(n1+1) = sqrt(pi./(2*alpha*rad)).*besselj(n1+0.5,alpha*rad);
sph_hana1(n1+1) = sqrt(pi./(2*alpha*rad)).*besselh(n1+0.5,1,alpha*rad);
sph_hanb1(n1+1) = sqrt(pi./(2*beta*rad)).*besselh(n1+0.5,1,beta*rad);

end

for n1=0:N+1
sph_besa(n1+1,:) = sqrt(pi./(2*alpha*raddef)).*besselj(n1+0.5,alpha*raddef);
sph_hana(n1+1,:) = sqrt(pi./(2*alpha*raddef)).*besselh(n1+0.5,1,alpha*raddef);
sph_hanb(n1+1,:) = sqrt(pi./(2*beta*raddef)).*besselh(n1+0.5,1,beta*raddef);
end
% ============= POTENTIAL AMPLITUDES ===============
for n2 = 0:N
    
E3 = -(1i)^n2*(2*n2 + 1)*((n2^2 - n2 - 0.5*beta^2*rad^2)*sph_besa1(n2+1) + 2*alpha*rad*sph_besa1(n2 + 2));
E4 = -(1i)^n2*(2*n2 + 1)*((n2 - 1)*sph_besa1(n2+1) - alpha*rad*sph_besa1(n2 + 2));

E31 = (n2^2 - n2 - 0.5*beta*rad^2)*sph_hana1(n2+1) + 2*alpha*rad*sph_hana1(n2 + 2);
E32 = -n2*(n2 + 1)*((n2 - 1)*sph_hanb1(n2 + 1) - beta*rad*sph_hanb1(n2 + 2));
E41 = (n2 - 1)*sph_hana1(n2 + 1) - alpha*rad*sph_hana1(n2 + 2);
E42 = -(n2^2 - 1 - 0.5*beta^2*rad^2)*sph_hanb1(n2 + 1) - beta*rad*sph_hanb1(n2 + 2);

M = [E31 E32;E41 E42]\[E3;E4];

An(:,n2+1) = M(1);%Amplitudes
Bn(:,n2+1) = M(2);

clear M 
end

% ============= LEGENDRE POLYNOMIALS =============
P=zeros(N,length(theta));
P(1,:) = 1;
P(2,:) = cos(theta);
for n3=2:N
    P(n3+1,:) = ((2*n3-1)/n3*cos(theta).*P(n3,:) - (n3 - 1)/n3.* P(n3-1,:));
end
%%
for n = 1:length(raddef)
    for zz = 1:length(theta)
        t1=0;t2=t1;t3=t1;t4=t1;t5=t1;
        for nn = 0:N    
    
t1=t1+(2*nn+1)*(1i)^nn*sph_besa(nn+1,n)*P(nn+1,zz);%Incident P-wave
t2=t2+An(:,nn+1)*sph_hana(nn+1,n)*P(nn+1,zz);%Reflected P-wave
t3=t3+Bn(:,nn+1)*sph_hanb(nn+1,n)*P(nn+1,zz);%Reflected SV-wave
t4=t4+An(:,nn+1)*sph_hana(nn+1,n)*P(nn+1,zz);%Scattered P-wave
t5=t5+Bn(:,nn+1)*sph_hanb(nn+1,n)*P(nn+1,zz);%Scatterd SV-wave

        end
        phi1(n,zz)=t1;%Incident P-wave
        phi2(n,zz)=t2;%Reflected P-wave
        phis(n,zz)=t4;%Scattered P-wave
        phit(n,zz)=t1+t2+t4;%Total P-wave potential
        psi1(n,zz)=t3;%Reflected SV-wave
        psi2(n,zz)=t5;%Scatterd SV-wave
        psit(n,zz)=t3+t5;%Total SV-wave potential
    end
    waitbar(n/451);
end
close(h);%close waitbar

figure(1),p2=pcolor(x2,y,real(phi1));axis equal;
set(p2,'LineStyle','none');title('Incident P-wave');
figure(2),p8=pcolor(x2,y,real(phis));axis equal
set(p8,'LineStyle','none');title('Scattered P-wave');
figure(3),p6=pcolor(x2,y,real(phit));axis equal
set(p6,'LineStyle','none');title('Total Longitudinal Wave Field')
figure(4),p10=pcolor(x2,y,real(psi2));axis equal
set(p10,'LineStyle','none');title('Scattered SV-wave')
%% 3D spherical plots using open source sphere3d.m
figure(5),
sphere3d(real(phi1),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Incident P-wave');

figure(6),
sphere3d(real(phis),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Scatterd P-wave');

figure(7),
sphere3d(real(phit),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Total P-wave');

figure(8),
sphere3d(real(psit),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Scattered SV-wave');


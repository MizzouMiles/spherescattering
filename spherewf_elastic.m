%%
%**************************************************
%                                                 *
% Incident, reflected, and scattered wave field   *
%                                                 *
%                Miles Barnhart                   *
%                                                 *
%**************************************************
clear all;clc;format compact;close all;
% ========= MATERIAL PROPERTIES(medium) ==========
% Silicon Rubber
v1=0.48;%poissons 
E1=0.1E9;%Elastic modulus
rho1=1.8E3;%Medium density (Stainless Steel)
nu1=E1/(2*(1+v1));%Lame constant (Shear Modulus)
lamb1=v1*E1/((1+v1)*(1-2*v1));%Lame constant
% =================================================
% ========== MATERIAL PROPERTIES(cavity) ==========
% Stainless Steel
v2=0.3;%poissons 
E2=200E9;%Elastic modulus
rho2=7.86E3;%Medium density (Stainless Steel)
nu2=E2/(2*(1+v2));%Lame constant (Shear Modulus)
lamb2=v2*E2/((1+v2)*(1-2*v2));%Lame constant
% =================================================


w=2*pi*10000;%Frequency

beta1=sqrt(w^2*rho1/(lamb1));
alpha1=sqrt((w^2*rho1)/(lamb1+2*nu1));

beta2=sqrt(w^2*rho2/(lamb2));
alpha2=sqrt((w^2*rho2)/(lamb2+2*nu2));
lambda2=2*pi/alpha2;

p = nu2/nu1;



rad=0.5*lambda2;%Cavity Radius (Inclusion)



N = 75;%Iter count

theta = 0:pi/100:2*pi;

raddef = rad:rad/N:10*rad;
[th,Rad] = meshgrid(theta,raddef);
[x2,y] = pol2cart(th,Rad);
h = waitbar(0,'Running...');
% ======= SPHERICAL BESSEL AND HANKEL EQN's ========
for n1=0:N+1
sph_besa1(n1+1) = sqrt(pi./(2*alpha1*rad)).*besselj(n1+0.5,alpha1*rad);
sph_besb1(n1+1) = sqrt(pi./(2*beta1*rad)).*besselj(n1+0.5,beta1*rad);
sph_hana1(n1+1) = sqrt(pi./(2*alpha1*rad)).*besselh(n1+0.5,1,alpha1*rad);
sph_hanb1(n1+1) = sqrt(pi./(2*beta1*rad)).*besselh(n1+0.5,1,beta1*rad);

sph_besa2(n1+1) = sqrt(pi./(2*alpha2*rad)).*besselj(n1+0.5,alpha2*rad);
sph_besb2(n1+1) = sqrt(pi./(2*beta2*rad)).*besselj(n1+0.5,beta2*rad);
sph_hana2(n1+1) = sqrt(pi./(2*alpha2*rad)).*besselh(n1+0.5,1,alpha2*rad);
sph_hanb2(n1+1) = sqrt(pi./(2*beta2*rad)).*besselh(n1+0.5,1,beta2*rad);

end

for n1=0:N+1
sph_besa11(n1+1,:) = sqrt(pi./(2*alpha1*raddef)).*besselj(n1+0.5,alpha1*raddef);
sph_besb11(n1+1,:) = sqrt(pi./(2*beta1*raddef)).*besselj(n1+0.5,beta1*raddef);
sph_hana11(n1+1,:) = sqrt(pi./(2*alpha1*raddef)).*besselh(n1+0.5,1,alpha1*raddef);
sph_hanb11(n1+1,:) = sqrt(pi./(2*beta1*raddef)).*besselh(n1+0.5,1,beta1*raddef);

sph_besa22(n1+1,:) = sqrt(pi./(2*alpha2*raddef)).*besselj(n1+0.5,alpha2*raddef);
sph_besb22(n1+1,:) = sqrt(pi./(2*beta2*raddef)).*besselj(n1+0.5,beta2*raddef);

sph_hana22(n1+1,:) = sqrt(pi./(2*alpha2*raddef)).*besselh(n1+0.5,1,alpha2*raddef);
sph_hanb22(n1+1,:) = sqrt(pi./(2*beta2*raddef)).*besselh(n1+0.5,1,beta2*raddef);
end
% ============= POTENTIAL AMPLITUDES ===============
for n2 = 0:N
E1 = -(1i)^n2*(2*n2 + 1)*(n2*sph_besa1(n2 + 1) - alpha1*rad*sph_besa1(n2+2));
E2 = -(1i)^n2*(2*n2 + 1)*sph_besa1(n2 + 1);
E3 = -(1i)^n2*(2*n2 + 1)*((n2^2 - n2 - 0.5*beta1^2*rad^2)*sph_besa1(n2+1) + 2*alpha1*rad*sph_besa1(n2 + 2));
E4 = -(1i)^n2*(2*n2 + 1)*((n2 - 1)*sph_besa1(n2+1) - alpha1*rad*sph_besa1(n2 + 2));

E11 = n2*sph_hana1(n2+1)-alpha1*rad*sph_hana1(n2+2);
E21 = sph_hana1(n2+1);
E12 = -n2*(n2+1)*sph_hanb1(n2 + 1);
E22 = -(n2 + 1)*sph_hanb1(n2 + 1) + beta1*rad*sph_hanb1(n2 + 2);
E13 = n2*sph_besa2(n2+1)-alpha2*rad*sph_besa2(n2+1);
E23 = sph_besa2(n2+1);
E33 = (n2^2 - n2 - 0.5*beta2^2*rad^2)*sph_besa2(n2+1) + 2*alpha2*rad*sph_besa2(n2+2);
E43 = (n2 - 1)*sph_besa2(n2+1)-alpha2*rad*sph_besa2(n2+2);
E14 = -n2*(n2+1)*sph_besb2(n2+1);
E24 = -(n2+1)*sph_besb2(n2+1)+beta2*rad*sph_besb2(n2+2);
E34 = -n2*(n2+1)*((n2-1)*sph_besb2(n2+1) - beta2*rad*sph_besb2(n2+2));
E44 = -(n2^2 - 1 - 0.5*beta2^2*rad^2)*sph_besb2(n2+1) - beta2*rad*sph_besb2(n2+2);
E31 = (n2^2 - n2 - 0.5*beta1*rad^2)*sph_hana1(n2+1) + 2*alpha1*rad*sph_hana1(n2 + 2);
E32 = -n2*(n2 + 1)*((n2 - 1)*sph_hanb1(n2 + 1) - beta1*rad*sph_hanb1(n2 + 2));
E41 = (n2 - 1)*sph_hana1(n2 + 1) - alpha1*rad*sph_hana1(n2 + 2);
E42 = -(n2^2 - 1 - 0.5*beta1^2*rad^2)*sph_hanb1(n2 + 1) - beta1*rad*sph_hanb1(n2 + 2);

M = [E11 E12 E13 E14;E21 E22 E23 E24;E31 E32 p*E33 p*E34;E41 E42 p*E43 p*E44]\[E1;E2;E3;E4];

An(:,n2+1) = M(1);%Amplitudes
Bn(:,n2+1) = M(2);
Cn(:,n2+1) = M(3);
Dn(:,n2+1) = M(4);

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
        t1=0;t2=t1;t3=t1;t4=t1;t5=t1;t6=t1;t7=t1;
        for nn = 0:N    
    
t1=t1+(2*nn+1)*(1i)^nn*sph_besa11(nn+1,n)*P(nn+1,zz);%Incident P-wave
t4=t4+An(:,nn+1)*sph_hana11(nn+1,n)*P(nn+1,zz);%Scattered P-wave
t5=t5+Bn(:,nn+1)*sph_hanb11(nn+1,n)*P(nn+1,zz);%Scatterd SV-wave

t6 = t6 + Cn(:,nn+1)*sph_besa22(nn+1,n)*P(nn+1,zz);%Refracted P-wave
t7 = t7 + Dn(:,nn+1)*sph_besb22(nn+1,n)*P(nn+1,zz);%Refracted SV-wave

        end
        phi1(n,zz)=t1;%Incident P-wave
        phis(n,zz)=t4;%Scattered P-wave
        phit(n,zz)=t1+t4;%Total P-wave potential
        psi1(n,zz)=t3;%Reflected SV-wave
        psi2(n,zz)=t5;%Scatterd SV-wave
        psit(n,zz)=t3+t5+t7;%Total SV-wave potential
        phif(n,zz)=-t6;%Refracted P-wave
        psif(n,zz)=-t7;%Refracted SV-wave
    end
    waitbar(n/676);
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

figure(5),p1=pcolor(x2,y,real(phif));axis equal
set(p1,'LineStyle','none');title('Refracted P-wave')
figure(6),p5=pcolor(x2,y,real(psif));axis equal
set(p5,'LineStyle','none');title('Refracted SV-wave')

%% 3D spherical plots using open source sphere3d.m
figure(7),
sphere3d(real(phi1),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Incident P-wave');

figure(8),
sphere3d(real(phis),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Scatterd P-wave');

figure(9),
sphere3d(real(phit),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Total P-wave');

figure(10),
sphere3d(real(psit),0,2*pi,-pi/2,pi/2,rad,6,'surf');title('Scattered SV-wave');




%%%%%%%%%%%%%%% a simple script to plot the kernels %%%%%%%%%%%%%%%

% 1. load Frechet kernel
% pro0=load('proc000000_rhop_alpha_beta_kernel.dat');
% pro1=load('proc000001_rhop_alpha_beta_kernel.dat');
% pro2=load('proc000002_rhop_alpha_beta_kernel.dat');
% pro3=load('proc000003_rhop_alpha_beta_kernel.dat');
  % or load Hessian kernels Habc, where Habc= Ha + Hbs + Hbm + Hc, and
  % Hb=Hbm + Hbs. 
pro0=load('proc000000_rhop_alpha_beta_kernel_Habc.dat');
pro1=load('proc000001_rhop_alpha_beta_kernel_Habc.dat');
pro2=load('proc000002_rhop_alpha_beta_kernel_Habc.dat');
pro3=load('proc000003_rhop_alpha_beta_kernel_Habc.dat');


rho_kappa_mu_kernel=cat(1,pro0,pro1,pro2,pro3);
clear pro*;
  
x=rho_kappa_mu_kernel(:,1);
z=rho_kappa_mu_kernel(:,2);
rho_kernel=rho_kappa_mu_kernel(:,3);
kappa_kernel=rho_kappa_mu_kernel(:,4);
mu_kernel=rho_kappa_mu_kernel(:,5); %mu

xmin=min(x);
xmax=max(x);
zmin=min(z);
zmax=max(z);

[X,Z]=meshgrid([xmin:100:xmax],[zmin:100:zmax]); % 100
%rho_kernel_2D=griddata(x,z,rho_kernel,X,Z);
kappa_kernel_2D=griddata(x,z,kappa_kernel,X,Z);
mu_kernel_2D=griddata(x,z,mu_kernel,X,Z);

min_v=-0.5e-8; max_v=0.5e-8;
figure(1);
subplot(1,2,1)
imagesc(X(1,:)/1000,Z(:,1)/1000,kappa_kernel_2D); %kappa_kernel_2D
caxis([min_v max_v]); % traveltime
colorbar;colormap('jet');
set(gca,'fontsize',18)
xlabel('x [km]','FontSize',18)
ylabel('z [km]', 'FontSize',18)
ax = gca;
ax.YDir = 'normal';
h = colorbar;
%ylabel(h, 'K_{\alpha} [s/m^2]') % traveltime
ylabel(h, 'K_{\alpha} [s/m]') % waveform

% text(10,-15,'(a)','Color','black','FontSize',18)
% line([0;800],[-140;-140],'linestyle','--');
% text(60,-120,'interface','Color','black','FontSize',18)


subplot(1,2,2)
imagesc(X(1,:)/1000,Z(:,1)/1000,mu_kernel_2D);
caxis([min_v max_v]); % traveltime
%caxis([-1e-9 1e-9]); 
colorbar;colormap('jet');
set(gca,'fontsize',18)
xlabel('x [km]','FontSize',18)
ylabel('z [km]', 'FontSize',18)
ax = gca;
ax.YDir = 'normal';
h = colorbar;
%ylabel(h, 'K_{\beta} [s/m^2]') % traveltime
ylabel(h, 'K_{\beta} [s/m]') % waveform
% text(10,-15,'(b)','Color','black','FontSize',18)
% line([0;800],[-140;-140],'linestyle','--');
% text(60,-120,'interface','Color','black','FontSize',18)



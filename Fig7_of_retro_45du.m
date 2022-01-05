clc;   clear all;
% close all;

% th_i =   30/180*pi;  % angle of incidence
% th_ind= 30/180*pi;         % angle at which the structure was designed to operate in a specific way (designed angle of incidence)
% th_rd =  -30/180*pi;  % - 2*th_ind; %retroreflection
th_i =   0/180*pi;  % angle of incidence
th_ind= 0/180*pi;         % angle at which the structure was designed to operate in a specific way (designed angle of incidence)
th_rd =  70/180*pi;  % - 2*th_ind; %retroreflection

th_i =  45/180*pi;  % 58.6003 % angle of incidence
th_ind= 45/180*pi;         % angle at which the structure was designed to operate in a specific way (designed angle of incidence)
th_rd =   -45/180*pi;  % - 2*th_ind; %retroreflection

%%
R= 0.99995%*sqrt(cos(th_ind)/cos(th_rd)); % define the retro reflectivity
%     = 0.99995  for Ana Diaz-Rubio$Eq.(2) ;  
%     adding term  *sqrt(cos(th_ind)/cos(th_rd))   for Ana Diaz-Rubio$Eq.(3)

%%
rho0=1.29; % density of the free medium
c0=343; % sound speed of the free medium

freq=2e3;  omega=2*pi*freq;   % designed frequency
lambda=c0/freq;           % designed wavelength
k0=2*pi/lambda;           % wavenumber of the incidence
%
N = 4;      % defines the last order of Fourier harmonic FOR THE IMPEDANCE expansion (should be EVEN)
%
thetain = 0;  % thetain = (0.01:0.4:pi/2.1);
dz = lambda*1;   Eledz = 20;
z = 0: dz/(Eledz-1) : dz;

for iez = 1:1:length(z)
    
    for ee = 1:1:length(thetain)
        
        D = lambda/(2*sin(th_ind));              D = lambda/ abs(sin(th_ind) - sin(th_rd));% satisfy phase 2\pi
        Dxn = D;
        
        for point = 1:1:length(D)
            
            Z0= rho0 * c0;  % wave impedance of free space for normal incidence
            
            %%
            offset = 0.5*Dxn/(2*N+1); % 0.5 - Middle point
            xn=0+offset: Dxn/(2*N+1): Dxn- Dxn/(2*N+1)+offset;  % coordinate along the gradient
            
            %% define surface impedance
            
            P0=1;      % amplitude of the incident wave
            
            P_t=P0*(exp(-1i*k0*sin(th_ind)*xn) + R*exp(-1i*k0*sin(th_rd)*xn +1i*pi*0 )); %  total pressure filed
            v_t=P0/Z0*(exp(-1i*k0*sin(th_ind)*xn) * cos(th_ind) -  R*exp(-1i*k0*sin(th_rd)*xn +1i*pi*0 ) * cos(th_rd)); % total particle velocity filed
            %             P_t=P0*(exp(-1i*k0*sin(th_ind)*xn) + R*0.5*exp(-1i*k0*sin(th_rd)*xn)+ R*0.5*exp(1i*k0*sin(th_rd)*xn)); %  total pressure filed
            %             v_t=P0/Z0*(exp(-1i*k0*sin(th_ind)*xn) * cos(th_ind) -  R*0.5*exp(-1i*k0*sin(th_rd)*xn ) * cos(th_rd)...
            %                                                                                        -  R*0.5*exp( 1i*k0*sin(th_rd)*xn ) * cos(th_rd)); % total particle velocity filed
            Zs = P_t./v_t; % required surface impedance
% % % % % %        Zs = -1i * Z0/cos(th_ind) * cot( (sin(th_ind) -  sin(th_rd)) *k0 *(xn-D/2) /2 );  % conventional design
             Phi_r = (sin(th_ind) -  sin(th_rd)) .*k0 .*xn;
%             Zs =  1i * Z0/cos(th_ind) .* cot( Phi_r ./2 );  % Ana Diaz-Rubio$Eq.(1); V.S. Asadchy2017 Eq(4)
%             Zs =  Z0 * (1+exp(1i*Phi_r))./(cos(th_ind)-cos(th_rd).*exp(1i.*Phi_r));  % Ana Diaz-Rubio$Eq.(2)
            Zs =  Z0/sqrt(cos(th_ind)*cos(th_rd)) .* (sqrt(cos(th_rd))+sqrt(cos(th_ind)).*exp(1i.*Phi_r))./...
                           (sqrt(cos(th_ind))-sqrt(cos(th_rd)).*exp(1i.*Phi_r));  % Ana Diaz-Rubio$Eq.(3)


            P_t=P0*(exp(-1i*k0*sin(th_ind)*xn) + R*exp(-1i*k0*sin(th_rd)*xn+  1i*pi * 0 )); % V.S. Asadchy2017 Eq.(5). total pressure filed
            %
            v_t=P0/Z0*(exp(-1i*k0*sin(th_ind)*xn) * cos(th_ind) -  R*exp(-1i*k0*sin(th_rd)*xn+  1i*pi * 0 ) * cos(th_rd)); % total particle velocity filed
            %             v_t=P0/Z0*(exp(-1i*k0*sin(th_ind)*xn) -  R*exp(-1i*k0*sin(th_rd)*xn) ); % total particle velocity filed
            %
            Zs = P_t./v_t; % required surface impedance


            %   in comsol          -i*Z0/cos(theta_i)*cot(k*(sin(theta_r) - sin(theta_i))*(x-D/2)/2)
            %% xx is the order of the harmonics
            for ii=1:N+1                        % auxiliary variable
                xx(ii)=ii;
            end
            xx= -(-xx+N/2+1);
            
            %% Finding Fourier harmonics of the impedance
            zp=(fftshift(fft(Zs)/(2*N+1)));
            
            %% Composing the impedance based on the determined harmonics, Z is impedance matrix
            
            Zs_Fourier=zeros(N+1,N+1);
            for in1=1:N+1     %row
                for in2=1:N+1 % colummn
                    Zs_Fourier(in1,in2)=zp(N+1+in2-in1);  % int2 and int1 ?: Zs_Fourier(in1,in2)=zp(N+1+in1-in2);
                end
            end
            
            Ya=zeros(N+1,N+1);
            for nn=1:N+1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                krx(nn) =   k0*sin(th_i) + 2*pi/D*(nn-1-N/2);
                krz(nn) = sqrt(k0^2- krx(nn)^2);     realkrz(nn) = real(krz(nn));  %     (krx./k0)*abs(b)
                %
                
                krx(nn) = isreal(krz(nn)) .* krx(nn);          %%%%%%%%%% important !%%%%%%
                krz(nn) = isreal(krz(nn)) .* krz(nn);          %%%%%%%%%% important !%%%%%%
                
                %
                zan(nn,1) =  rho0*omega/krz(nn) ; % Ya is a diagonal matrix with its main entry representing the admittance of each space harmonic
                theta_r(nn) = asin(krx(nn)/k0);
                
            end
            Za = diag(zan);  %%%
            Ya = (Za)^-1;  %%%
            
            I=eye(N+1); % I is the identity matrix
            REF = pinv(Zs_Fourier*Ya + I )*(Zs_Fourier*Ya - I );            % 'pinv' function works better than 'inv'
            
            Pin_m = zeros(N+1,1);
            Pin_m(N/2+1,1) = P0;  % specify here what incident harmonics excite the surface. by default only one propagating wave is incident.
            Pr_m = REF*Pin_m;                          % Prm gives the amplitudes of the reflected harmonics
            %%
            
            RekrzM(point,:) = realkrz/k0;
            theta_rM(point,:) = real(theta_r);
            
        end
        if ee==1
            Akrz = RekrzM;
            Btheta = theta_rM;
        else
            Akrz = horzcat(Akrz,RekrzM);
            Btheta = horzcat(Btheta,theta_rM);
        end
    end
    a = Pin_m;   b = Pr_m;
    
    %     b(N/2-2) = 0;
    %     b(N/2-1) = 0;
    %     b(N/2) = 0;
    % %     b(N/2+1) = 0;
    %     b(N/2+2) = 0;
    %     b(N/2+3) = 0;
    %     b(N/2+4) = 0;
    
    
    Elex = length(xn);
    
    for j = 1: Elex
        p_ref(j) =  ((exp(1i*krz*z(iez)*1)).*a.'*0 + (exp(-1i*krz*z(iez)*1)).*b.') ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j))  ).' ; % Eq.(6.41)
        v_ref(j) =  (( zan.^-1)' .*( (exp(1i*krz*z(iez)*1)).*a.'*0 - (exp(-1i*krz*z(iez)*1)).*b.' ) ) ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j)) ).';

        %% velocity direction
        v_ref_x(j) =  ((Z0* zan.^-1)' .*( (exp(1i*krz*z(iez)*1)).*a.'*0 + (exp(-1i*krz*z(iez)*1)).*b.'   .* krx./k0 ) ) ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j)) ).' /Z0;
        v_ref_y(j) =  ((Z0* zan.^-1)' .*( (exp(1i*krz*z(iez)*1)).*a.'*0 + (exp(-1i*krz*z(iez)*1)).*b.'   .* krz./k0 ) ) ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j)) ).' /Z0;

        %% total field
        p_tot(j) =  ((exp(1i*krz*z(iez)*1)).*a.'*1 + (exp(-1i*krz*z(iez)*1 )).*b.') ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j))  ).' ; % Eq.(6.41)
        v_tot(j) =  ((Z0* zan.^-1)' .*( (exp(1i*krz*z(iez)*1)).*a.'*1 - (exp(-1i*krz*z(iez)*1 )).*b.' ) ) ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j)) ).' /Z0;
        % total x y velocity
        v_tot_x(j) =  ((Z0* zan.^-1)' .*( (exp(1i*krz*z(iez)*1)).*a.'*sin(th_i) + (exp(-1i*krz*z(iez)*1)).*b.'   .* krx./k0 ) ) ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j)) ).' /Z0;
        v_tot_y(j) =  ((Z0* zan.^-1)' .*(  - (exp(1i*krz*z(iez)*1)).*a.'*cos(th_i) + (exp(-1i*krz*z(iez)*1)).*b.'   .* krz./k0 ) ) ...
            *( sqrt(1/1) .*exp(-1i*krx*xn(j)) ).' /Z0;
    end

    %%
    P_ref(iez, 1 :1: Elex )  = p_ref;       P_tot(iez, 1 :1: Elex )  = p_tot;
    V_ref(iez, 1 :1: Elex )  = v_ref;       V_tot(iez, 1 :1: Elex )   = v_tot;
    
    V_ref_x(iez, 1 :1: Elex )  = v_ref_x;         V_ref_y(iez, 1 :1: Elex )  =  v_ref_y;
    V_tot_x(iez, 1 :1: Elex )  = v_tot_x;         V_tot_y(iez, 1 :1: Elex )  =  v_tot_y;
    
    
end

% plot something

% Incident_Inten_tot = 0.5*real(P_t.*conj(v_t*Z0));   
Incident_Inten_tot = 0.5*real(P0*exp(-1i*k0*sin(th_i)*xn).*conj(P0/Z0*(exp(-1i*k0*sin(th_i)*xn) * cos(th_i))))

ReInten_tot = 0.5*real(P_ref.*conj(V_tot));   

ReInten_ref_x =  0.5*real(P_ref.*conj(V_ref_x));        ReInten_ref_y =  0.5*real(P_ref.*conj(V_ref_y));
ReInten_tot_x =  0.5*real(P_tot.*conj(V_tot_x));        ReInten_tot_y  =  0.5*real(P_tot.*conj(V_tot_y));
%%
% angle_Intensity = atan(ReInten_y./ReInten_x)*180/pi;
Real_V_ref = real(V_ref);
ReInten_ref = 0.5*real(R*exp(-1i*k0*sin(th_rd)*xn ).*conj(V_ref));

% [gx, gz] = gradient(Real_V_ref);

% Harmonics that propagating into far field.
iii = 1;
for nn=1:N+1
    if krz(nn)~=0
        harmonics(iii) = nn-1-N/2;
        marker_xx(iii) = nn;
        iii = iii+1;
    end
end

%% subplot
fontsz = 17; linewd = 1.5; markersz = 8; fontszgca = 3;
X = 0 :Dxn/(Elex-1): Dxn;   %z = 0: dz/(Eledz-1) : dz;  % z = z(end:-1:1);
[X, z] = meshgrid(X,  z);


beishu0 = 0.95;
% figure;
factor = dz/Dxn;   %*factor
n1 = 2; % number of rows
n2 = 3; % number of columns
figure('Position',[0 0 250*n2 250*n1*factor]+50);
nw = 0.6; % normalized width
nh = 0.6;  % normalized height
%%
for k1 = 1:n1
    for k2 = 1:n2
        h=subplot(n1, n2, (k1-1)*n2 + k2);
        
        %
        switch k1
            case 1
                switch k2
                    case 1
                        surf( X/D, z/D, real(P_ref), 'LineStyle','none');
                        view(2);  colormap(jet(80));   shading interp;  c2 = caxis;
                        colorbar('Ticks',[-2:1:2]);  
                         caxis([-2 2]);
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:1.4);
                        title('$\mathrm{Re}(p_{\mathrm{r}})$','Interpreter','latex','FontSize',fontsz+3);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
                        ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 2
                        surf( X/D, z/D, real(V_ref), 'LineStyle','none');
                        view(2);  colormap(jet(80));   shading interp;
                        colorbar('Ticks',[-0.002:0.001:0.002]);%colorbar('Ticks',[-10:0.5:10]);   
                        caxis([-0.002 0.002]);
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:1.4);
                        title('$\mathrm{Re}(v_{\mathrm{r}})$','Interpreter','latex','FontSize',fontsz+3);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
                        %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 3
                        surf( X/D, z/D, real(P_tot), 'LineStyle','none');
                        view(2);  colormap(jet(80));   shading interp;
                        colorbar('Ticks',[-10:1:10]);   caxis([-2 2]);   c2 = caxis;
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        title('$\mathrm{Re}(p_{\mathrm{t}})$','Interpreter','latex','FontSize',fontsz+3);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
                        %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 4
                        plot(xn/Dxn,imag(Zs)/Z0,'-o','color','b','LineWidth',2,'MarkerFaceColor','b');
                        xlabel('$x/D$', 'Interpreter','LaTeX','FontSize',fontsz);
                        ylabel('${\rm Imag}(Z_s/Z0)$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        % xlim([-20 20])
                        legend('');   h = legend;   h.Visible = 'off';
                        set(legend,'Interpreter','latex','edgecolor','none','Location','northwest','FontSize',fontsz);
                        grid on
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 5
                        plot( xn/Dxn, - real(vpa(ReInten_tot_y(1,:)/ Incident_Inten_tot(1), 5)),'-d','color','k','LineWidth',2,'MarkerFaceColor','b');hold on;
                        plot( xn/Dxn,real(vpa(ReInten_tot_y(iez,:)/ Incident_Inten_tot(1), 5)),'--o','color','b','LineWidth',2,'MarkerFaceColor','b');hold on;
                        xlabel('$x/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        ylabel('$I_y$/$I_{\mathrm{i}}$', 'Interpreter', 'LaTeX','FontSize',fontsz);
%                                  ylim([-3e-4 3e-4]);  yticks(-3e-4:3e-4/3:3e-4); %colorbar('Ticks',[-90:30:90]);    caxis([-90 90]);
                        %                         ylabel('Ticks',[-4.46 -4.44]);
%                         colorbar('Ticks',[-10:2:10]);               caxis([-4.46 -4.44]);
                        legend('at $y=0$','at $y=d_z$');   h = legend;  % h.Visible = 'off';
                        set(legend,'Interpreter','latex','edgecolor','none','Location','northoutside','FontSize',fontsz);
                        grid on;
                        %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
% % %                         plot( xn/Dxn, real(vpa(ReInten_tot_x(1,:),5)),'-d','color','k','LineWidth',2,'MarkerFaceColor','b');hold on;
% % %                         plot( xn/Dxn,real(vpa(ReInten_tot_x(iez,:),5)),'--o','color','b','LineWidth',2,'MarkerFaceColor','b');hold on;
% % %                         xlabel('$x/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % %                         ylabel('$I_x$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % %                             ylim([-0e-3 1e-3]);  yticks(-0e-3:1e-3/4:1e-3); %colorbar('Ticks',[-90:30:90]);    caxis([-90 90]);
% % %                         %                         ylabel('Ticks',[-4.46 -4.44]);
% % % %                         colorbar('Ticks',[-10:2:10]);               caxis([-4.46 -4.44]);
% % %                         legend('at $y=0$','at $y=d_z$');   h = legend;  % h.Visible = 'off';
% % %                         set(legend,'Interpreter','latex','edgecolor','none','Location','northoutside','FontSize',fontsz);
% % %                         grid on;
% % %                         %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % %                         set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 6
                        plot( z/dz,   real(vpa(ReInten_tot_x(:,1)/Incident_Inten_tot(1),5)),'-d','color','k','LineWidth',2,'MarkerFaceColor','b');hold on;
                        plot( z/dz, real(vpa(ReInten_tot_x(:,Elex- 0)/Incident_Inten_tot(1),5)),'--o','color','b','LineWidth',2,'MarkerFaceColor','b');hold on;
                        xlabel('$y/d_z$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        ylabel('$I_x$/$I_{\mathrm{i}}$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                                    ylim([-2e-3 2e-3]);  yticks(-2e-3:1e-3:2e-3); %colorbar('Ticks',[-90:30:90]);    caxis([-90 90]);
                        %                         ylabel('Ticks',[-4.46 -4.44]);
%                         colorbar('Ticks',[0:0.5:2]*0.001);               
caxis([0 2]*0.001);
                        legend('at $x=0$','at $x= D$');  % h = legend;  % h.Visible = 'off';
                        set(legend,'Interpreter','latex','edgecolor','none','Location','northoutside','FontSize',fontsz);
                        grid on;
                        %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
% % %                         plot( z/dz, real(vpa(ReInten_tot_y(:,1),5)),'-d','color','k','LineWidth',2,'MarkerFaceColor','b');hold on;
% % %                         plot( z/dz,real(vpa(ReInten_tot_y(:,Elex),5)),'--o','color','b','LineWidth',2,'MarkerFaceColor','b');hold on;
% % %                         xlabel('$y/d_z$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % %                         ylabel('$I_y$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % %                                 ylim([-3e-4 3e-4]);  yticks(-3e-4:3e-4/3:3e-4); %colorbar('Ticks',[-90:30:90]);    caxis([-90 90]);
% % %                         %                         ylabel('Ticks',[-4.46 -4.44]);
% % % %                         colorbar('Ticks',[-10:2:10]);               caxis([-4.46 -4.44]);
% % %                         legend('at $x=0$','at $x= D$');  % h = legend;  % h.Visible = 'off';
% % %                         set(legend,'Interpreter','latex','edgecolor','none','Location','northoutside','FontSize',fontsz);
% % %                         grid on;
% % %                         %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % %                         set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                        %%
                        row = 1;   col = row;
                        conservation = - sum(ReInten_tot_x(1+row:iez-row,1+col)) ...
                                              + sum(ReInten_tot_x(1+row:iez-row, Elex-col)) ...
                                              - sum(ReInten_tot_y(1+row,1+col:Elex-col)) ...
                                             + sum(ReInten_tot_y(iez-row,1+col:Elex-col))
                end
                %%
                
                
            case 2
                switch k2
                    case 1
                        surf( X/D, z/D, real(V_tot), 'LineStyle','none');
                        view(2);  colormap(jet(80));   shading interp;
                        colorbar('Ticks',[-20:10:20]*0.0001);   caxis([-20 20]*0.0001);  c3 = caxis;
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        title('$\mathrm{Re}(v_{\mathrm{t}})$','Interpreter','latex','FontSize',fontsz+3);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
                                                ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz); 
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 2
                        contour(X/D, z/D, ReInten_ref,10,'ShowText','off', 'LineWidth',1);
                        %colorbar('Ticks',[-10:0.5:10]);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
%                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        hold on;
                      % step = 3;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  cos(ReInten(1:step:iez,1:step:Elex)).*sign(ReInten(1:step:iez,1:step:Elex)) .* abs(Inten(1:step:iez,1:step:Elex)),  sin(ReInten(1:step:iez,1:step:Elex)) .* abs(Inten(1:step:iez,1:step:Elex)), 'MaxHeadSize', 0.15 ,'AutoScaleFactor', 2 ,'AutoScale','on');
                        step = 3;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  ...
                              ReInten_ref_x(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
                              ReInten_ref_y(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), 'MaxHeadSize', 0.15 ,'AutoScaleFactor', 0.5 ,'AutoScale','on');
                        hold off;
                        t = title('$I_{\mathrm{r}}$','Interpreter','latex','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 3
                        contour(X/D, z/D, ReInten_tot,10,'ShowText','off', 'LineWidth',1);
                        %colorbar('Ticks',[-10:0.5:10]);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
%                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        hold on;
                        %%
                        step = 1;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  ...
                                      ReInten_tot_x(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
                                      ReInten_tot_y(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
                                      'MaxHeadSize', 0.15 ,'AutoScaleFactor', 2 ,'AutoScale','on');
                        %%
                        hold off;
                        t = title('$I_{\mathrm{t}}$','Interpreter','latex','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 4
                        surf( X/D, z/D, real(V_tot), 'LineStyle','none');
                        view(2);  colormap(jet(80));   shading interp;
                        colorbar('Ticks',[-2:1:2]*0.0001);   caxis([-2 2]*0.0001);  c3 = caxis;
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        title('$v_{\mathrm{t}}$','Interpreter','latex','FontSize',fontsz+3);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
                        %                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);  ;
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 5
                        contour(X/D, z/D, Real_V_ref,10,'ShowText','off', 'LineWidth',1);
                        %colorbar('Ticks',[-10:0.5:10]);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
%                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        hold on;
                      % step = 3;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  cos(ReInten(1:step:iez,1:step:Elex)).*sign(ReInten(1:step:iez,1:step:Elex)) .* abs(Inten(1:step:iez,1:step:Elex)),  sin(ReInten(1:step:iez,1:step:Elex)) .* abs(Inten(1:step:iez,1:step:Elex)), 'MaxHeadSize', 0.15 ,'AutoScaleFactor', 2 ,'AutoScale','on');
                        step = 3;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  ...
                              ReInten_ref_x(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
                              ReInten_ref_y(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), 'MaxHeadSize', 0.15 ,'AutoScaleFactor', 2 ,'AutoScale','on');
                        hold off;
                        t = title('$I_{\mathrm{r}}$','Interpreter','latex','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                    case 6
                        contour(X/D, z/D, Real_V_ref,10,'ShowText','off', 'LineWidth',1);
                        %colorbar('Ticks',[-10:0.5:10]);
                        xlabel('$x/D$','Interpreter', 'LaTeX','FontSize',fontsz);
%                         ylabel('$y/D$', 'Interpreter', 'LaTeX','FontSize',fontsz);
                        xlim([0 Dxn/D]); ylim([0 dz/D]);  xticks(0:0.5:10); yticks(0:0.5:10);
                        hold on;
                        %%
                        step = 1;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  ...
                                      ReInten_tot_x(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
                                      ReInten_tot_y(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
                                      'MaxHeadSize', 0.15 ,'AutoScaleFactor', 2 ,'AutoScale','on');
                        %%
                        hold off;
                        t = title('$I_{\mathrm{t}}$','Interpreter','latex','FontSize',fontsz);
                        set(gca,'position',[(1-nw)/n2/2 + (k2-1)/n2+0.03, (1-nh)/n1/2 + 1-k1/n1,  nw/n2, nh/n1]*beishu0);
                end
        end
        
    end
end
% % % % 
% % % % %
% % % % subplot(2, n2, [1 2])
% % % % % hp = get(subplot(2,n2,[1 2]),'Position');
% % % % % h = colorbar('Position', [hp(1)+hp(3)/2  hp(2)*1  0.015  hp(2)+hp(3)*5]);
% % % % 
% % % % yyaxis left;
% % % % % title('Plots with Different y-Scales');
% % % % grid on;       %   grid minor;
% % % % hold on;
% % % % PP = plot(xx,abs(Pr_m),'-o','LineWidth',2,'MarkerFaceColor','b' ); %,'color','b'
% % % % hold on;
% % % % PP2 = plot(xx,abs(Pr_m),'o', 'MarkerIndices',marker_xx, 'MarkerFaceColor','g', 'MarkerSize',10);
% % % % xlabel('Harmonics $n$','Interpreter', 'LaTeX','FontSize',fontsz);
% % % % ylabel('Amplitude ($|A_n|/P_{\rm in}$)', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % % % xlim([-4 4]);  xticks(-10:1:10); yticks(-10:0.2:10);  %ylim([-90 90]);
% % % % legend('');   h = legend;   h.Visible = 'off';
% % % % set(legend,'Interpreter','latex','edgecolor','none','Location','northwest','FontSize',fontsz);
% % % % %
% % % % yyaxis right;
% % % % grid on;       %grid minor;
% % % % hold on;
% % % % PP = plot(xx,theta_r*180/pi,'--s','LineWidth',2,'MarkerFaceColor','r' ); % ,'color','b'
% % % % hold on;
% % % % PP2 = plot(xx,theta_r*180/pi,'s', 'MarkerIndices',marker_xx, 'MarkerFaceColor','g', 'MarkerSize',10);
% % % % %                         xlabel('Harmonics $n$','Interpreter', 'LaTeX','FontSize',fontsz);
% % % % ylabel('Reflected angle $\theta_{\rm r}$', 'Interpreter', 'LaTeX','FontSize',fontsz);
% % % % xlim([-4 4]); ylim([-90 90]); xticks(-10:1:10); yticks(-90:15:90); %colorbar('Ticks',[-90:30:90]);    caxis([-90 90]);
% % % % legend('');   h = legend;   h.Visible = 'off';
% % % % set(legend,'Interpreter','latex','edgecolor','none','Location','northwest','FontSize',fontsz);
% % % % set(gca,'position',[0.13*0.5    0.5838    0.2866*0.9    0.3412]*beishu0);
% % % % hold off

print(gcf, '-djpeg', '-r300', './Fig2_retro_45du20211219.jpg');
print(gcf, '-depsc2', '-r300', './Fig2_retro_45du20211219.eps');

% figure; More Clear
% % % figure('Position',[0 0 2000*Dxn 2000*dz]+50);
% % % step = 1;  hh = quiver(X(1:step:iez,1:step:Elex)/D, z(1:step:iez,1:step:Elex)/D,  ...
% % %     ReInten_tot_x(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
% % %     ReInten_tot_y(1:step:iez,1:step:Elex) .* abs(ReInten_tot(1:step:iez,1:step:Elex)), ...
% % %     'MaxHeadSize', 0.15 ,'AutoScaleFactor', 2 ,'AutoScale','on');
% % % %
% % % % hold off;
% % % t = title('$I_{\mathrm{t}}$','Interpreter','latex','FontSize',fontsz);


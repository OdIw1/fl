unction varargout=RK4IP_localerr(alpha,beta,gamma,steep,t_raman,...
    T0,P0,m_shape,C0,...
    T,nt,L,nz,...
    nt_plot,nz_plot,atol,rtol)
% NLSE - nonlinear Schrodinger equation solver via RK4IP method
% local error step adjustement
% comments here
tTotal=tic;
switch length(beta)
    case 1
        beta2=beta(1);
        beta3=0;
    otherwise
        beta2=beta(1);
        beta3=beta(2);
end

if ( (nz>30)||(nt>30) )
    error('NLSE: too large time or space grid, both nz and nt as powers of 2 must not exceed 30\n');
end

nz=2^nz;                                                            % nz input is power of 2
nt=2^nt;                                                            % nt input is power of 2
if nz_plot>nz
    nz_plot=nz;
end;
if nt_plot>nt
    nt_plot=nt;
end;
if nt_plot==0
    nt_plot=nt;
end

% time and frequency grids
dt=(2*T0*T)/nt;
t=(-nt/2:nt/2-1)*dt;
w=pi/(T*T0)*[(0:nt/2-1) (-nt/2:-1)];

% initial pulse shape, m_shape=0 - secant, 1 - Gauss, > 1 - super Gauss
% u - time domain, U - frequency domain
if  ( m_shape == 0 )
        u = sqrt(P0)*sech(t/T0).*exp(-0.5i*C0*(t/T0).^2);		
else
        u = sqrt(P0)*exp(-0.5*(1+1i*C0).*(t/T0).^(2*m_shape));
end

u0=u;
U0=fftshift(ifft(u)).*(2*T*T0)/sqrt(2*pi);

% initial pulse energy and soliton order
E0=trapz(abs(u).^2)*dt;
N2=sqrt(gamma*P0*T0^2/abs(beta2));
N3=sqrt(gamma*P0*T0^3/abs(beta3));
fprintf('L_D 2nd order is %g\n',T0^2/abs(beta2));
fprintf('L_D 3rd order is %g\n',T0^3/abs(beta3));
fprintf('L_N is %g\n',1/(gamma*P0));
fprintf('Soliton order (2nd disp) is %g\n',N2);
fprintf('Soliton order (3rd disp) is %g\n',N3);

%preparing plot data
dzPlot=L/(nz_plot-1);
tPlotIndices=round(linspace(1,nt,nt_plot));
uPlotData=zeros(nz_plot,nt_plot);
uPlotData(1,:)=u0(tPlotIndices);
UPlotData=zeros(nz_plot,nt_plot);
UPlotData(1,:)=U0(tPlotIndices);
zOutput=1;                                                      % index of the next dzPlot point for output

h=L/nz;
disp1=exp((1i*0.5*beta2*w.^2-alpha/2+1i*w.^3*beta3/6)*h/2);     % dispersion operator, b3 sign is under question
disp2=exp((1i*0.5*beta2*w.^2-alpha/2+1i*w.^3*beta3/6)*h);
%disp2=sqrt(disp1);

% main loop
z=0;
nSteps=0;                                                       % number of successful steps
nRejected=0;                                                    % number of rejected steps
Steps=[];
err_prev=1;
% PI control terms
ae=0.7;
be=0.4;
while z<L
    % one h step
    uI1=fft(disp1.*ifft(u));
    k11=fft(disp1.*ifft(N(u,h)));
    k12=N(uI1+k11/2,h);
    k13=N(uI1+k12/2,h);
    k14=N(fft(disp1.*ifft(uI1+k13)),h);
    u1=fft(disp1.*ifft(uI1+(k11+2*(k12+k13))/6))+k14/6;
    
    %two h/2 steps
    uI1=fft(disp2.*ifft(u));
    k11=fft(disp2.*ifft(N(u,h/2)));
    k12=N(uI1+k11/2,h/2);
    k13=N(uI1+k12/2,h/2);
    k14=N(fft(disp2.*ifft(uI1+k13)),h/2);
    u2=fft(disp2.*ifft(uI1+(k11+2*(k12+k13))/6))+k14/6;
    
    uI2=fft(disp2.*ifft(u2));
    k21=fft(disp2.*ifft(N(u2,h/2)));
    k22=N(uI2+k21/2,h/2);
    k23=N(uI2+k22/2,h/2);
    k24=N(fft(disp2.*ifft(uI2+k23)),h/2);
    u2=fft(disp2.*ifft(uI2+(k21+2*(k22+k23))/6))+k24/6;
    
    %step control
    %errscale=atol+rtol*max(norm(u1),norm(u2));
    %err=norm(u2-u1)/(sqrt(nt*errscale));
    
    errscale=atol+rtol*max(u1,u2);
    %err=norm(abs(u1-u2)./errscale)/sqrt(nt);    
    err=abs(max(abs(u1-u2)./errscale));
    if err>1
        nRejected=nRejected+1;
        h=0.8*h*max(err^(-ae/5)*err_prev^(be/5),1/5);
        disp1=exp((1i*0.5*beta2*w.^2-alpha/2+1i*w.^3*beta3/6)*h/2);     
        disp2=sqrt(disp1);
    else 
        
        % parse plot data
        while ((min(z+h,L)>=zOutput*dzPlot)&&(zOutput*dzPlot>z))
            theta=zOutput*dzPlot-z;
            if theta<h/2
                b1=theta-3/2*theta^2+2/3*theta^3;
                b23=theta^2-2/3*theta^3;
                b4=-(theta^2)/2+2/3*theta^3;

                utemp=fft(disp1.*ifft(uI1+(b1*k11+b23*(k12+k13))))+k14*b4;
            else
                theta=theta-h/2;
                b1=theta-3/2*theta^2+2/3*theta^3;
                b23=theta^2-2/3*theta^3;
                b4=-(theta^2)/2+2/3*theta^3;

                utemp=fft(disp1.*ifft(uI2+(b1*k21+b23*(k22+k23))))+k24*b4;
            end;                
            uPlotData(zOutput+1,:)=utemp(tPlotIndices);
            Utemp=fftshift(ifft(utemp)).*(2*T*T0)/sqrt(2*pi);
            UPlotData(zOutput+1,:)=Utemp(tPlotIndices);

            zOutput=zOutput+1;
        end
        
        % complete the step
        z=z+h;
        nSteps=nSteps+1;
        Steps=[Steps h];
        h=0.8*h*min(10,err^(-ae/5)*err_prev^(be/5));
        disp1=exp((1i*0.5*beta2*w.^2-alpha/2+1i*w.^3*beta3/6)*h/2);     
        disp2=sqrt(disp1);
        err_prev=err;
        u=u2;
    end; 
end

uL=utemp;
UL=Utemp;
EL=trapz(abs(uL).^2)*dt;                                            % final pulse energy
fprintf('Initial pulse energy = %g\n', E0);
fprintf('Final pulse energy = %g\n', EL);
fprintf('Delta energy is %g\n',abs((EL-E0)/E0));
fprintf('Successful steps number is %g\n',nSteps);
fprintf('Rejected steps number is %g\n',nRejected);

% finalizing plot data
uPlotData(nz_plot,:)=uL(tPlotIndices);
UPlotData(nz_plot,:)=UL(tPlotIndices);
uPlotData=abs(uPlotData).^2;
UPlotData=abs(UPlotData).^2;
[tPlot zPlot]=meshgrid(t(tPlotIndices),linspace(0,L,nz_plot));     % mesh for 3D plots

%passing output
varargout{1}=t/T0;
varargout{2}=fftshift(w)/(2*pi)*T0;
varargout{3}=abs(u0).^2;
varargout{4}=abs(U0).^2;
varargout{5}=angle(u0);
varargout{6}=abs(uL).^2;
varargout{7}=abs(UL).^2;
varargout{8}=angle(uL);
varargout{9}=tPlot./T0;
varargout{10}=zPlot;
varargout{11}=uPlotData./P0;
varargout{12}=UPlotData;
varargout{13}=Steps;

fprintf('Total time is %g   ------------------------------------------------------\n',toc(tTotal));

function r=N(a,step)
%nonlinear operator in the interaction picture
    temp2=abs(a).^2;
    if ((steep==0) && (t_raman==0))
        r=1i*step*gamma*temp2.*a;
    elseif (steep==0)
            r=a.*(1i*step*(gamma*temp2-gamma*t_raman*df(temp2,dt)));
        else
            r=a.*(step*(1i*gamma*temp2-...
                (gamma*steep+1i*gamma*t_raman)*df(temp2,dt)-gamma*steep*df(a,dt).*conj(a)));
    end

%     temp2 = a .* conj(a);
%     temp = conj(a) .* df(a, dt);
%         
%     r=a.*(step*gamma*(1i*temp2-steep*temp-2*(steep+1i*t_raman)*real(temp)));


end

end %NLSE

function d=df(a,dt)
% first order derivative by d u_i=(u_i+1 - u_i-1)/2*h +O(h^2)
    d=(a([2:end 1]) - a([end 1:end-1]))/(2*dt);
    d(1) = d(2) + (d(2)-d(3));
    d(end) = d(end-1) + (d(end-1)-d(end-2));
end











 
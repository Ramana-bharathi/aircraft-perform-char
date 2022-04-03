% % propeller engine

%given constant parameters
m=750;
S=12;
b=10;
C_Do=0.036;
C_lmax=2.7;
e=0.87;

% altitude and air density
h=(0:100:20000)';
sig=sigma(h);
rho=1.225.*sig;

%derived constant parameters
W=m*9.81;
AR=(b^2)/S;
K=(pi*e*AR)^-1;
v_stall=(2*W./(S.*rho.*C_lmax)).^0.5;
v_stall_eq=v_stall.*(sig.^0.5);


%Power model
P_sl=100*745.699872; %hp to W conversion
P=P_sl.*(sig.^0.5);

% vector to store velocity
v_P=zeros(length(h),2);

figure;
hold on;
grid on;
%solving the biquadratic equation and plotting thrust plots
syms v;
for i=1:1:length(h)
    v_=solve(0.25*(rho(i)^2)*S*C_Do*v^4-0.5*rho(i)*P(i)*v+(K*W^2)/S == 0,v,'Real',true);
    if isempty(v_)
        break;
    end
    v_P(i,:)=v_;
    
    if mod(h(i),3000)==0
        v_space=min(v_):0.025*(max(v_)-min(v_)):max(v_);
        D_=drag(rho(i),v_space,S,C_Do,K,W);
        P_=D_.*v_space;
        P_ex=P(i)-P_;
        plot(v_space,P_ex);
        
    end
end

hold off;
xlim([0 80]);
ylim([0 60000]);
xlabel("v (m/s)");
ylabel("P (W)");
legend({'h=0 m','h=3000 m','h=6000 m','h=9000 m','h=12000 m'},'Location','northeast','NumColumns',2);
title('Excess Power at different altitudes');

Ec_max=((3*C_Do/K)^0.75)/(4*C_Do);
P_rmin=(Ec_max^-1)*(2*W^3./(rho.*S)).^0.5;
figure;
grid on;
hold on;
plot(P_rmin,h);
plot(P,h);
plot(P-P_rmin,h);
hold off;
xlim([0 80000]);
ylim([0 20000]);

xlabel("P (W)");
ylabel("h (m)");
legend({'P_{r,min}','P_a','P_{ex,max}'},'Location','northeast');
title('Power at different altitudes');

%equivalent airspeed
v_P_eq=v_P.*(sig).^0.5;

%ceiling calculations
sig_ceil=((Ec_max*P_sl)^-1)*(2*W^3/(1.225*S))^0.5;
h_ceil=siginv(sig_ceil);
rho_ceil=1.225*sig_ceil;
P_ceil=P_sl*(sig_ceil^0.5);
v_ceil=(2*W/(S*rho_ceil*(3*C_Do/K)^0.5)).^0.5;
v_ceil_eq=v_ceil*(sig_ceil)^0.5;
v_ceil_2=vpasolve(0.25*(rho_ceil^2)*S*C_Do*v^4-0.5*rho_ceil*P_ceil*v+(K*W^2)/S == 0,v);

%plotting envelope
v_stall_plot=v_stall;
h1=h;
v_P(i:end,:)=[];
h1(i:end)=[];
v_stall_plot(i:end)=[];
h2_plot=cat(1,h1,h_ceil,flip(h1));
v_P_plot=cat(1,v_P(:,1),v_ceil,flip(v_P(:,2)));

figure;
plot(v_P_plot,h2_plot);
hold on;
plot(v_stall_plot,h1);
plot([0 v_ceil],max(h1)*ones(2,1),'b--');
plot(v_ceil*ones(2,1),[0 h_ceil],'m--');
hold off;
grid on;
xlim([0 100]);
ylim([0 17500]);
xlabel("v (m/s)");
ylabel("h (m)");
legend({'h-v envelope','V_{stall}','h_{ceil}','v_{ceil}'},'Location','northwest','NumColumns',2);
title('Envelope for propeller engine');

%plotting envelope for equivalent airspeed
v_stall_eq_plot=v_stall_eq;
h1=h;
v_P_eq(i:end,:)=[];
h1(i:end)=[];
v_stall_eq_plot(i:end)=[];
h1_plot=cat(1,h1,h_ceil,flip(h1));
v_eq_plot=cat(1,v_P_eq(:,1),v_ceil_eq,flip(v_P_eq(:,2)));

figure;
plot(v_eq_plot,h1_plot);
hold on;
plot(v_stall_eq_plot,h1);
plot([0 v_ceil_eq],h_ceil*ones(2,1),'b--');
plot(v_ceil_eq*ones(2,1),[0 h_ceil],'m--');
hold off;
grid on;
xlim([0 100]);
ylim([0 17500]); 
xlabel("v (m/s)");
ylabel("h (m)");
legend({'h-v envelope','V_{stall}','h_{ceil}','v_{ceil}'},'Location','northwest','NumColumns',2);
title('Envelope for propeller engine (equivalent airspeed)');

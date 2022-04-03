% % jet engine

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

%Thrust model
T_sl=1140;
T=T_sl.*(sig.^(1/3));

% condition for real roots
check=(T/W).^2 >= 4*K*C_Do ;

% vector to store velocity
v_T=zeros(length(h),2);

figure;
hold on;
grid on;
%solving the biquadratic equation and plotting thrust plots
syms v;
for i=1:1:length(h)
    v_=vpasolve(0.25*(rho(i)^2)*S*C_Do*v^4-0.5*rho(i)*T(i)*v^2+(K*W^2)/S == 0,v,[0 Inf]);
    if isempty(v_)
        break;
    end
    v_T(i,:)=v_;
    
    if mod(h(i),3000)==0
        v_space=min(v_):0.025*(max(v_)-min(v_)):max(v_);
        D_=drag(rho(i),v_space,S,C_Do,K,W);
        T_ex=T(i)-D_;
        plot(v_space,T_ex);
%        plot(v_space,D_); %toggle comment for second graph
        
    end
end

% % toggle comment for second graph
% T_rmin=2*W*(C_Do*K)^0.5;
% plot([0 100],T_rmin*ones(2,1),'b--'); 
hold off;
xlim([0 100]);
ylim([0 650]);
% ylim([500 1350]);
xlabel("v (m/s)");
ylabel("T (N)");
legend({'h=0 m','h=3000 m','h=6000 m','h=9000 m','h=12000 m','h=15000 m'},'Location','northeast','NumColumns',2);
title('Excess thrust at different altitudes');
% % toggle comment to produce second graph
% legend({'h=0 m','h=3000 m','h=6000 m','h=9000 m','h=12000 m','h=15000 m','T_{r,min}'},'Location','northeast','NumColumns',3);
% title("Thrust required at different altitudes");        

%equivalent airspeed        
v_T_eq=v_T.*(sig).^0.5;

%ceiling calculations
T_ceil=W*(4*C_Do*K)^0.5;
sig_ceil=(T_ceil/T_sl)^3;
h_ceil=siginv(sig_ceil);
rho_ceil=1.225*sig_ceil;
v_ceil=(T_ceil/(rho_ceil*C_Do*S))^0.5;
v_ceil_eq=v_ceil*(sig_ceil)^0.5;

%plotting envelope
v_stall_plot=v_stall;
h1=h;
v_T(i:end,:)=[];
h1(i:end)=[];
v_stall_plot(i:end)=[];
h1_plot=cat(1,h1,h_ceil,flip(h1));
v_plot=cat(1,v_T(:,1),v_ceil,flip(v_T(:,2)));

figure;
plot(v_plot,h1_plot);
hold on;
plot(v_stall_plot,h1);
plot([0 v_ceil],h_ceil*ones(2,1),'b--');
plot(v_ceil*ones(2,1),[0 h_ceil],'m--');
hold off;
grid on;
xlim([0 100]);
ylim([0 19500]);
xlabel("v (m/s)");
ylabel("h (m)");
legend({'h-v envelope','V_{stall}','h_{ceil}','v_{ceil}'},'Location','northwest','NumColumns',2);
title('Envelope for jet engine');

%plotting envelope for equivalent airspeed
v_stall_eq_plot=v_stall_eq;
h1=h;
v_T_eq(i:end,:)=[];
h1(i:end)=[];
v_stall_eq_plot(i:end)=[];
h1_plot=cat(1,h1,h_ceil,flip(h1));
v_eq_plot=cat(1,v_T_eq(:,1),v_ceil_eq,flip(v_T_eq(:,2)));

figure;
plot(v_eq_plot,h1_plot);
hold on;
plot(v_stall_eq_plot,h1);
plot([0 v_ceil_eq],h_ceil*ones(2,1),'b--');
plot(v_ceil_eq*ones(2,1),[0 h_ceil],'m--');
hold off;
grid on;
xlim([0 100]);
ylim([0 19500]);
xlabel("v (m/s)");
ylabel("h (m)");
legend({'h-v envelope','V_{stall}','h_{ceil}','v_{ceil}'},'Location','northwest','NumColumns',2);
title('Envelope for jet engine (equivalent airspeed)');



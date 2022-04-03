% % propeller engine

%given constant parameters
m=750;
S=12;
b=10;
C_Do=0.036;
C_lmax=2.7;
e=0.87;

% altitude and air density
h=(0:25:20000)';
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
eta=1;

% climb  calculations
v_roc_max=(2*W./(rho.*S*(3*C_Do/K)^0.5)).^0.5;
Y_roc_max=asind((eta.*P./(v_roc_max.*W))-(0.5*rho.*v_roc_max.^2*S*C_Do/W)-((2*K*W)./(rho.*S.*v_roc_max.^2)));
roc_max=v_roc_max.*sind(Y_roc_max);
v_Y_max=(4*W^2*K)./(P.*rho.*eta*S);
Y_max=asind((eta.*P./(v_Y_max.*W))-(0.5*rho.*v_Y_max.^2*S*C_Do/W)-((2*K*W)./(rho.*S.*v_Y_max.^2)));

% ceiling calculations
roc_max_f=@(hq) interp1(h,roc_max,hq);
roc_abs_ceil=0;
roc_ser_ceil=100*0.00508;

abs_ceil=fzero(@(hq) roc_max_f(hq)-roc_abs_ceil,[10000 16000]);
ser_ceil=fzero(@(hq) roc_max_f(hq)-roc_ser_ceil,[8000 16000]);

% plotting
figure;
plot(roc_max,h);
hold on;
plot([0 8],abs_ceil*ones(2,1),'m--');
plot([0 roc_ser_ceil],ser_ceil*ones(2,1),'b--');
plot(roc_ser_ceil*ones(2,1),[0 ser_ceil],'g--');
hold off;
grid on;
xlim([0 8]);
ylim([0 19500]);
xlabel("Maximum rate of climb (m/s)");
ylabel("h (m)");
legend({'ROC_{max} vs h','absolute ceiling','service ceiling','ROC_{max} for service ceiling'},'Location','northeast','NumColumns',2);
title('maximum rate of climb vs altitude for propeller engine');

figure;
plot(roc_max,v_roc_max);
grid on;
xlim([0 8]);
xlabel("Maximum rate of climb (m/s)");
ylabel("velocity at maximum rate of climb (m/s)");
legend('ROC_{max} vs V_{ROC_{max}}','Location','northeast','NumColumns',2);
title('maximum rate of climb vs velocity for propeller engine');

figure;
plot(Y_max,v_Y_max);
grid on;
 xlim([0 30]);
xlabel("Maximum climb angle (deg)");
ylabel("velocity at maximum climb angle (m/s)");
legend('\gamma_{max} vs V_{\gamma_{max}}','Location','northeast','NumColumns',2);
title('maximum climb angle vs velocity for propeller engine');

v_P=zeros(length(h),2);

% solving the biquadratic equation and plotting climb schedule
figure;
hold on;
grid on;
syms v;
Ps=0;
for j=0:1:5
    Ps=j;
    for i=1:1:length(h)
        v_=vpasolve(P(i)*eta-Ps*W-0.5*C_Do*rho(i)*v^3*S-(2*K*W^2)/(rho(i)*v*S) == 0,v,[0 Inf]);
        if isempty(v_)
            break;
        end
        v_P(i,:)=v_;
    end

%plotting envelope

    h1=h;
    v_P1=v_P;
    v_P1(i:end,:)=[];
    h1(i:end)=[];

    h1_plot=cat(1,h1,flip(h1));
    v_plot=cat(1,v_P1(:,1),flip(v_P1(:,2)));

    if j==0
        v_stall_plot=v_stall;
        v_stall_plot(i:end)=[];
        plot(v_stall_plot,h1);
    end
    plot(v_plot,h1_plot);
end

% plotting climb schedule
hold off;
xlim([0 100]);
ylim([0 19500]);
xlabel("v (m/s)");
ylabel("h (m)");
legend({'v_{stall}','P_s='+string(0)+'m/s','P_s='+string(1)+'m/s','P_s='+string(2)+'m/s','P_s='+string(3)+'m/s','P_s='+string(4)+'m/s','P_s='+string(5)+'m/s'},...
'Location','northwest','NumColumns',2);
title('Climb schedule for propeller engine');
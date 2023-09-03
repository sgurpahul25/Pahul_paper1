clear
close all

n=10;

sigma1 = [0 1;1 0];   
sigma2 = [0 complex(0,-1);complex(0,1) 0];
sigma3 = [1 0;0 -1];

coeff = rand(1,n);
coeff = coeff/sum(coeff);
coeff = sqrt(coeff);

vecs = eye(n);

psi = zeros(n,1);
for i = 1:n
    psi = psi + coeff(i)*vecs(:,i);
end

%psi = 1/sqrt(2)*([1;0] + [0;1]);
rho = psi*psi';
rho_0 = rho;

N=1000;
ts = linspace(0,10,N);
h = 0.05;
ti = 6;
tf = 8;
initial_state_prob = zeros([1 N]);
orthogonal_state_prob = zeros([1 N]);
initial_state_prob(1,1) = trace(rho*rho); 
orthogonal_state_prob(1,1) = trace(rho*(eye(n)-rho)); 

prob_eigenvals = zeros(n,N);
if n==2
    coord =zeros([N 3]);
    coord(1,1) = trace(rho*sigma1);
    coord(1,2) = trace(rho*sigma2);
    coord(1,3) = trace(rho*sigma3);
end
prob_first_eigenval(1,1) = trace(vecs(:,1)'*rho*vecs(:,1));
prob_second_eigenval(1,1) = trace(vecs(:,2)'*rho*vecs(:,2));

eigenvals = 10*rand(1,n);
% epsilon = eigenvals(1);
H = zeros(n,n);
%m = prob_func(coeff);
for k=1:n
     H = H +eigenvals(k)*(vecs(:,k)*vecs(:,k)');
end
%H = diag([1 2 3 4]);
gamma = 3;
P = vecs(:,1)*vecs(:,1)';
%P = diag([4 3 2 1]);
for i = 1:N
%     if i < ti || i>tf
%         M = expm(-complex(0,1)*H*h)*(rho*expm(complex(0,1)*H'*h));
%     else 
%         H(m,m) = complex(0,1)*epsilon;
%         %H(2,2) = -complex(0,1)*0.1*epsilon;
%         M = expm(-complex(0,1)*H*h)*rho*expm(complex(0,1)*H'*h);
%         %M = (cos(p2*t(i)))^2*rho + complex(0,1)*cos(p2*t(i))*sin(p2*t(i))*(rho*H'-H*rho)/p2 + (sin(p2*t(i)))^2*H*rho*H'/p2^2;
%     end
    Ht = H + complex(0,1)*(gamma/2)*(tanh(gamma*(ts(i)-ti)) - tanh(gamma*(ts(i)-tf)))*P;
    M = expm(-complex(0,1)*Ht*h)*(rho*expm(complex(0,1)*Ht'*h));
    rho_t = M/trace(M);
    disp(Ht)
    if n==2
        coord(i,1) = trace(rho_t*sigma1);
        coord(i,2) = trace(rho_t*sigma2);
        coord(i,3) = trace(rho_t*sigma3);
    end
    initial_state_prob(1,i) = trace(rho_t*rho_0);
    orthogonal_state_prob(1,i) = trace(rho_t*(eye(n)-rho_0));
    for k=1:n
        prob_eigenvals(k,i) = trace(vecs(:,k)'*rho_t*vecs(:,k));
    end
    rho = rho_t;

end

if n==2
    figure;
    hold on
    plot3(coord(:,1), coord(:,2), coord(:,3), LineWidth=2.5, color = 'blue')
    plot3(0,0,1, '.', color = 'red', MarkerSize= 15)
    plot3(0,0,-1, '.', color = 'green', MarkerSize= 15)
    [X,Y,Z] = sphere;
    surf(X, Y,Z, FaceColor="none", EdgeLighting="flat")
    xlabel("x","FontSize",12)
    ylabel("y","FontSize",12)
    zlabel("z","FontSize",12)
    %title("\gamma_m = "+gamma)
    set(gca, 'CameraPosition', [16.240247520153392 -4.499944292577442 4.000607684821305]);
    set(gcf, 'Position', [100,100,250,250])
    hold off
end
figure;
hold on
plot(ts, initial_state_prob, LineWidth=2)
plot(ts, orthogonal_state_prob, LineWidth=2)
plot([ti ti], [0 1], Color='black', LineWidth=1)
% plot([tf/100 tf/100],[0 1], Color='black')
xlabel("t", FontSize=12)
ylabel("Transition Probability", FontSize=12)
lgd1 = legend("Initial state", "Orthogonal state");
lgd1.FontSize = 11.5;
set(gcf, 'Position', [400,400,470.6,312.6])
hold off
figure; 
hold on 
for i=1:n
    plot(ts, prob_eigenvals(i,:), LineWidth=2)
end
plot([ti ti],[-0.03 1.03], Color ='black', LineWidth=1)
%plot([tf tf],[-0.03 1.03], Color ='black', LineWidth=0.7)
lgd = legend("1", "2", "3");
lgd.FontSize = 11.5;
ylabel("Eigenstate Probability", FontSize=12)
xlabel("t", FontSize = 12)
%t = text(5.6,0.07,'t_i', 'FontSize',12);
%t2 = text(8.1,0.07,'t_f', 'FontSize',12); 
xlim([0.04 10])
set(gcf, 'Position', [400,400,470.6,312.6])
ylim([-0.03 1.03])
%box on
%title("\gamma = "+gamma)
hold off


% figure;
% Y = fft(initial_state_prob(1:ti));
% P2 = abs(Y/ti);
% P1 = P2(1:ti/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = 100*(0:(ti/2))/ti;
% plot(f(1:200),P1(1:200)) 


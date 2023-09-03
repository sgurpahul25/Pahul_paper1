sigma1 = [0 1;1 0];   
sigma2 = [0 complex(0,-1);complex(0,1) 0];
sigma3 = [1 0;0 -1];

% gamma = 0;
% K = sigma1 -complex(0,1)*gamma*sigma3;
% [vec, val] = eig(K);
% vecr1 = vec(:,1);
% vecr2 = vec(:,2);
vecr1 = [1 0]';
vecr2 = [0 1]';

N=1000;
h = 0.01;
gamma1 = 0.5; 

ts = zeros(N,1);
ti = 7;
tf = 8;
dist = zeros(N-1,1);
rx = 0;
rz = 0; 
ry = sqrt(1-rx^2-rz^2);
rho = 0.5*(eye(2) + rx*sigma1 + ry*sigma2 + rz*sigma3);
rho_0 = rho;
coord =zeros([N 3]);
coord(1,1) = trace(rho*sigma1);
coord(1,2) = trace(rho*sigma2);
coord(1,3) = trace(rho*sigma3);
dist_from_start = zeros(N,1);
prob_first_eigenval(1,1) = trace(vecr1'*rho*vecr1);
%prob_first_eigenval(1,1) = vecr1'*rho*vecr1 + (vecr1'*rho*vecr2)*(vecr1'*vecr2) + (vecr2'*vecr1)*(vecr1'*rho*vecr1)*(vecr2'*vecr1) + (vecr2'*vecr1)*(vecr1'*rho*vecr2);
prob_second_eigenval(1,1) = trace(vecr2'*rho*vecr2);
initial_state_prob(1,1) = trace(rho*rho); 
final_state_prob(1,1) = trace(rho*(eye(2)-rho)); 

for i = 2:N
    ts(i,1) = h*(i-1);
    p1 = sqrt(1-(gamma1^2/2)*(tanh(gamma1*(ts(i)-ti))-tanh(gamma1*(ts(i)-tf))));
    H = p1*[1 0;0 -1];
    disp(p1)

%         p2 = sqrt(gamma^2 - 1);
%         beta = 2*(p2*sqrt(2*gamma*(gamma-p2))+gamma*(gamma+p2));
%         p2 = complex(0,1)*p2;
%         H = [p2 beta;0 -p2];
        %H = [-p2 0;-beta p2];
    M = expm(-complex(0,1)*H*h)*rho*expm(complex(0,1)*H'*h);
    rho_t = M/trace(M);
    coord(i,1) = trace(rho_t*sigma1);
    coord(i,2) = trace(rho_t*sigma2);
    coord(i,3) = trace(rho_t*sigma3);
    dist(i-1,1) = acos(coord(i-1,1)*coord(i,1) + coord(i-1,2)*coord(i,2) + coord(i-1,3)*coord(i,3));
    dist_from_start(i,1) = acos(rx*coord(i,1) + ry*coord(i,2) + rz*coord(i,3));
    initial_state_prob(1,i) = trace(rho_t*rho_0);
    final_state_prob(1,i) = trace(rho_t*(eye(2)-rho_0));
        %prob_first_eigenval(1,i) = vecr1'*rho*vecr1 + (vecr1'*rho*vecr2)*(vecr1'*vecr2) + (vecr2'*vecr1)*(vecr1'*rho*vecr1)*(vecr2'*vecr1) + (vecr2'*vecr1)*(vecr1'*rho*vecr2);
    prob_first_eigenval(1,i) = trace(vecr1'*rho_t*vecr1);
    prob_second_eigenval(1,i) = trace(vecr2'*rho_t*vecr2);
    rho = rho_t;

%     if j>2 && dist(i-1,1)/h <10^(-2)
%       break 
%     end
end

% 
figure;
hold on
plot3(coord(:,1), coord(:,2), coord(:,3),linewidth =3 , color='blue')
% plot3(fixed_points(1,1),fixed_points(2,1), fixed_points(3,1), '.', MarkerSize=22, Color='red')
% plot3(fixed_points(4,1),fixed_points(5,1), fixed_points(6,1), '.', MarkerSize=22, Color='red')
[X,Y,Z] = sphere;
surf(X, Y,Z, FaceColor="none", EdgeLighting="flat")
hold off

figure;
hold on
plot(ts(2:end,1), dist(:,1)/h, LineWidth=2)
legend("\gamma = "+gamma1)
xlabel("t")
ylabel("speed")

figure;
hold on
plot(ts,initial_state_prob, LineWidth=1.5)
plot(ts, final_state_prob, LineWidth=1.5)
plot([ti ti],[0 1],Color='black')
plot([tf tf],[0 1], Color='black')
xlabel("t")
ylabel("Probability")
legend("Probability of being in the initial state", "Probability of being in the orthogonal state")
title("r_x ="+ rx + " r_y = "+ry)
hold off
% hold off
% figure;
% hold on
% plot(ts(:,1), dist_from_start(:,1), LineWidth=1.5)
% legend("r_x = "+rx)
% xlabel("t")
% ylabel("distance")
% title("\gamma = "+gamma1)
% hold off
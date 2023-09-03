%Pauli matrices
sigma1 = [0 1;1 0];   
sigma2 = [0 complex(0,-1);complex(0,1) 0];
sigma3 = [1 0;0 -1];

sigma = [sigma1 sigma2 sigma3];
gamma = 1.5;
K = sigma3 + complex(0,1)*gamma*(eye(2)+sigma3)/2;
[vec, val] = eig(K);
vecr1 = vec(:,1);
vecr2 = vec(:,2);

rx = -0.6;
ry = 0.5; 
rz = sqrt(1-rx^2-ry^2);
rho = 0.5*(eye(2) + rx*sigma1 + ry*sigma2 + rz*sigma3);
p = sqrt(1-gamma^2);
rho_0 = rho;
R = [1 0 -complex(0,1)*gamma];
r = [rx ry rz];

N=1000;
h = 0.01;
t = h*ones(N-1);
ts = zeros(N);
ti = 1000;
tf = 750;
initial_state_prob = zeros([1 N]);
final_state_prob = zeros([1 N]);
prob_first_eigenval = zeros([1 N]);
prob_second_eigenval = zeros([1 N]);
coord =zeros([N 3]);
coord(1,1) = trace(rho*sigma1);
coord(1,2) = trace(rho*sigma2);
coord(1,3) = trace(rho*sigma3);
prob_first_eigenval(1,1) = trace(vecr1'*rho*vecr1);
%prob_first_eigenval(1,1) = vecr1'*rho*vecr1 + (vecr1'*rho*vecr2)*(vecr1'*vecr2) + (vecr2'*vecr1)*(vecr1'*rho*vecr1)*(vecr2'*vecr1) + (vecr2'*vecr1)*(vecr1'*rho*vecr2);
prob_second_eigenval(1,1) = trace(vecr2'*rho*vecr2);
initial_state_prob(1,1) = trace(rho*rho); 
final_state_prob(1,1) = trace(rho*(eye(2)-rho)); 
prob_trans = zeros([1 N]);
prob_trans(1,1) = trace(rho*rho); 
%disp(coord(1,1:3))
for i = 2:N
    %rho_t = (exp(-complex(0,1)*K*t(i))*rho*exp(complex(0,1)*K'*t(i)))/trace(exp(-complex(0,1)*K*t(i))*rho*exp(complex(0,1)*K'*t(i)));
    ts(i) = h*(i-1);
    if i < ti || i>tf
        M = expm(-complex(0,1)*K*h)*rho*expm(complex(0,1)*K'*h);
        N = expm(-complex(0,1)*K*ts(i))*rho_0*expm(complex(0,1)*K'*ts(i));
    else 
        gamma = 3;
        p2 = complex(0,1)*sqrt(gamma^2 - 1);
        H = p2*[0 p1;p1 0];
        %M = (cos(p2*t(i)))^2*rho + complex(0,1)*cos(p2*t(i))*sin(p2*t(i))*(rho*H'-H*rho)/p2 + (sin(p2*t(i)))^2*H*rho*H'/p2^2;
    end
    rho_t = M/trace(M);
    coord(i,1) = trace(rho_t*sigma1);
    coord(i,2) = trace(rho_t*sigma2);
    coord(i,3) = trace(rho_t*sigma3);
    initial_state_prob(1,i) = trace(rho_t*rho_0);
    final_state_prob(1,i) = trace(rho_t*(eye(2)-rho_0));
    %prob_first_eigenval(1,i) = vecr1'*rho*vecr1 + (vecr1'*rho*vecr2)*(vecr1'*vecr2) + (vecr2'*vecr1)*(vecr1'*rho*vecr1)*(vecr2'*vecr1) + (vecr2'*vecr1)*(vecr1'*rho*vecr2);
    prob_first_eigenval(1,i) = trace(vecr1'*rho_t*vecr1);
    prob_second_eigenval(1,i) = trace(vecr2'*rho_t*vecr2);
    rho = rho_t;
    %prob_trans(1,i) = 0.5+(0.5/trace(N))*((complex(0,1)*cos(p*ts(i))*sin(p*ts(i))/p)*(r*R'-R*r')+ cos(p*ts(i))*cos(p*ts(i)) + (sin(p*ts(i))*sin(p*ts(i))/p^2)*(complex(0,1)*cross(R,conj(R))*r'+(R*r')*(r*R')-dot(cross(R,r),cross(conj(R),r))));
    prob_trans(1,i) = rx/trace(N);
end

%figure;
hold on
plot3(coord(:,1), coord(:,2), coord(:,3), LineWidth=2.5, Color="blue")
[X,Y,Z] = sphere;
surf(X, Y,Z, FaceColor="none", EdgeLighting="flat")
title("\gamma = "+ gamma)
xlabel("x", FontSize=13)
ylabel("y", FontSize=13)
zlabel("z", FontSize=13)
hold off
figure;
hold on
plot(ts,initial_state_prob, LineWidth=1.5)
plot(ts, final_state_prob, LineWidth=1.5)
plot([ti/100 ti/100],[0 1],Color='black')
plot([tf/100 tf/100],[0 1], Color='black')
xlabel("t")
ylabel("Probability")
legend("Probability of being in the initial state", "Probability of being in the orthogonal state")
title("r_x ="+ rx + " r_y = "+ry)
hold off
figure; 
hold on 
plot(ts, prob_first_eigenval, LineWidth=1.5)
plot(ts, prob_second_eigenval, LineWidth=1.5)
plot([ti/100 ti/100],[0 1], Color='black')
plot([tf/100 tf/100],[0 1], Color ='black' )
xlabel("t")
ylabel("Probability")
legend("Probability of being in 1st eigenstate", "Probability of being in second eigen state")
title("r_x ="+ rx + " r_y = "+ry + "\gamma = "+ gamma)
hold off
% % 
% figure;
% plot(ts, coord(:,3));
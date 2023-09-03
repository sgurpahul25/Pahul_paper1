n=8;

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

rho = psi*psi';
rho_0 = rho;

N=3000;
h = 0.01;
ts = zeros(N);
ti = 2000;
tf = 750;
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

%eigenvals = 1:10;
eigenvals = 10*rand(1,n);
epsilon = eigenvals(1);
H = zeros(n,n);
m = prob_func(coeff);
for k=1:n
    H = H +eigenvals(k)*(vecs(:,k)*vecs(:,k)');
end

for i = 2:N
    if i < ti 
        M = expm(-complex(0,1)*H*h)*(rho*expm(complex(0,1)*H'*h));
    else 
        H(m,m) = complex(0,1)*epsilon;
        %H(2,2) = -complex(0,1)*0.1*epsilon;
        M = expm(-complex(0,1)*H*h)*rho*expm(complex(0,1)*H'*h);
        %M = (cos(p2*t(i)))^2*rho + complex(0,1)*cos(p2*t(i))*sin(p2*t(i))*(rho*H'-H*rho)/p2 + (sin(p2*t(i)))^2*H*rho*H'/p2^2;
    end
    rho_t = M/trace(M);
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
    ts(i) = h*(i-1);

end

if n==2
    figure;
    hold on
    plot3(coord(:,1), coord(:,2), coord(:,3), LineWidth=3)
    [X,Y,Z] = sphere;
    surf(X, Y,Z, FaceColor="none", EdgeLighting="flat")
    hold off
end
% figure;
% hold on
% plot(ts,initial_state_prob, LineWidth=1.5)
% plot(ts, orthogonal_state_prob, LineWidth=0.1)
% plot([ti/100 ti/100],[0 1],Color='black')
% % plot([tf/100 tf/100],[0 1], Color='black')
% xlabel("t")
% ylabel("Transition Probability")
% legend("Probability of being in the initial state", "Probability of being in the orthogonal state")
% hold off
figure; 
hold on 
for i=1:n
    plot(ts, prob_eigenvals(i,:), LineWidth=1.5)
end
plot([ti/100 ti/100],[0 1], Color='black')
% plot([tf/100 tf/100],[0 1], Color ='black' )
ylabel("Eigenstate Probability")
xlabel("t")
xlim([0.04 10])
set(gcf, 'Position', [400,400,417.6000000000001,268.8])
hold off



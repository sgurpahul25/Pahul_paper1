function m = prob_func(coeff)
    n = length(coeff);
    r = rand(1,1);
    disp(r)
    cumu_coeff= zeros([1 n+1]);
    for i=2:n+1
        cumu_coeff(1,i) = cumu_coeff(1,i-1)^2+coeff(i-1)^2;
        if cumu_coeff(1,i-1)<r && r<cumu_coeff(i)
            break
        end
    
    end
    m = i-1;
end
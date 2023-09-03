function T = tanh_t(t,g,ti,tf)
    T = tanh(g*(t - ti)) - tanh(g*(t-tf)); 
end
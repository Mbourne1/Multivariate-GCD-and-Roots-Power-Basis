function [f,g,u,v,m,n,n_t,m_t] = Examples(ex_num)

switch ex_num
    
    case 1
        
        f = [-25  30  5  -5  -6  0  0  1  0  0];
        g = [10  -12  5  2  -6  0  0  1  0  0];
        v = [2 0 1];
        u = [-5 0 1];

    case 2
        f = [12 -18 10 6 -15 2 0 5 -3 0 0 0 1 0 0];
        g = [-105 90 -56 8 48 -7 0 8 6 0 0 0 1 0 0];
        v = [35 5 7 0 1 0];
        u = [-4 2 -2 0 1 0];
        
    case 3
        f = [24 -36 -14 12 21 2 0 -7 -3 0 0 0 1 0 0];
        g = [-10 15 2 -5 -3 0 0 1 0 0];
        v = [-5 0 1];
        u = [12 0 -7 0 0 1];

    case 4
        f = [-12 22 18 -12 -33 -6 2 18 11 0 0 -3 -6 0 0 0 0 1 0 0 0];
        g = [-30 55 36 -30 -66 -6 5 36 11 0 0 -6 -6 0 0 0 0 1 0 0 0];
        u = [-2 0 1];
        v = [-5 0 1];
        

        
end

rt_f = roots([1 1 -(2*length(f))]);
rt_f = rt_f(rt_f>=0);
m = rt_f-1;

rt_g = roots([1 1 -(2*length(g))]);
rt_g = rt_g(rt_g>=0);
n = rt_g-1;

rt_v = roots([1 1 -(2*length(v))]);
rt_v = rt_v(rt_v>=0);
n_t = rt_v-1;

rt_u = roots([1 1 -(2*length(u))]);
rt_u = rt_u(rt_u>=0);
m_t = rt_u-1;
end


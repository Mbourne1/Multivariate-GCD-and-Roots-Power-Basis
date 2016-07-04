function [] = Untitled()
syms x y

f0 = (x+1) * (x+2)^2 * (x+3)^3 * (y+1) * (y+2)^2 * (y+3)^3;
f1 = (x+2) * (x+3)^2 * (y+2) * (y+3)^2;
f2 = (x+3) * (y+3);

m0 = double(feval(symengine, 'degree', f0));
m1 = double(feval(symengine, 'degree', f1));
m2 = double(feval(symengine, 'degree', f2));

f0 = double(rot90(coeffs(f0,[x,y],'All'),2));
f1 = double(rot90(coeffs(f1,[x,y],'All'),2));
f2 = double(rot90(coeffs(f2,[x,y],'All'),2));

arr_f{1} = f0;
arr_f{2} = f1;
arr_f{3} = f2;

vDeg_fxy = [m0; m1; m2];

hi = Deconvolve_Bivariate_Batch(arr_f,vDeg_fxy);

end
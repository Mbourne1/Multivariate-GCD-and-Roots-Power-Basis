
% Build Subresultant of noisy unprocessed polynomials
S_Unproc = BuildSylvesterMatrix(fxy,gxy,0,0,1,1,1);
vSingularValues_S_Unproc = svd(S_Unproc);
vSingularValues_S_Unproc  = normalise(vSingularValues_S_Unproc);

% Plot singular values of fxy_working and fxy_output
figure('name','Singular Values of S with and without SNTLN')
hold on

switch LOW_RANK_APPROXIMATION_METHOD
    case {'Standard SNTLN', 'Standard STLN'}
        vSingularValues_S_LowRank = svd(S_LowRankApprox);
        vSingularValues_S_LowRank = normalise(vSingularValues_S_LowRank);
        plot(log10(vSingularValues_S_LowRank),'-*','DisplayName','S(f(x,y),g(x,y)) SNTLN')
    case {'None'}
    otherwise
        error('err')
        
end

switch BOOL_PREPROC
    case 'y'
        vSingularValues_S_Preproc = svd(S_Preproc);
        vSingularValues_S_Preproc = normalise(vSingularValues_S_Preproc);
        plot(log10(vSingularValues_S_Preproc),'-o','DisplayName','S(f(\omega_{1},\omega_{2}),g(\omega_{1},\omega_{2})) Preprocessed')
    case 'n'
    otherwise
        error('err')
end

plot(log10(vSingularValues_S_Unproc),'-*','DisplayName','S(f(x,y),g(x,y)) Noisy')
legend(gca,'show');
ylabel('log_{10} \sigma_{i} / \sigma_{1}')
xlabel('i')


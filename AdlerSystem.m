function dPhi = AdlerSystem(Phi, ParamVec)

D = [pi/2 - ParamVec.Sigma, pi/2 + ParamVec.Sigma];

dPhi = zeros(2, 1);

%Разрывная функция связи
% dPhi(1) = ParamVec.Gamma(1) - sin(Phi(1)) - ParamVec.d * Ffunc(Phi(2), D);
% dPhi(2) = ParamVec.Gamma(2) - sin(Phi(2)) - ParamVec.d * Ffunc(Phi(1), D);

% Непрерывная функция связи
k = -500;
Alpha = D(1);
Delta = 2 .* ParamVec.Sigma;
dPhi(1) = ParamVec.Gamma(1) - sin(Phi(1)) - ParamVec.d .* F2func(Phi(2), k, Delta, Alpha);
dPhi(2) = ParamVec.Gamma(2) - sin(Phi(2)) - ParamVec.d .* F2func(Phi(1), k, Delta, Alpha);

% Sigma = ParamVec.Sigma;
% 
% Center = ParamVec.Gamma(1) - ParamVec.d / 2;
% if Center > 1
%     Center = pi / 2;
% elseif Center < -1
%         Center = 3 * pi / 2;
% else
%     Center = asin(Center);
% end
% 
% D = [Center - Sigma, Center + Sigma];
% 
% dPhi = zeros(2, 1);
% 
% %Разрывная функция связи
% % dPhi(1) = ParamVec.Gamma(1) - sin(Phi(1)) - ParamVec.d * Ffunc(Phi(2), D);
% % dPhi(2) = ParamVec.Gamma(2) - sin(Phi(2)) - ParamVec.d * Ffunc(Phi(1), D);
% 
% % Непрерывная функция связи
% k = -500;
% Alpha = D(1);
% Delta = 2 .* ParamVec.Sigma;
% dPhi(1) = ParamVec.Gamma(1) - sin(Phi(1)) - ParamVec.d .* F2func(Phi(2), k, Delta, Alpha);
% dPhi(2) = ParamVec.Gamma(2) - sin(Phi(2)) - ParamVec.d .* F2func(Phi(1), k, Delta, Alpha);

end


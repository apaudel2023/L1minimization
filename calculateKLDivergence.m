% To calculate kl divergence
function klDiv = calculateKLDivergence(xtrue, xpred)
% This function calculates the Kullback-Leibler Divergence
[P, xi] = ksdensity(xtrue);
[Q, ~] = ksdensity(xpred, xi);

% P and Q are two discrete probability distributions
% Ensure P and Q are normalized (sum to 1)
P = P / sum(P);
Q = Q / sum(Q);

% Ensure no zero elements; KL is undefined for P(i) = 0 or Q(i) = 0
epsilon = 1e-10;
P(P == 0) = epsilon;
Q(Q == 0) = epsilon;

% Calculate KL Divergence
klDiv = sum(P .* log(P ./ Q));
end
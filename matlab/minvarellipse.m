% Minimum variance ellipsoid estimator for the mean and covariance
% Jussi Tohka 1st March 2002
% Input : data is the n x p matrix containing n samples for each p variates. 
% Algorithm (approximative)is described in page 259 onwards in
% P.J. Rousseeuw and A.M. Leroy: Robust Regression and Outlier Detection
% John Wiley & Sons 1987

function [mu,sigma] = minvarellipse(data);

    sz = size(data);
    if sz(1) < sz(2)
      data = data';
      sz = size(data);
    end

    n = sz(1);
    p = sz(2);
    
  
    % select number of trials if not given.

    trials = min(n*2,10000);
    score = zeros(trials,1);
    scale = zeros(trials,1);
    samples = zeros(trials,p + 1);

    for i = 1:trials
      % draw a random subsample
      % samples(i,:) = unidrnd(n,1,p + 1);
      samples(i,:) = drawsample(n,p + 1);
      % calculate the mean and covariance of the sample
      m = mean(data(samples(i,:),:));
      c = cov(data(samples(i,:),:));
      if rcond(c) < 3*eps
        score(i) = Inf;
      else
        % deflate or inflate to contain half of the data points
        cinv = inv(c);
        standata = data - repmat(m,n,1);
        scale(i) = median(sum(standata*cinv.*standata,2));
        % calculate the score i.e volume of the ellipsoid
        score(i) = sqrt(det(scale(i)*c));
      end
    end

    [minimum,index] = min(score);
    mu = mean(data(samples(index,:),:));
    sigma = scale(index)*cov(data(samples(index,:),:))*((1/chi2inv(0.5,p)));
    

    function s = drawsample(datalen,sz)
        s = zeros(1,sz);
        s(1) = unidrnd(datalen);
        for i = 2:sz
          s(i) = unidrnd(datalen - i);
          s(i) = s(i) + sum((s(i) >= s(1:i -1)));
        end
          

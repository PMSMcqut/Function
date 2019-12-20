function [X_rms]  = RMS_cal(X)
le = length(X);
X_sum = 0;
for i = 1:le
    X_sum = X_sum + X(i)^2;
end
X_rms = (X_sum/le)^0.5;
    
  
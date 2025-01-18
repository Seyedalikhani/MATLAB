function y = truncated_sinc(t, fsamp, M)
    y = sinc(t * fsamp); % sinc(x) = sin(pi*x)/(pi*x) in MATLAB
    y(abs(t) > M / fsamp) = 0; % Truncate outside [-M*T_s, M*T_s]
end
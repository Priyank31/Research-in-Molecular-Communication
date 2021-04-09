function ans = Q_fun(lambda, n)
    u1 = factorial(n);
    u3 = lambda.^n;
    u2 = exp(-lambda)*u3;
    ans = u2/u1;
end
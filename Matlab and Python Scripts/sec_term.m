function ans1 = sec_term(Co1, avg1)
    t3 = Co1/avg1;
    t2 = log(1+t3); 
    t1 = Co1/t2;  % tao
    
    ans1 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=0:1:t1
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans1 = ans1 + (d1-d2)/factorial(jj);
    end
end
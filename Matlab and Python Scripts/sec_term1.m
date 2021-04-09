function ans2 = sec_term1(Co1, avg1, final_tao)
    ans2 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=0:1:final_tao
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans2 = ans2 + (d1-d2)/factorial(jj);
    end
end
%Import lukji, forwardrow and backwardrow from linear_systems_direct.m
function[x]=MEG(A,b)
        [n,~]=size(A);
        LU=lukji(A);
        L=tril(LU);
        for i=1:n 
            L(i,i)=1; 
        end
        y=forwardrow(L,b);
        x=backwardrow(triu(LU),y);
end

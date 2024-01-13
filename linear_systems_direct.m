function[x]=forwardrow(L,b)
%Forward substitution for lower triangular matrices
%Input: Matrix L and vector b
%Output: solution x
%passo 1: declaration and initialization of x
n=height(L);
if (istril(L)) 
    x=zeros(n,1);
    %passo 2: algorithm start
    x(1)=b(1)/L(1,1);
    for i=2:n
        x(i)=(b(i)-L(i,1:i-1)*x(1:i-1))/L(i,i);
    end
else 
    error('Matrice non è triangolare inferiore')
end
end

function [x]= backwardrow(U,b)
%Input: Matrix U upper triangular, vector b
%Output: solution x with backward substitution
%Passo 1: declaration of x
if (istriu(U))
    n=height(U);
    x=zeros(n,1);
    %Passo 2: algorithm start
    x(n)=b(n)/U(n,n);
    for i=n-1:-1:1
        x(i)=(b(i)-(U(i,i+1:n)*x(i+1:n)))/U(i,i);
    end
else 
    error('Matrice non è triangolare superiore') 
end
end


function [H]=chol2(H)
%function implementing cholesky factorization
if (issymmetric(H)) 
    H(1,1)=sqrt(H(1,1));
    for j=2:height(H)
                if ((H(j,j)<=0)) 
            error ('Error: non positive element on the diagonal')
                end
        for i=1:j-1
            H(i,j)=((H(i,j)-(H(1:i-1,i)'*H(1:i-1,j)))/H(i,i));
        end
        H(j,j)=sqrt(H(j,j)-(H(1:j-1,j)'*H(1:j-1,j)));
    end
else
error('Error: non symmetric matrix')
end
end

function [A] = lukji(A)
%Input: Matrix A
%Output: LU in place factorization
%inferiore stretta di A)
n=height(A);
for k=1:n-1
    if (isequal(A(k,k),0)) 
        error ('Error: at least one element on the diagonal is zero')
    end
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    for j=k+1:n
        for i=k+1:n
            A(i,j)=A(i,j)-A(i,k)*A(k,j);
        end
    end
end
end

function[A]=thomas(A)
%Function implementing thomas algorithm for tridiagonal matrices
    if (isbanded(A,1,1))
        n=height(A);
        for i=1:n-1
            A(i+1,i)=A(i+1,i)/A(i,i);
            A(i+1,i+1)=A(i+1,i+1)-(A(i+1,i)*A(i,i+1));
        end
    else
        error('Error: the matrix is not tridiagonal')
    end
end

function[x]=backwardbid(A,b)
%variant of backwards substitution for bidiagonal matrix
    if(isbanded(A,0,1)) %Controllo sia bidiagonale superiore
        n=height(A);
        x=zeros(n,1); %Prealloco spazio in memoria per velocizzare
        x(n)=b(n)/A(n,n);
        for i=n-1:-1:1
            if (isequal(A(i,i),0))
                error('Elemento diagonale nullo')
            end
            %Applico algoritmo
            x(i)=(b(i)-A(i,i+1)*x(i+1))/A(i,i);
        end
    else 
        error('La matrice non è bidiagonale superiore')
    end
end

function[x]=forwardbid(A,b)
%variant of forwards substitution for bidiagonal matrix
    if (isbanded(A,1,0)) %Controllo che la matrice sia bidiagonale inferiore
        n=height(A);
        x=zeros(n,1); %Prealloco spazio in memoria per velocizzare
        x(1)=b(1)/A(1,1);
        for i=2:n
            if (isequal(A(i,i),0)) 
                error ('Elemento diagonale matrice nullo')
            end
            %A(i,i-1)=beta(i)
            %Algoritmo generale quindi divido tutto per L(i,i), applicabile
            %non solo a Thomas
            x(i)=(b(i)-A(i,i-1)*x(i-1))/A(i,i);
        end
    else
        error('La matrice non è bidiagonale inferiore')
    end
end

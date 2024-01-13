function[x,iter,res]=jacobi(A,b,toll,indx,varargin)
    list={'Relative residue control (euclidean norm)',...
        'Increment control (infinity norm)'};
    tf=0;
    if (indx~=1&&indx~=2)
        while tf==0
            [indx,tf]=listdlg('PromptString',{'Select control type'},...
                'ListString',list);
        end
    end
    [n,m]=size(A);
    if (n~=m)
        error('Error: the system has too many or too few arguments')
    end
    %%Alternative version where tolerance is a user input
    %prompt='Insert tolerance;
    %toll=input(prompt);
    %toll=1e-7;
    x=zeros(n,1);
    xprec=x;
    iter=1;
    %Selecting initial residue
    switch indx
        case 1
             res(iter)=norm(b-A*x);
             r0=res(iter);
        case 2
            %Initial value higher than toll
            res(iter)=10;
    end
    while res(iter)>toll
        xprec=x;
        iter=iter+1;  %I want to print how many iterations it took to get to the result
        for i=1:n
            s=0;
                for j=1:n
                    if (j~=i) s=s+A(i,j)*xprec(j); end 
            x(i)=(b(i)-s)/A(i,i);
                end
        end
        switch indx
            case 1
                %Relative residue control
                res(iter)=norm(A*x-b)/r0;
            case 2
                %Increment control
                res(iter)=norm(x-xprec,Inf);
        end
    end
    %fprintf('Iterations=  %d',iter);
end

%%function implementing SOR (successive over relaxation) and Gauss-Seidel when omega=1

function[x,iter,res]=SOR(A,b,omega,toll,indx,varargin)

      list={'Controllo residuo relativo norma euclidea',...
        'Controllo incremento norma infinito'};
    tf=0;

    if (nargin~=5)
        while tf==0
            [indx,tf]=listdlg('PromptString',{'Selezionare tipo di controllo'},...
                'ListString',list);
        end
    end
    [n,m]=size(A);
    if (n~=m)
        error('Sistema sovra o sotto dimensionato')
    end
    %%Memory allocation
    %% Alternative version where tolerance is a user input
    %   prompt='Insert tolerance value';
    %   toll=input(prompt);
    %   prompt='insert relaxation parameter value'
    %   omega=input(prompt)
    %toll=1e-7;
    x=zeros(n,1);
    xprec=x;
    iter=1;
    %res has different values depending on user choice
    switch indx
        case 1
             res=b-A*x;
             r0=norm(res);
        case 2
            %Valore iniziale più alto di toll
            res=Inf;
    end
   %algorithm start
   while res(iter)>toll
       xprec=x;       
       iter=iter+1;
       for i=1:n
           s=0;
           for j=1:n
               if (j~=i) s=s+A(i,j)*x(j); end
           end
           x(i)=(omega*(b(i)-s))/A(i,i)+(1-omega)*xprec(i);
       end
       switch   indx
           case 1
               res(iter)=norm(A*x-b)/r0;
           case 2
               res(iter)=norm(x-xprec,Inf);
       end

   end
end

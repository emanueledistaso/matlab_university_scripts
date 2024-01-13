function[u]=eulesp(t0,h,tf,u0,f)
    %f function handle di variabili t,y tc dy/dt=f
    %t0 tf estremi temporali, h passo
    %u0=y0 è il valore al tempo 0 della funzione
    y=u0;
    u=y;
    for t=t0:h:tf-h
        fn=f(t,y);
        y=y+h*fn;
        u=[u;y];
    end
end 

 function[u]=eulmod(t0,h,tf,u0,f)
    %f function handle di variabili t,y tc dy/dt=f
    %t0 tf estremi temporali, h passo
    %u0=y0 è il valore al tempo 0 della funzione
    y=u0;
    u=y;
    for t=t0:h:tf-h
        fn=f(t+h/2,y+h*f(t,y)*0.5);
        y=y+h*fn;
        u=[u;y];
    end
    end

    function[u]=rk4(t0,h,tf,u0,f)
    %Implementazione metodo Runge-kutta a 4 stadi per risoluzione EDO
    %f function handle di variabili t,y tc dy/dt=f
    %t0 tf estremi temporali, h passo
    %u0=y0 è il valore al tempo 0 della funzione
    y=u0;
    u=y;
    for t=t0:h:tf-h
        K1=f(t,y);
        K2=f(t+h/2,y+h*K1/2);
        K3=f(t+h/2,y+h*K2/2);
        K4=f(t+h,y+h*K3);
        y=y+h*(K1+2*K2+2*K3+K4)/6;%y è uno scalare che salva un
        u=[u;y];
    end
    end

    function[u]=cranknic(t0,h,tf,u0,f,tol,maxiter)
    %f function handle di variabili t,y tc dy/dt=f
    %t0 tf estremi temporali, h passo
    %u0=y0 è il valore al tempo 0 della funzione
    %tol=tolleranza iterazione, maxiter=iterazioni massime
    %u vettore riga
    t=t0:h:tf;
    n=length(t);
    u_old=zeros(n,1);
    err=Inf;
    u=u_old;
    u_old(1)=u0; u(1)=u0; iter=0;
    while (err>=tol && iter<=maxiter)
        for i=1:n-1
             fnew=f(t(i),u(i))+f(t(i+1),u(i+1));
             u_old(i+1)=u(i+1);
             u(i+1)=u(i)+h*fnew*0.5; 
        end
        iter=iter+1;
        err=max((abs(u-u_old)),[],'all');
    end
    end

    function[u]=heun(t0,h,tf,u0,f)
    %f function handle di variabili t,y tc dy/dt=f
    %t0 tf estremi temporali, h passo
    %u0=y0 è il valore al tempo 0 della funzione
    y=u0;
    u=y;
    for t=t0:h:tf-h
        fn=f(t,y)+f(t+h,y+h*f(t,y));
        y=y+h*fn*0.5;
        u=[u;y];
    end
    end

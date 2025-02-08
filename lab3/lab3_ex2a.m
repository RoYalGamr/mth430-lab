function [E, A] = lab3_ex2a(P, N, k, e)
   E = zeros(P*N+1,1);
   A = zeros(P*N+1,1);

   V_ini = [1-e, 0, 0 ,sqrt((1+e)/(1-e))];
    h = 2*pi/N;
    V = zeros(P*N+1,4);
    V(1,:) = V_ini;
   #  Adam Bashford 1 step method
    if (k == 1)
        for i = 1:N*P
            x = V(i,1);
            y = V(i,2);
            xdd = -x/((x^2 + y^2)^(3/2));
            ydd = -y/((x^2 + y^2)^(3/2));
            D = [V(i,3),V(i,4),xdd,ydd];
            V(i+1,:) = V(i,:) + h*D;
        end
    #  Adam Bashford 2 step method
    elseif (k==2)
        D = [V(1,3),   V(1,4),  -V(1,1)/((V(1,1)^2 + V(1,2)^2)^(3/2)),   -V(1,1)/((V(1,1)^2 + V(1,2)^2)^(3/2))];
        V(2,:) = V(1,:) + h*D;

        for i = 2:N*P
            D1 = [V(i-1,3),  V(i-1,4),   -V(i-1,1)/((V(i-1,1)^2 + V(i-1,2)^2)^(3/2)) ,  -V(i-1,2)/((V(i-1,1)^2 + V(i-1,2)^2)^(3/2))] ;
            D2 = [V(i,3),  V(i,4),   -V(i,1)/((V(i,1)^2 + V(i,2)^2)^(3/2)) ,  -V(i,2)/((V(i,1)^2 + V(i,2)^2)^(3/2))] ;
            V(i+1,:) = V(i,:) + h*(3*D2/2 - D1/2);
        end
    #  Adam Bashford 3 step method
    elseif (k==3)
        D = [V(1,3),   V(1,4),  -V(1,1)/((V(1,1)^2 + V(1,2)^2)^(3/2)),   -V(1,1)/((V(1,1)^2 + V(1,2)^2)^(3/2))];
        V(2,:) = V(1,:) + h*D;
        i = 2;
        D1 = [V(i-1,3),  V(i-1,4),   -V(i-1,1)/((V(i-1,1)^2 + V(i-1,2)^2)^(3/2)) ,  -V(i-1,2)/((V(i-1,1)^2 + V(i-1,2)^2)^(3/2))] ;
        D2 = [V(i,3),  V(i,4),   -V(i,1)/((V(i,1)^2 + V(i,2)^2)^(3/2)) ,  -V(i,2)/((V(i,1)^2 + V(i,2)^2)^(3/2))] ;
        V(i+1,:) = V(i,:) + h*(3*D2/2 - D1/2);

        for i = 3:N*P
            D1 = [V(i-2,3),  V(i-2,4),   -V(i-2,1)/((V(i-2,1)^2 + V(i-2,2)^2)^(3/2)) ,  -V(i-2,2)/((V(i-2,1)^2 + V(i-2,2)^2)^(3/2))] ;
            D2 = [V(i-1,3),  V(i-1,4),   -V(i-1,1)/((V(i-1,1)^2 + V(i-1,2)^2)^(3/2)) ,  -V(i-1,2)/((V(i-1,1)^2 + V(i-1,2)^2)^(3/2))] ;
            D3 = [V(i,3),  V(i,4),   -V(i,1)/((V(i,1)^2 + V(i,2)^2)^(3/2)) ,  -V(i,2)/((V(i,1)^2 + V(i,2)^2)^(3/2))] ;

            V(i+1,:) = V(i,:) + h*(23*D3/12 - 4*D2/3 + 5*D1/12);
        end

    end

    for i = 1:N*P+1
        E(i) = (V(i,3)^2 + V(i,4)^2)/2;
        E(i) = E(i) - 1/(((V(i,1)^2 + V(i,2)^2)^(1/2)));
        A(i) = V(i,1)*V(i,4) - V(i,2)*V(i,3);
    end

    t = [1:P*N+1];

    figure;
    plot(t,V(:,1),"r");
    xlabel("time");
    ylabel("X Location");

    figure;
    plot(t,V(:,2),"o");
    xlabel("time");
    ylabel("Y Location");

    figure;
    plot(V(:,1),V(:,2),"b");
    xlabel("X Location");
    ylabel("Y Location");

    figure;
    plot(t,E,"-");
    xlabel("time");
    ylabel("Energy");

    figure;
    plot(t,A,"o");
    xlabel("time");
    ylabel("Angular Momentum");


end

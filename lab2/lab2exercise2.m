clc;
clearvars;
clear;

function [] = lab2exercise2(N, T , a1, a2, b1, b2, y10, y20)
    E = zeros(N+1,2);
    E(1,1) = y10;
    E(1,2) = y20;

    for n = 1:N
        t = n*T/N;
        h = T/N;
        Enw = zeros(2,2);

#        Generating Initial value by euler method
        Enw(1,1) = E(n,1) + h*(E(n,1)*(a1 - b1*E(n,2)));
        Enw(1,2) = E(n,2) + h*(E(n,2)*(-a2 +b2*E(n,1)));

        a = E(n,1);
        b = E(n,2);
        while 1
            printf("Goinng...\n");
            J = [ (1-h*a1+h*b1*Enw(1,2)),  h*b1*Enw(1,1) ; -h*b2*Enw(1,2), 1+h*a2-h*b2*Enw(1,1) ];
            F = [ Enw(1,1)-h*a1*Enw(1,1)+h*b1*Enw(1,2)*Enw(1,1)-a  ; Enw(1,2)+h*a2*Enw(1,2)-h*b2*Enw(1,2)*Enw(1,1)-b ];
            Enw(2,:) = Enw(1,:) - (J \ F)';


            % Compute the difference without using norm
            diff = Enw(2,:) - Enw(1,:);
            diff_norm = sqrt(diff(1)^2 + diff(2)^2);  % Equivalent to norm(diff)

            if ( diff_norm < 10^-6)
                break;
            end
            Enw(1,:) = Enw(2,:);
        end
        E(n+1,:) = Enw(2,:);
    end
    disp(E);

    figure;
    plot(1:N+1,E(:,1),"r",1:N+1,E(:,2),"b");
    xlabel("time");
    ylabel("population");

    figure;
    plot3(1:N+1,E(:,1),E(:,2));
    xlabel("time");
    ylabel("prey");
    zlabel("predator");
    rotate3d;
    grid on;
end


lab2exercise2(4000,250,1,0.5,0.1,0.02,100,10);

function A = matrix_A(k12, k13, k23, k3)
    A = [-(k12 + k13)+1, 0, 0 ;
        k12, -k23+1, 0;
        k13, k23, -k3+1;];
end


function [Mo, r] = isobs(A, C)
    n = size(A, 1);
    Mo = [];
    for i = 0:n-1
        Mo = [Mo; C*A^i];
    end

    Mo = Mo';

    r = rank(Mo);
end

function Mr = isctr(A, B)
    n = size(A, 1);
    Mr = [];
    for i = 0:n-1
        Mr = [Mr, A^i*B];
    end  
end


function n_obtained = PBH_obs_test(A, C, n_needed)
    lambdas = eig(A);
    for i = 1:length(lambdas)
        M = [lambdas(i)*eye(n_needed)-A ;
            C];
        n_obtained = rank(M);
        if n_obtained ~= n_needed
            fprintf('Unobservable mode at lambda = %.4f\n', lambdas(i));
        end
    end
end


function n_obtained = PBH_rch_test(A, B, n_needed)
    lambdas = eig(A);
    for i = 1:length(lambdas)
        M = [lambdas(i)*eye(n_needed)-A B];
        n_obtained = rank(M);
        if n_obtained ~= n_needed
            fprintf('Unreachable mode at lambda = %.4f\n', lambdas(i));
        end
    end
end


function [V, Z, iTr, Tr, hA, hB] = T_ctrl(Mr, A, B)
    V = orth(Mr);
    Z = null(Mr');

    iTr = [V, Z];
    Tr = inv(iTr);

    hA = Tr * A * iTr;
    hB = Tr*B;

end



function [V, Z, iTo, To, hA, hC] = T_obs(Mo, A, C)
    V = null(Mo');
    Z = orth(Mo);

    iTo = [Z, V];
    To = inv(iTo);

    hA = To * A * iTo;
    hC = C*iTo;

end

function rk = r(M, name)
    rk = rank(M);
    fprintf('rank of %s = %d\n', name, rk);
end

function unobs_eig(A, C)
    lambdas = eig(A);
    for i = 1:length(lambdas)
        M = lambdas(i)*eye(size(A, 1)) - A;
        Obs = [M ; C];
        if rank(Obs)< size(A, 1)
            fprintf('Unobservable for lambdas = %.4f\n', lambdas(i));
        end
    end
end
h = 1;

k12 = 0.5/h;
k13 = 0.5/h;
k23 = 0.5/h;
k3  = 0.5/h;

Vol = 2;

A = matrix_A(k12, k13, k23, k3);

B = [0; 1; 0];

Mr = isctr(A, B);
r(Mr, 'Mr')


[V, Z, iTr, Tr, hA, hB] = T_ctrl(Mr, A, B);
V, Z, iTr, Tr, hA, hB

B1=[1 0;0  1;0  0];

Mr = isctr(A, B1);
r(Mr, 'Mr')

x3 = [0.5 1 1]';

u = Mr1\x3



k23 = 0;

A = matrix_A(k12, k13, k23, k3);

C = [0 0 1/Vol];

Mo = isobs(A, C);
r(Mo, 'Mo {k23 = 0}');

k23 = 0.5;
A = matrix_A(k12, k13, k23, k3);
Mo = isobs(A, C);
r(Mo, 'Mo {k23 = 0.5}');
unobs_eig(A, C);


[V, Z, iTo, To, hA, hC] = T_obs(Mo, A, C);
V, Z, iTo, To, hA, hC


X_no = V;
M = [Mo' ; C*A^3];
y = M*X_no;
plot(0:3, y);



k23 = 1;
A = matrix_A(k12, k13, k23, k3);

Mo = isobs(A, C);

y = [0; 0.75; 0.625]

Mo_T = Mo';
x0 = Mo_T\y
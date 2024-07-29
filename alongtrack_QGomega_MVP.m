function w=alongtrack_QGomega_MVP(u,b,N2,f,dx,dz,max_iterations,tolerance)
f2=f.^2;
ux=u2rho_2d((u(:,2:end)-u(:,1:end-1))./dx);
bx=u2rho_2d((b(:,2:end)-b(:,1:end-1))./dx);
Q=ux.*bx;
Qx=u2rho_2d((Q(:,2:end)-Q(:,1:end-1))./dx);


r1=N2./(dx.^2);
r2=f2./(dz.^2);




    % max_iterations: 
    % tolerance: 
    % omega: SOR方法的松弛参数 就随便取1.5

    % param

    [nz, nx] = size(Qx);
    

    %bc.rigid lid
    w = 0.*Qx;
%     divQ= 0.*w;

    % error
    err = inf;
    iteration = 0;

    % SOR intera
    while err > tolerance && iteration < max_iterations
        iteration = iteration + 1;
        w_old = w;

        for i = 2:nx-1
            for k = 2:nz-1
                w(k,i)=1./(r1(k,i)+r2).*Qx(k,i)+r1(k,i)./(2*(r1(k,i)+r2)).*(w(k,i+1)+w(k,i-1))+r2./(2*(r1(k,i)+r2)).*(w(k+1,i)+w(k-1,i));
            end
        end

        % error
        err = max(abs(w(:) - w_old(:)));
        error(iteration)=err;
        if mod(iteration,100)==0
            disp(num2str(iteration))
        end
    end

    fprintf('finish，inter:%d error：%e\n', iteration, err);
return


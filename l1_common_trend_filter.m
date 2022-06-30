function [ahat,xhat,lambda]=l1_common_trend_filter(Y,p,R,c)

tol=10^(-6); N=100; [T,n]=size(Y); 
I=eye(T); D=diff(I,2); Pi=[ones(T,1),(1:T)'];
xhat=Pi*inv(Pi'*Pi)*Pi'*Y(:,1);
for i=1:N
   atilde=Y'*xhat/(xhat'*xhat);
   ahat=atilde-R'*inv(R*R')*(R*atilde-c);
   Yast=Y*ahat/(ahat'*ahat);
   if i==1
      psi=(ahat'*ahat)*(2*sin(pi/p))^(-4);
      xhat_hp=(I+(psi/(ahat'*ahat))*D'*D)\Yast;
      cvx_begin
         variables xhat(T)
         minimize(norm(D*xhat,1))
         subject to
            sum((Yast-xhat).^2)<=sum((Yast-xhat_hp).^2)
      cvx_end
      lambda=2*(ahat'*ahat)*norm((D*D')\(D*(Yast-xhat)),inf);
   else
      lam=lambda/(ahat'*ahat);
      cvx_begin
         variables xhat(T)
         minimize(sum_square(Yast-xhat)+lam*norm(D*xhat,1))
      cvx_end
   end
   optval=norm(Y-xhat*ahat','fro')^2+lambda*norm(D*xhat,1);
   if i==1
      optval_old=optval;
   else
      if abs(optval-optval_old)<tol
         break
      else
         optval_old=optval;
      end
   end
end

end

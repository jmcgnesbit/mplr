function f = irfpdf(Atilde,const,EvecB,NT,nuT,ST,vecB)

% Purpose:
% This code computes the log of the joint posterior density of structural impulse responses up to constant
% input :
% Atilde: Obtained from rotating matrix A by matrix U, i.e., Atilde=A*U
% const : 1 if the intercept term is included; 0 otherwise
% EvecB : Posterior mean of vec(B)
% NT    : N0+X'*X
% nuT   : degrees of freedom for the inverse-Wishart distribution
% ST    : Scale matrix for the inverse-Wishart distribution
% vecB  : Posterior draw of vec(B) = vec([c Phi1 ... Phip]')
% output:
% y     : Log of the joint posterior density up to constant
%
% Record of revisions:
% Date        Programmer            Description of change
% 02/15/2011  Atsushi Inoue         Original code
% 05/11/2011  Atsushi Inoue         Modified
% 02/26/2014  Jonas Arias' corrections incorporated


n = size(Atilde,1);           % The number of variables
p = (size(vecB,1)/n-const)/n; % The order of the VAR model
if const==1
    index = kron(ones(n,1),[0;ones(n*p,1)]);
else
    index = ones(n*n*p,1);
end

% Communication Matrix
Kn=[];
for i=1:n
    for j=1:n
        E      = zeros(n,n);
        E(i,j) = 1;
        Kn     = [Kn reshape(E,n*n,1)];
    end
end
n2 = n^2;
Kn2=[];
for i=1:n2
    for j=1:n2
        E      = zeros(n2,n2);
        E(i,j) = 1;
        Kn2     = [Kn2 reshape(E,n2*n2,1)];
    end
end

% Duplication Matrix
Dn = [];
for j=1:n
    for i=j:n
        E      = zeros(n,n);
        E(i,j) = 1;
        E(j,i) = 1;
        Dn     = [Dn reshape(E,n*n,1)];
    end
end
Dnplus = inv(Dn'*Dn)*Dn';

% Duplication Matrix" for Non-Symmetric Matrices (such as A)
Dbarn = [];
for j=1:n
    for i=j:n
        E      = zeros(n,n);
        E(i,j) = 1;
        Dbarn     = [Dbarn reshape(E,n*n,1)];
    end
end

% Eh and Ek Matrices
Eh = [];
for j=1:n
    for i=j:n
        E      = zeros(n,n);
        E(i,j) = 1;
        Eh     = [Eh reshape(E,n*n,1)];
    end
end
Ek = [];
for j=2:n
    for i=1:j-1
        E      = zeros(n,n);
        E(i,j) = 1;
        Ek     = [Ek reshape(E,n*n,1)];
    end
end

% Compute reduced-form impulse responses
B     = (reshape(vecB,const+n*p,n))';
Phi   = [B(:,1+const:end);eye(n*(p-1)) zeros(n*(p-1),n)];
Theta = eye(n);
for i = 1:p
    Phii  = Phi^i;
    Theta = [Theta Phii(1:n,1:n)];
end

% Compute the Jacobian of vec(Theta) with respect to vech(A), veck(U) and vec(Phi)
Sigma       = Atilde*Atilde';
A           = chol(Sigma)';
U           = inv(A)*Atilde;
%Ju          = Ek-Eh*pinv(Dnplus*(kron(U,eye(n))+kron(eye(n),U)*Kn)*Eh)*Dnplus*(kron(U,eye(n))+kron(eye(n),U)*Kn)*Ek;
%J1          = [kron(eye(n),A)*Ju kron(U',eye(n))*Dbarn];
%J3          = Dnplus*(kron(A,eye(n))+kron(eye(n),A)*Kn)*Dbarn;
J           = Eh'*(kron(A,eye(n))+kron(eye(n),A)*Kn)*Eh;

%logdetJ     = -log(abs(det(J1)))-0.5*n*p*log(det(Sigma))+log(abs(det(J3))); % The second term is the log-determinant of the Jacobian w.r.t. B
logdetJ     = log(abs(det(J)))+0.5*n*p*log(det(Sigma)); % The second term is the log-determinant of the Jacobian w.r.t. B
VvecB       = kron(Sigma,inv(NT));
VvecBinv    = inv(VvecB(index==1,index==1));
logdetVvecB = logdet(VvecB(index==1,index==1));
logdetVvecB = n*p*log(det(Sigma)); % Alternative if the previous line is numericaly unstable
f           = logdetJ-0.5*logdetVvecB-0.5*(vecB(index==1,1)-EvecB(index==1,1))'*VvecBinv*(vecB(index==1,1)-EvecB(index==1,1));
f           = f-0.5*(nuT+n+1)*log(det(Sigma))-0.5*trace(nuT*ST*inv(Sigma));



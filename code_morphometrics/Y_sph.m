function [ Y ] = Y_sph(I,N,theta,phi)
% Compute the spherical harmonics
Y=zeros(I,(N+1)^2);
for i=1:I
    k=1;
    for n=0:N
        P=legendre(n,cos(theta(i))); 
        for m=-n:n  
            if m>0
                Y(i,k)=(sqrt(((2*n+1)*factorial(n-abs(m)))/((4*pi)*factorial(n+abs(m)))))*P(abs(m)+1)*sqrt(2)*cos(m*phi(i));

            elseif m == 0
                Y(i,k)=(sqrt(((2*n+1)*factorial(n))/((4*pi)*factorial(n))))*P(1);
                
            else
                Y(i,k)=(sqrt(((2*n+1)*factorial(n-abs(m)))/((4*pi)*factorial(n+abs(m)))))*P(abs(m)+1)*sqrt(2)*sin(abs(m)*phi(i));

            end
            k=k+1;
            
        end
    end
end
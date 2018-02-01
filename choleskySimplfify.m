function [G]=choleskySimplify(matDimM,matDimN,dx,w,lambda,velocity)
	%Obtención de C
	C=exp(-(lambda*w*dx)/(2*pi*velocity));
	G=zeros(matDimM,matDimN);
	%Obtención de Sw y G basado en https://www.researchgate.net/publication/245286240_Simulation_of_Stochastic_Wind_Velocity_Field_on_Long-Span_Bridges?enrichId=rgreq-258f8163b360e69a5eee767aa7941079-XXX&enrichSource=Y292ZXJQYWdlOzI0NTI4NjI0MDtBUzozODcxMzk0NzM4MjE2OTZAMTQ2OTMxMjY1MTIyMw%3D%3D&el=1_x_2&_esc=publicationCoverPdf
	for j=1:matDimM
		for m=1:j
			if m==1 && m<=j
				G(j,m)=C^(abs(j-m));
            else
				G(j,m)=C^(abs(j-m))*sqrt(1-C^2);
			end
		end
	end
end
%%Anexar control de errores cuando m>n

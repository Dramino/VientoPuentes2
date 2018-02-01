function [x]=WindSimMult(Sw,numSim,dx)
	tic
	numFreq=2;
	Dw=0.4;
	Coh=0;
	t0=0;
	tf=600;
	Ndatos=61;
	t=linspace(t0,tf,Ndatos);
	w0=0.4;
    w=w0;
	x=zeros(numSim,length(t));
	fi=[4.3982 3.7699 1.8850; 0.6283 2.5132 1.2566; 0.6283 4.3982 5.0265];
	for i=1:numFreq
		%Coh=exp(-(10*w*dx)/(2*pi*velocity));
		%Cholesky decomposition
		C=[0.6376 0.4066 0.1541];
		G=choleskySimplify(numSim,C(i));
		a=sqrt(2*Sw(i)*Dw);
		
		for j=1:numSim
			for k=1:numSim
				x(k,:)=x(k,:)+G(k,j).*a.*cos(w.*t+fi(j,i));
			end
        end
        w=w+Dw;
	end
	toc
end

function [G]=choleskySimplify(matDimM,Coh)
	%Obtención de C
	C=Coh;
	G=zeros(matDimM);
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
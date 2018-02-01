function [x]=WindSimMult(Sw,numSim,dx)
	tic
	numFreq  = 2;
	Dw       = 0.4;
	Coh      = 0;
	t0       = 0;
	tf       = 600;
	Ndatos   = 61;
	velocity = 40
	t        = linspace(t0,tf,Ndatos);
	w0       = 0.4;
    w        = w0;
	x        = zeros(numSim,length(t));
	fi       = 2 * pi * rand(numSim,numFreq);
	for i = 1:numFreq
		C = exp( -(10 * w * dx) / (2 * pi * velocity));
		%Cholesky decomposition
		G = choleskySimplify(numSim,C(i));
		
		a = sqrt(2 * Sw(i) * Dw);
		
		for j=1:numSim
			for k=1:numSim
				x(k,:) = x(k,:) + G(k,j) .* a .* cos( w .* t + fi(j,i));
			end
        end
        w = w + Dw;
	end
	toc
end

function [G]=choleskySimplify(matDimM,Coh)
	%Obtención de C
	C = Coh;
	G = zeros(matDimM);
	for j=1:matDimM
		for m = 1:j
			if m == 1 && m<=j
				G(j,m) = C ^ (abs(j-m));
            else
				G(j,m) = C ^ (abs(j-m)) * sqrt(1-C ^ 2);
			end
		end
	end
end

%%Nota: Esta versión funciona para matlab R2016 en adelante

function [re_val]=kw(nu,mat_A,N_s,t)
 if mod (nu ,t* 3)<0.1
     re_val =0; % return zero for harmonic of order n*3
 else
      X = 0;
      Y = 0;
      R = 0;
      N = 0;
      ase = 2 * pi / N_s * nu;
      val=0;
      for i=1:N_s
        val = mat_A(i);
        angle = ase * i;
        if(val<0)
            val = -val;
            angle = angle+pi;
        end
        X = X + val * cos(angle);
        Y = Y + val * sin(angle);
        N = N + val;
      end
      R = sqrt(X*X+Y*Y);
      re_val = R/N;
      if(re_val<1e-8)
          re_val = 0;
      end
 end

    
function output= halfrootfilter(input,beta,fs,kind)

switch kind,
    case 'trans'
        for i=1:length(input)
            if input(i)==0
                output(i)=fs*(1+beta*(4/pi -1));
            elseif input(i) == 1/(fs*4*beta) | input(i) == -1/(fs*4*beta)
               output(i)= (fs/sqrt(2))*beta*((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/4*beta));
            else
                output(i)=fs*(sin(pi*(1-beta)*input(i)*fs)+4*beta*fs*input(i)*cos(pi*fs*(1+beta)))/(pi*fs*(1-(4*beta*input(i)*fs)^2));
            end
        end
        
    case 'res'
          for i=1:length(input)
            if input(i)==0
                output(i)=fs*(1+beta*(4/pi -1));
            elseif input(i) == 1/(fs*4*beta) | input(i) == -1/(fs*4*beta)
               output(i)= (fs/sqrt(2))*beta*((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/4*beta));
            else
                output(i)=fs*(sin(pi*(1-beta)*(-1)*input(i)*fs)+4*beta*fs*(-1)*input(i)*cos(pi*fs*(1+beta)))/(pi*fs*(1-(4*beta*(-1)*input(i)*fs)^2));
            end
          end
end
%FFTn_coulomb(dr, WF1, WF2, option)
%
%option = 'potential' or 'integral' or 'exchange potential' or 'exchangeintegral' 
%
%dr in nm
%WF1 normalized so sum(sum(abs(WF1)^2)) = 1
%
%Output in meV
%
%
%Note: dr = dx=dy=dz

function [out1, out2] = FFTn_coulomb(dr, WF1, WF2, option)
WF1 = WF1/sqrt(sum(sum(sum(WF1.^2))));
WF2 = WF2/sqrt(sum(sum(sum(WF2.^2))));

    %persistent H;
    H = [];
    
    epsilon_r = 13.214;
    e = 1.60219E-19;
    c = 2.997925E8;
    epsilon = 1E7/(4*pi*c^2);
    eV = e;
    meV = eV/1000;
    
    exchange = 0;
    
    if nargin<3
        WF2 = WF1;
    end

    if nargin<4
        option = '';
    end
    
    if length(option>7)
        if strcmp(option(1:8), 'exchange');
            exchange = 1;
            option(1:8) = '';
        end
    end
    
    [sy, sx, sz] = size(WF1);
    dim = length(find([sy sx sz]>1));
            
    if ((size(H, 1)~=size(WF1, 1))|(size(H, 2)~=size(WF1, 2))|(size(H, 3)~=size(WF1, 3))|(nargout < 1))
        x = (1:sx)*dr; x = x - mean(x);
        y = (1:sy)*dr; y = y - mean(y);
        z = (1:sz)*dr; z = z - mean(z);
        [X, Y, Z] = meshgrid(x , y, z); R = sqrt(X.^2+Y.^2+Z.^2);
        %disp('update H...')
        if dim == 2
            R(find(R==0)) = 1/(4/dr);
        elseif dim == 3
            R(find(R==0)) = 1/(3/dr);
        end
        invR = 1./R;
        H = fftn(invR);
                    
    if (nargout >= 1)

        CC = e^2/(4*pi*epsilon*epsilon_r*1E-9)/meV;
        
        if exchange
            F = fftn(conj(WF2).*WF1);
        else
            F = fftn(conj(WF1).*WF1);
        end
        
        C = CC*ifftn(F.*abs(H));
        
        if (nargout==2)|strcmp(option, 'integral')
            
            if exchange
                I = sum(sum(sum(conj(WF1).*WF2.*C)));
            else
                I = sum(sum(sum(conj(WF2).*WF2.*C)));
            end
            
        end

        if nargout==2
            out1 = C;
            out2 = I;
        elseif strcmp(option, 'integral')
            out1 = I;
        else
            out1 = C;
        end
        
    end
    
end

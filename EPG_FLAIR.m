%F0 is complex as EPG stores transverse magnetisation as complex numbers
%and imaginary part contaons phase info and imaginary part is essentially 0
function [F0,Fn,Zn,F] = EPG_FLAIR(theta,ESP,T1,T2,TI,varargin);

for ii=1:length(varargin)
    if strcmpi(varargin{ii},'kmax')  
        kmax = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'zinit')
        zinit=varargin{ii+1};
    end
    
end


np = length(theta); 
if ~exist('kmax','var') 
    kmax = 2*(np - 1); 
end 

if isinf(kmax) 
    allpathways = true; 
    kmax = 2*(np - 1);
else
    allpathways = false; 
end


if allpathways
    kmax_per_pulse = 2*(1:np); 
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax; 
else
    kmax_per_pulse = 2*[1:ceil(np/2) (floor(np/2)):-1:1]+1; 
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax; 
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end


N=3*(kmax+1); 
alpha = abs(theta); 
phi=angle(theta); 
phi(2:end) = phi(2:end) + pi/2; 
S = EPG_shift_matrices(kmax); 
S = sparse(S); 
E1 = exp(-0.5*ESP/T1); 
E2 = exp(-0.5*ESP/T2);
E = diag([E2 E2 E1]); 
b = zeros([N 1]); 
b(3) = 1-E1; 


if exist('diff','var') 
    E = E_diff(E,diff,kmax,N);
else
    E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N); 
end
    
SE=S*E; 
SE=sparse(SE); 
T = sparse(zeros(N,N));
i1 = []; 
for ii=1:3  
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3))); 
end


F = zeros([N np-1]); 
FF = zeros([N 1]); 
if exist('zinit','var') 
    FF(3)=zinit; 
else
    FF(3)=1; 
end
   

%FLAIR inversion recovery

    function FF = FLAIR_Prep( FF, TI, T1, T2)

    Ainv = RF_rot(pi,0); % creates RF rotation matrix for 180 instead of 90 and phase 0
    FF(1:3) = Ainv * FF(1:3); % before inversion F+, F-, Z0 and after F+, F- and -Z0

    E1_TI = exp(-TI/T1); 
    E2_TI = exp(-TI/T2); % E1 and E2 calculate relaxation factors
    FF(1) = E2_TI * FF(1); % transverse magnetisation disappears during inversion recovery
    FF(2) = E2_TI * FF(2); % same for F-
    FF(3) = E1_TI * FF(3) + (1-E1_TI); % Bloch equation written in EPG variables

    end

if TI > 0
    FF = FLAIR_Prep( FF, TI, T1, T2);
end

%%%%% 

A = RF_rot(alpha(1),phi(1)); 
FF(1:3) = A*FF(1:3); 
kidx=1:6; 
FF(kidx) = SE(kidx,kidx)*FF(kidx)+b(kidx);

for jj=2:np 
    A = RF_rot(alpha(jj),phi(jj)); 
    build_T(A);
    kidx = 1:3*kmax_per_pulse(jj); 
    FF(kidx)=T(kidx,kidx)*FF(kidx); 
    F(kidx,jj-1) = SE(kidx,kidx)*FF(kidx)+b(kidx); 
    F(1,jj-1)=conj(F(1,jj-1)); 
    
    if jj==np
        break 
    end
    FF(kidx) = SE(kidx,kidx)*F(kidx,jj-1)+b(kidx); 
    FF(1)=conj(FF(1));  
    
end

F0 = F(1,:); %remove *1i
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)]; 
kvals = -kmax:kmax; 
idx(1:2)=[];
kvals(1:2)=[]; 


Fn = F(idx,:); 
Fn(kvals<0,:)=conj(Fn(kvals<0,:)); 
Zn = F(3:3:end,:); 


    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

    function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end




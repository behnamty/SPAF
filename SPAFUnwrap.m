%Smart phasor average filtering (SPAF)

% 1- SPAF is a filtering method to reduce
% Number of the phase residues of a wrapped phase.

% 2- The filtered image can be unwrapped using different unwrapping methods.
% A simple way is employing 1D MATLAB unwrapping function, but we found that
% combination of SPAF and 2D MNF unwrapping can be more robust to noise.

% 3- SPAF finds the phase residues and only filters a small patch around those  
% points with minimum filter size. In this case the image degradation will be 
% restricted to those small patches. SPAF will preserve both high-resolution
% and noise (except residue noise)information in the image. It is
% appropriate for microscopic or any high resolution applications.

% 4- variables:
% Phasor    = Complex amplitude, A*exp(i*phi); 
% maxFilter = maximum number of pixels for averaging, e.g. 6 means average of 6 neighbor pixels;
% maxPaf    = maximum patch size, the filter is employed only in this patch around a residue;

% 5- The code was written by Behnam Tayebi on August 2020 for 
% "Smart Filtering of Phase Residues in Noisy Wrapped Holograms",
% please send an email to behnamty@gmail.com for any questions.


% Main Function
function [SPAFphasor,unwrappedPhase]=SPAFUnwrap(phasor,maxFilter,maxPatch)

if nargin<2;maxFilter=floor(length(phasor)/50)+3;end           %filter size
if nargin< 3;maxPatch=floor(length(phasor)/50)+3;end           %patch size

    x=phasor./abs(phasor);                      %normalization
    patchi=1;

for ix=2:maxPatch
    for iy=2:maxPatch
         
        if sum(patchi(:))==0                    %break if there is no residue
           break
        end
        
        singlePatch=ones(ix,iy);                 %patch matrix 
        
        phi=angle(x);                           %phase
        
        resi=phaseResidues(phi)~=0;             %find location of residue
        
        patchi=conv2(resi,singlePatch,'same')>0;%build binary patch around residues

       for jy=2:maxFilter
          for jx=2:maxFilter

              x=conv2(x,ones(jx,jy)./jx./jy,'same').*patchi+...
                  (1-patchi).*x;                %local filtering

          end
       end
       
    end
end

    SPAFphasor=x;                               %final phasor without phase residue
    unwrappedPhase=matlabUnwrap(angle(x));      %final unwrapped phase by 1D unwrapping function in matlab

end



%%

% Function #S1  phase Residues Locations
function resi=phaseResidues(phase)
    
    % different directions
    down =circshift(phase,[-1 0]);  
    right =circshift(phase,[0 -1]); 
    downright =circshift(phase,[-1 -1]); 

    % Residues
    res1=mod(phase - down + pi, 2*pi);     
    res2=mod(down - downright + pi, 2*pi);  
    res3=mod(downright - right + pi, 2*pi) ; 
    res4=mod(right - phase + pi, 2*pi) ;     

	totResi=res1+res2+res3+res4-4*pi;           %Total Residues     

    resi=(totResi>=6) - (totResi<=-6);  
    resi(:,end)=0; resi(end,:)=0;        
end


% Function #S2  Matlab 1D phase unwrap
function uPhase=matlabUnwrap(x)

    phaseUnwrapX=unwrap((x));                   % Matlab 1D phase unwrap
    
    rotatePhaseUnwrapX=rot90(phaseUnwrapX,1);   % rotating the matrix 90 degree
    
    phaseUnwrapY=unwrap(rotatePhaseUnwrapX);    % Matlab 1D phase unwrap
    
    uPhase=rot90(phaseUnwrapY,1);               % rotating the matrix 90 degree
    
end



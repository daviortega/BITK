function [pwX_wf,p_wf] = project_aln(algn,wX,W)
% usage: [pwX]=project_aln(algn,wX,W); 

% This function computes the projected weighted alignment (pwX_wf) given
% the alignment in ascii format (algn, to compute the frequencies), the
% weighted 3-D alignment tensor (wX, generated by weight_aln.m), and the
% weight matrix (W, generated by weight_aln.m). 

% The projection vector (p_wf) is calculated as the weighted frequency of
% each amino acid normalized over that of all amino acids at each position
% (see Note xx). Where does this calculation come from? As discussed in
% Note xx, the projection vector conceptually derives from an interesting
% property of the singular value decomposition (SVD) of each 20 x 20 amino
% acid correlation matrix for each pair of amino acids i and j.  The
% findings are (1) that the top singular value sufficies to capture the
% information content of this matrix for each i and j, and (2) that the
% first left singular vectors for position i over all other positions j is
% nearly invariant.  That is, the weights for amino acids that contributes
% to correlations for a given position i is not strongly dependent on the
% identity of position j. Thus we can take the mean left singular vectors
% for each position i over all j as a "projection vector" that indicates
% the weights for each of the 20 amino acids at position i of relevance.
% This also implies that the projection should be computable from just
% knowledge of the amino acid composition at each position without regard
% to correlations at all.  Indeed, the calculation made here - the
% normalized weighted frequency of each amino acid at position i - produces
% a projection vector that is an excellent approximation to that derived
% from SVD's of each 20 X 20 amino acid correlation matrix for position i
% over all j (see Note xx).

% This function is called in sca5.m

% ***********************
% Authors: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%          Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 8/2011
%
% Copyright R.Ranganathan 1999-2011


[N_seq,N_pos,N_aa]=size(wX);

% first calculate frequencies from alignment; the function freq.m is
% embedded here.

freq_aln=freq(algn);

% We define a projection matrix from the weighted normalized frequency at
% each position

p_wf=zeros(N_pos,N_aa);
wf=zeros(N_aa,N_pos);

for i=1:N_pos;
    for a=1:N_aa
    wf(a,i)=(W(i,a)*freq_aln(a,i));
    end
    if norm(wf(:,i))>0; p_wf(i,:)=wf(:,i)./norm(wf(:,i));end
end

pwX_wf=zeros(N_seq,N_pos);

for i=1:N_pos
    for a=1:N_aa 
        pwX_wf(:,i)=pwX_wf(:,i)+p_wf(i,a)*wX(:,i,a); 
    end 
end

end

function [f]=freq(algn)
% usage: [freq_matrix]=freq(alignment)
%
% Amino acid frequencies for all positions of a multiple sequence
% alignments (algn).  Given an alignment of M sequences by N positions,
% returns a 20 X N positions matrix of aa frequencies.
%
% initalize variables and preliminaries

[N_seq,N_pos]=size(algn);
N_aa=20;
alg=lett2num(algn);
Abin=zeros(N_seq,N_pos,N_aa);
for a=1:N_aa
    Abin(:,:,a)=(alg==a);
end
f=zeros(N_aa,N_pos);

% determine frequencies
for aa = 1:N_aa
    f(aa,:) = sum(Abin(:,:,aa))./N_seq;
end


end

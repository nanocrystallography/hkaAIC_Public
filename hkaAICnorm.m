%hkaAIC but with normalized amplitudes
%Andrew H. Dempsey

clear;close all; clc; 

% Stores names of all 17 plane groups and 21 settings thereof following the
% ordering/numbering convention of CRISP
GroupNames = {'p1'; 'p2'; 'pm m_|_x'; 'pm m_|_y'; 'pg g_|_x';'pg g_|_y';...
              'cm m_|_x';'cm m_|_y'; 'p2mm'; 'p2mg m_|_x'; 'p2mg m_|_y';...
              'p2gg'; 'c2mm'; 'p4'; 'p4mm'; 'p4gm'; 'p3'; 'p3m1';'p31m';...
              'p6';'p6mm'}; 

% Stores Multiplicities of General Positions for use later
% the nth element of this vector corresponds to the multiplicity of the
% general position in the nth space group, following CRISP numbering
% convention for the 21 settings of the 17 2D space groups
MGP = [1 2 2 2 2 2 4 4 4 4 4 4 8 4 8 8 3 6 6 6 12]';

% Generates pop-up window to import hka file into workspace
% Collects path and file names to be used to import data
[FileName,PathName] = uigetfile({'*.hka',...
    'All Image Files';'*.*','All Files' },...
    'Import Image File','C:\','MultiSelect','on');

% Fixes issue where program would break if only one file was imported. In
% such cases, the file name is stored as a character vector instead of a
% cell array of strings. As such, the length is however many characters is
% in the file name. This checks to see if FileName is a cell 
% and if it is not, converts it to a cell so length is equal to one. 
if iscell(FileName) == 0
    FileName = {FileName};
end

% Pre-allocated array to store file IDs
fileID = zeros(1,length(FileName));

% Pre-allocated cell array to store all important information about hka
% files

C = {zeros(length(FileName),17)};

% Import pertinent data from hka files
for m = 1:length(FileName)



% Opens the file selected above 
fileID(m) = fopen(cell2mat(strcat(PathName,FileName(m))));

% Extracts and stores space group number (per CRISP convention) in column
% 17 of C
C(m,17) = textscan(fileID(m), '%*10s %f',1,'headerlines',1,'delimiter','\n');
frewind(fileID(m));

% Extracts and stores real-space lattice parameters in column 19 of C as 
% [a; b; gamma]
C(m,19) = textscan(fileID(m), '%*10s %f',3,'headerlines',2,'delimiter','\t');
frewind(fileID(m));

% Imports data from selected .hka file into a cell (C) whose contents 
% are C(m,1): h; C(m,2): k; C(m,3): Amp; C(m,4): AmpS; C(m,5): Pha; 
% C(m,6): PhaS; C(m,7): Err
% I've assumed all .hka files have 7 header lines. 
C(m,1:7) = textscan(fileID(m),'%f %f %f %f %f %f %s','HeaderLines',7, ... 
    'Delimiter','\t');

% Stores model name of imported file in 8th column of C
C(m,8) = {GroupNames(C{m,17})};

% Stores Multiplicity of General Position of given model in 9th column of C
C(m,9) = {MGP(C{m,17})};


end



%Stores file names in column 12 of C
C(:,12) = FileName';

% Sorts rows of C in ascending order of plane group number
C = sortrows(C,17);


% Get rid of "high resolution" peaks

    % Calculate number of pixels along diagonals of each unit cell. Divide
    % that length by 3 and round down to nearest integer (such that 0.99 = 0
    % instead of 1) to determine maximum h,k.
    for m = 1:size(C,1)

    Diagonals(m,1) = sqrt(C{m,19}(1)^2 + C{m,19}(2)^2 - ... 
                          2*C{m,19}(1)*C{m,19}(2)*sin((C{m,19}(3) - 90)*(pi/180)));
                  

    % Calculates maximum number of full oscillations of a sin function that 
    % can occur along the pixels lying on the diagonal, uses this integer as
    % the maximum values of h and k
    MaxRes(m) = floor(Diagonals(m,1)/3);

    Max_hk(m,1) = floor(C{m,19}(1)/3);

    Max_hk(m,2) = floor(C{m,19}(2)/3);

    
end
MinRes = min(min(Max_hk));
MaxRes = min(MaxRes);


% The goal of this loop is to ensure approximately the same number of peaks
% is used in calculations for each plane symmetry group

% for m = 1:size(C,1)
%     for n = 1:length(C{m,1})
%         
%         % Gets rid of peaks that: 
%         %   -fall outside a circular area of radius smaller than MaxRes OR
%         %    
%         %   -have observed amplitude less than 1% of maximum observed
%         %    amplitude (i.e. 100) OR
%         %
%         %   -have peaks with miller indices h k that are greater than the
%         %    max h k found in p2
%         if C{m,1}(n)^2+C{m,2}(n)^2 > (MaxRes)^2 ...
%                 || C{m,3}(n)<100 ... 
%                 || C{m,1}(n)>max(C{1,1}) ...
%                 || C{m,2}(n) > max(C{1,2})...        
%             C{m,3}(n) = 0;
%             C{m,4}(n) = inf;
%             C{m,5}(n) = 0;
%             C{m,6}(n) = 0;
%             C{m,10}(n) = 0;
%             C{m,11}(n) = 0;
%             
% 
%         else
%         end
%     
%          
%         % applies reflection conditions for pg, pmg, etc. 
%         if mod(C{m,1}(n),2) == 1 && C{m,2}(n) == 0 ...
%               || mod(C{m,2}(n),2) == 1 && C{m,1}(n) == 0 ...
%               || mod(C{m,1}(n)+C{m,2}(n),2) == 1 
%  
%             C{m,3}(n) = 0;
%             C{m,4}(n) = inf;
%             C{m,5}(n) = 0;
%             C{m,6}(n) = 0;
%             C{m,10}(n) = 0;
%             C{m,11}(n) = 0;
%         else
%         end
%    
%     
%         Minimum_Amplitude(m) = min(C{m,4}(C{m,4}>0));
%         
%             if C{m,4}(n) == inf
%                 C{m,4}(n) = 0;
%             else
%             end
%         
%     end
%        
% end



Minimum_Amplitude = Minimum_Amplitude';

NormC=cell(size(C,1),2);

% normalizes observed and symmetrized amplitudes by largest
% observed/symmetrized amplitude
for m = 1:size(C,1)
    NormC{m,1} = C{m,3}/max(C{m,3});
    NormC{m,2} = C{m,4}/max(C{m,3});
end

% Combine amplitudes and phases of Fourier coefficients, convert to
% rectangular form, and calculate RSS. 
for m = 1:size(C,1)
    % Raw Fourier coefficients
    [a,b] = pol2cart(C{m,5}*(pi/180),NormC{m,1});
    C{m,10} = a + 1i*b;
    
    % Symmetrized Fourier coefficients
    [c,d] = pol2cart(C{m,6}*(pi/180),NormC{m,2});
    C{m,11} = c + 1i*d;
    
    
    
end






for m = 1:size(C,1)
   % Stores number of peaks used for residual calculations
   C{m,18} = nnz(C{m,4}); 
   
   % calculates RSS for complex fourier coefficients.
 
   C{m,13} = sum((C{m,10}-C{m,11}).*conj(C{m,10}-C{m,11}));
   % Calculate RSS for amplitudes only. 

    C{m,14} = sum((NormC{m,1} - NormC{m,2}).^2);
    
   % Calculate R_sym for mth model and store result in column 14 of C
   C{m,15} = 100*sum(abs(abs(C{m,3}) - abs(C{m,4})))/sum(C{m,4});
   
   % Calculate phase residuals for mth model and store in column 15 of C 
   C{m,16} = sum(C{m,3}.*abs(abs(C{m,5}) - abs(C{m,6})))/sum(C{m,3});
end


    
   





% Calculate ratio of residuals 
ResRatio = zeros(size(C,1));
KRatio = zeros(size(C,1));
for m = 1:size(C,1)
    for n = 1:size(C,1)
        if C{m,9} > C{n,9}
            ResRatio(m,n) = C{m,13}/C{n,13};
            KRatio(m,n) = 1+(2*(C{m,9}-C{n,9})/(C{m,9}*(C{n,9}-1)));
        else
            ResRatio(m,n) = inf;
            KRatio(m,n) = inf;
        end
    end
end

% Determines if the inequality ResRatio < KRatio is true
% LogicalRatio = ResRatio < KRatio;

% Determines the best model in the set by checking each column and row of
% LogicalRatio and finding the ones with all zeros
% BestModelIdx = 0;
% for m = 1:size(LogicalRatio,2)
%     if nnz(LogicalRatio(:,m)) == 0 && nnz(LogicalRatio(m,:)) == 0
%         BestModelIdx = m;
%     elseif nnz(LogicalRatio(:,m)) == 0
%         BestModelIdx = m;
%     else
%             
%     end
% end

exist MaxRes;
if ans == 0;
    MaxRes = inf;
else
end


FC_RSS = cell2mat(C(:,13));
Amplitude_RSS = cell2mat(C(:,14));
R_sym = cell2mat(C(:,15));
Phi_Res = cell2mat(C(:,16));
Number_data_points = cell2mat(C(:,18));
% Minimum_Amplitude = MinAmp;
for m = 1:size(C,1)
    Enforced_Symmetry(m,1) = [C{m,8}];
end
File_Name = C(:,12);

% Generates table to assist with data inerpretation, and for possible
% export to Excel
T = table(Enforced_Symmetry, FC_RSS, Amplitude_RSS, R_sym, ...
    Phi_Res, Minimum_Amplitude, Number_data_points, File_Name)

% Construct a questdlg with two options, prompt user if they want to export
% to Excel
choice = questdlg('Would you like to export data to .xls?', ...
	'Export Data?', ...
	'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
        
        FilePrompt = {'Maximum Resolution: ', MaxRes};
        FilePrompt = sprintf('%s %d \n\nPlease enter file name:',FilePrompt{:});
        ExpFileName = inputdlg(FilePrompt);
        ExpFileName = char(strcat(ExpFileName,'.xls'));
        writetable(T,ExpFileName)
    case 'No'
        
        
    
end

residuals = FC_RSS;

for m = 1:size(C,1); normcoeff(m)=max(C{m,3});end



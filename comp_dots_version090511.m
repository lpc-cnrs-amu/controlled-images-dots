% This programme requires installation of the psychtoolbox:
% http://psychtoolbox.org/PTB-2/

% This programme generates stimuli for studies using non-symbolic number
% stimuli. The different visual properties of the stimuli are documented in
% an output file. Therefore post hoc analyses can be conducted to
% investigate whether a relation exists between numerical distance and the
% difference in visual properties of the stimuli or between the difference
% in visual properties of the stimuli and the numerical distance effect
% present in the behavioral or neuroimaging data.

% For most designs the program generates stimuli that are not confounded
% with visual cues. Nevertheless a post hoc analyses to verify whether this
% is indeed the case is recommended. Especially when small numerosities and
% large number distances are used, it is unavoidable that strong relations
% between number and area subtended or circumference arise.

% The programme generates (1) .bmp images of the stimuli, (2) a file called
% 'picturenames.txt' that contains all the names of the images. This file
% can be of help to load the bitmap images in E-Prime, (3) a file called
% 'all.txt' which contains all the information of each stimulus and (4) a
% file called 'regr.txt' which contains the information about the
% differences in visual and numerical properties of each stimulus pair. This
% file can be used to directly verify the relation between the difference
% in visual properties and numerical distance.

% To get started: you have to insert the following values below:
% 'design', 'S1', 'S2', 'nrofblocks', 'repmin', 'repmax', 'images', 'balance',
% 'video', 'scalingfactor'.

%possible designs are:
%1 = all numbers are compared
%2 = comparison to a fixed number
%3 = comparison between specific number combinations
%4 = same as design 2 but now in the habituation format

%examples of the different designs are:
% design 1: all posible combinations are made (e.g. 1-1, 1-2, 1-3, etc)
% S1 = [1 2 3 4 5];
% S2 = 0;
%
% design 2: all combinations include S2 (e.g. 8-16, 10-16, 13-16, etc)
% S1 = [8 10 13 20 24 32];
% S2 = 16;
%
% design 3: all combinations are specified (e.g. 1-2, 3-4, 6-6, etc)
% S1 = [1 3 6 9 3 4 8 9 4 4 9 9];
% S2 = [2 4 6 8 1 2 6 7 1 1 6 6];

% design 4: S2 is the standard, S1 the deviant number (e.g. 16 16 16 16 8 16 16 16 etc)
% S1 = [8 10 13 20 24 32];
% S2 = 16;

clear all
commandwindow;

design = 2; % number 1, 2, 3, or 4 see above for explanation
S1 = [12 16 18 20 22 26 29 32 36 48]; %a vector of numbers
S2 = 24; %a single number for designs 1, 2, 3 and a vector of numbers for design 4
nrofblocks = 1; %number of repetitions of all stimuli
repmin = 3; %minimum number of repetitions (this is only taken into account for design 4)
repmax = 7; %maximum number of repetitions (this is only taken into account for design 4)
images = 2; %2 to have the program generate .bmp files of the stimuli otherwise 1
balance = 1; %2 to have S1 - S2 as well as S2 - S1 stimulus pairs otherwise 1
video = 2; %2 to have the stimuli presented on the screen using psychtoolbox otherwise 1
scalingfactor = 1; %a number in the range [0.5,->] can be entered to decrease or increase the image


% Note, .bmp files can only be created when images are presented on the
% screen.
if images == 2;
    video = 2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start

stimulusTime = 0.3; %duration of the stimulus

stimColor = [230 230 230];
windowColor = [50 50 50];


if length(S1)~=length(S2);
    if length(S2) == 1;
        if design == 1;
            S2 = repmat(S1,1,length(S1));%makes all possible combinations
            S1 = Expand(S1,length(S1),1);
        elseif design == 2 || design == 4;
            S2 = repmat(S2,1,length(S1));%increases length of S2 array to length of S1
        end
    end
end

if max(S1) > max(S2);
    maxnr = max(S1);
else
    maxnr = max(S2);
end




[width, height]=Screen('WindowSize', 0);
Wc = width/2;
Hc = height/2;

smallest = round(scalingfactor*5);
largest = smallest*7;

BallSizes = [smallest:largest];

maxball = max(BallSizes);
maxsurf = pi*(maxball/2)^2*maxnr;
maxarea = sqrt((3*maxsurf)/pi)*1.15;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Here the condition matrix is created.

if design == 4 || design ==1;
    balance = 1;
end


trialnr = 0;
for AreaOrder = 1:2;
    for SizeOrder = 1:2;% 1 = small-large / 2 = large-small
        for NumberOrder = 1:balance;% 1 = small-large / 2 = large-small
            for S1Trials = 1:length(S1);
                trialnr = trialnr+1;
                CondMatrix(1,trialnr) = AreaOrder;
                CondMatrix(2,trialnr) = SizeOrder;
                CondMatrix(3,trialnr) = NumberOrder;
                CondMatrix(4,trialnr) = S1(1,S1Trials);
                CondMatrix(5,trialnr) = S2(1,S1Trials);
                if CondMatrix(3,trialnr) == 1;
                    if S1(1,S1Trials)>S2(1,S1Trials); CondMatrix(7,trialnr)=1; CondMatrix(6,trialnr) = (S1(1,S1Trials)-S2(1,S1Trials))/S2(1,S1Trials);
                    else  CondMatrix(7,trialnr) = 2; CondMatrix(6,trialnr) = (S2(1,S1Trials)-(S1(1,S1Trials)))/(S1(1,S1Trials));
                    end
                else
                    if S1(1,S1Trials)>S2(1,S1Trials); CondMatrix(7,trialnr)=2;CondMatrix(6,trialnr) = (S1(1,S1Trials)-S2(1,S1Trials))/S2(1,S1Trials);
                    else  CondMatrix(7,trialnr) = 1;CondMatrix(6,trialnr) = (S2(1,S1Trials)-(S1(1,S1Trials)))/(S1(1,S1Trials));
                    end
                end
            end
        end
    end
end

value = nrofblocks*length(CondMatrix);

%%%%%%%%%%%%%
%skewness: the dots are drawn from a skewed distribution. In this manner
%all dotsizes can appear but some are more likely to appear than others. In
%half the trials, more large dots will be included in the distribution
%whereas in the other half of the trials more small dots are included.

DM = [10 33 10 20];
DM = round(DM*scalingfactor);
skewness = 0.5;

for i = 1:4;
    name = int2str(i);
    DiameterMax = DM(1,i);
    x = repmat(DiameterMax,1,length(BallSizes));
    y = betapdf(skewness,BallSizes,x);
    z = round(y);

    trialnr = 0;
    for t = 1:length(BallSizes);
        value = z(1,t);
        if value > 0;
            for s = 1:value;
                trialnr = trialnr+1;
                if i == 1;
                    matrix1(1,trialnr)= BallSizes(1,t);
                elseif i == 2;
                    matrix2(1,trialnr)= BallSizes(1,t);
                elseif i == 3;
                    matrix3(1,trialnr)= BallSizes(1,t);
                elseif i == 4;
                    matrix4(1,trialnr)= BallSizes(1,t);
                end
            end
        else
            t = t+1;
        end
    end
end


if maxnr > 40;
    matrix1 = repmat(matrix1,1,2);
    matrix2 = repmat(matrix2,1,2);
    matrix3 = repmat(matrix3,1,2);
    matrix4 = repmat(matrix4,1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% possible locations for the dots are defined


maxpos = ceil(maxnr*1.5);
n=round(sqrt(maxpos));
coords = zeros((n+1)*(n+1),2);
a=0;

for i = 0:n
    for j=0:n
        a=a+1;
        coords(a,1)=(i-n/2 + (mod(j,2)-0.5)/4 ) / (n/2);
        coords(a,2)=(j-n/2 + (mod(i,2)-0.5)/4 ) / (n/2);
    end
end

dist = sqrt( coords(:,1).*coords(:,1) + coords(:,2).*coords(:,2) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the screen for presentation

if video == 2;
    [w, scherm] = Screen('OpenWindow',0, windowColor);
    Screen('FillRect', w, windowColor);
    onset = Screen('Flip',w);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the program
trialnr = 0;

clear MatrixCondOut;
for block = 1:nrofblocks;

    TrialListCond = randperm(length(CondMatrix));
    check = 0;
    for j = 1:length(CondMatrix);

        AR = CondMatrix(1,TrialListCond(1,j));
        %number and its order
        W = CondMatrix(6,TrialListCond(1,j));
        if CondMatrix(3,TrialListCond(1,j)) == 1 && design ~=4;
            F = CondMatrix(4,TrialListCond(1,j));
            S = CondMatrix(5,TrialListCond(1,j));
        else
            F = CondMatrix(5,TrialListCond(1,j));
            S = CondMatrix(4,TrialListCond(1,j));
        end



        if  CondMatrix(7,TrialListCond(1,j)) == 1;
            if CondMatrix(2,TrialListCond(1,j)) == 1;
                matrixff = matrix3;
                matrixll = matrix4;
            else
                matrixff = matrix2;
                matrixll = matrix1;
            end
        else
            if CondMatrix(2,TrialListCond(1,j)) ==1;
                matrixff = matrix1;
                matrixll = matrix2;
            else
                matrixff = matrix4;
                matrixll = matrix3;
            end
        end

        if design == 4;
            NrOfRepetitionsF = ceil(repmin-1 + (repmax-repmin-1).*rand(1,1));
        else NrOfRepetitionsF = 1;
        end

        % first number

        % if it is a habituation design, area and diameter and thus
        % density, surface and circumference have to change randomly for
        % the standards. The last standard and the following deviant
        % are controlled in a similar manner as for design 1 to 3.

        for repetition = 1:NrOfRepetitionsF;
            trialnr = trialnr +1;
            if design ~= 4;
                codenr = 1;%deviant number
            else
                if repetition == 1;
                    codenr = 2;%first standard
                elseif repetition == NrOfRepetitionsF;
                    codenr = 3;%last standard
                else codenr = 0;
                end
            end



            if repetition < NrOfRepetitionsF;

                repvar = (ceil(0 + (4-0).*rand(1,1)));

                %variation in visual properties is created for the standards

                if repvar <=3; matrixf = matrixff;

                else matrixf = matrixll;
                end

                areaF = ((0.75 +(1-0.75).*rand(1,1))*maxarea);
                coordsF = coords*areaF;
                distsF = sqrt( coordsF(:,1).*coordsF(:,1) + coordsF(:,2).*coordsF(:,2) ); %dist from coord to centre
                locsF = coordsF(distsF<(areaF+15),:);
                TrialListLocf = randperm(length(locsF));
            else

                if AR==1;
                    areaF = 0.75*maxarea;
                    areaS = maxarea;
                else
                    areaF = maxarea;
                    areaS = 0.75*maxarea;
                end

                coordsF = coords*areaF;
                distsF = sqrt( coordsF(:,1).*coordsF(:,1) + coordsF(:,2).*coordsF(:,2) ); %dist from coord to centre
                locsF = coordsF(distsF<(areaF+15),:);
                TrialListLocf = randperm(length(locsF));
                matrixf = matrixff;
            end



            TrialListBallSize = randperm(length(matrixf));
            for m=1:F;
                Loc = TrialListLocf(1,m);
                BallSize = matrixf(1,TrialListBallSize(1,m));
                a = locsF(Loc,1)+Wc;
                b = locsF(Loc,2)+Hc;
                c = a+BallSize;
                d = b+BallSize;
                sizesF(:,m) = [a b c d];
                surfaceF(1,m) = (pi*((BallSize/2)^2));
                circumfF(1,m) = 2*pi*(BallSize/2);
                diameterF(1,m) = BallSize;
            end
            totSurfaceF = sum(surfaceF(1,1:F));
            avgDiameterF = (sum(diameterF(1,1:F)))/F;
            totCircumfF = sum(circumfF(1,1:F));

            if video == 2;
                Screen('FillOval', w, stimColor,sizesF);
                onset = Screen('Flip',w,onset+stimulusTime);
            end


            %calculating shortest contour around all the dots
            sizesFch(1,1:2*F) = repmat(sizesF(1,:),1,2);
            sizesFch(1,((2*F)+1):4*F) = repmat(sizesF(3,:),1,2);
            sizesFch(2,1:F) = sizesF(2,:);
            sizesFch(2,(F+1):(2*F)) = sizesF(4,:);
            sizesFch(2,(2*F+1):(3*F)) = sizesF(2,:);
            sizesFch(2,(3*F+1):(4*F)) = sizesF(4,:);

            [k,v] = convhull(sizesFch(1,:),sizesFch(2,:));
            CHareaF = v;
            CHdensityF = v/totSurfaceF;


            MatrixCondOut(trialnr,1) = block;
            MatrixCondOut(trialnr,2) = trialnr;
            MatrixCondOut(trialnr,3) = codenr;
            MatrixCondOut(trialnr,4) = F;
            MatrixCondOut(trialnr,5) = W;
            MatrixCondOut(trialnr,6) = CHareaF;
            MatrixCondOut(trialnr,7) = CHdensityF;
            MatrixCondOut(trialnr,8) = totSurfaceF;
            MatrixCondOut(trialnr,9) = avgDiameterF;
            MatrixCondOut(trialnr,10) = totCircumfF;

            clear sizesFch;
            clear sizesF;
            clear matrixf;


            if images == 2;
                xf1 = (width/2)-maxarea-50;
                xf2 = (width/2)+maxarea+50;
                yf1 = (height/2)-maxarea-50;
                yf2 = (height/2)+maxarea+50;

                imageArray=Screen('GetImage', w , [xf1 yf1 xf2 yf2]);
                imageName = strcat((int2str(block)),'_',(int2str(trialnr)),'_1.bmp');
                imwrite(imageArray, imageName);
                picturename{trialnr,1} = imageName;
            end
        end




        %second number
        trialnr = trialnr +1;
        if design ~= 4;
            codenr = 2;%second nr
        else
            codenr = 1;%deviant number
        end

        coordsS = coords*areaS;
        distsS = sqrt( coordsS(:,1).*coordsS(:,1) + coordsS(:,2).*coordsS(:,2) ); %dist from coord to centre
        locsS = coordsS(distsS<(areaS+15),:);
        TrialListLocs = randperm(length(locsS));


        TrialListBallSize = randperm(length(matrixll));
        for p=1:S;
            Loc = TrialListLocs(1,p);
            BallSize = matrixll(1,TrialListBallSize(1,p));
            a = locsS(Loc,1)+Wc;
            b = locsS(Loc,2)+Hc;
            c = a+BallSize;
            d = b+BallSize;
            sizesS(:,p) = [a b c d];
            surfaceS(1,p) = (pi*((BallSize/2)^2));
            circumfS(1,p) = 2*pi*(BallSize/2);
            diameterS(1,p) = BallSize;
        end
        totSurfaceS = sum(surfaceS(1,1:S));
        avgDiameterS = (sum(diameterS(1,1:S)))/S;
        totCircumfS = sum(circumfS(1,1:S));


        if video == 2;
            Screen('FillOval', w, stimColor,sizesS);
            onset = Screen('Flip',w,onset+stimulusTime);
        end

        %calculating shortest contour around array
        sizesSch(1,1:2*S) = repmat(sizesS(1,:),1,2);
        sizesSch(1,((2*S)+1):4*S) = repmat(sizesS(3,:),1,2);
        sizesSch(2,1:S) = sizesS(2,:);
        sizesSch(2,(S+1):(2*S)) = sizesS(4,:);
        sizesSch(2,(2*S+1):(3*S)) = sizesS(2,:);
        sizesSch(2,(3*S+1):(4*S)) = sizesS(4,:);

        [k,v] = convhull(sizesSch(1,:),sizesSch(2,:));
        CHareaS = v;
        CHdensityS = v/totSurfaceS;

        MatrixCondOut(trialnr,1) = block;
        MatrixCondOut(trialnr,2) = trialnr;
        MatrixCondOut(trialnr,3) = codenr;
        MatrixCondOut(trialnr,4) = S;
        MatrixCondOut(trialnr,5) = W;
        MatrixCondOut(trialnr,6) = CHareaS;
        MatrixCondOut(trialnr,7) = CHdensityS;
        MatrixCondOut(trialnr,8) = totSurfaceS;
        MatrixCondOut(trialnr,9) = avgDiameterS;
        MatrixCondOut(trialnr,10) = totCircumfS;

        clear sizesSch;
        clear sizesS;
        clear matrixl;


        if images == 2;
            xs1 = (width/2)-maxarea-50;
            xs2 = (width/2)+maxarea+50;
            ys1 = (height/2)-maxarea-50;
            ys2 = (height/2)+maxarea+50;

            imageArray=Screen('GetImage', w , [xs1 ys1 xs2 ys2]);
            imageName = strcat((int2str(block)),'_',(int2str(trialnr)),'_2.bmp');
            imwrite(imageArray, imageName);
            picturename{trialnr,1} = imageName;
        end


    end
end



if video == 2;

    Screen('CloseAll');
end


%%%%%%%%%%%%%%%%%%%%%%%
%calculating visual properties for absolute difference and weber fraction.


if design ~=4;
    j = 1;
    for i = 1:2:(length(MatrixCondOut)-1);
        ya(j,1) = MatrixCondOut(i,4)-MatrixCondOut(i+1,4);%absolute difference
        yw(j,1) = MatrixCondOut(i,5);%relative difference
        xx(j,1) = MatrixCondOut(i,6)-MatrixCondOut(i+1,6);
        xx(j,2) = MatrixCondOut(i,7)-MatrixCondOut(i+1,7);
        xx(j,3) = MatrixCondOut(i,8)-MatrixCondOut(i+1,8);
        xx(j,4) = MatrixCondOut(i,9)-MatrixCondOut(i+1,9);
        xx(j,5) = MatrixCondOut(i,10)-MatrixCondOut(i+1,10);
        i = i+1;
        j = j+1;
    end
else %comparing average of standards to the number deviant

    j = 1;
    counter = 1;
    value = MatrixCondOut(1,:);
    for i = 3:length(MatrixCondOut);
        if MatrixCondOut(i,4)==S2(1,1);
           value = value + MatrixCondOut(i,:);
           counter = counter+1;
        else
           MatrixCondOut2(j,:) = value/counter;
           counter = 0;
           value = zeros(1,10);
           MatrixCondOut2(j+1,:) = MatrixCondOut(i,:);
           j = j+2;
        end
    end


    j = 1;
    for i = 1:2:(length(MatrixCondOut2)-1);
        ya(j,1) = MatrixCondOut2(i,3)-MatrixCondOut2(i+1,3);%absolute difference
        yw(j,1) = MatrixCondOut2(i,4);%relative difference
        xx(j,1) = MatrixCondOut2(i,5)-MatrixCondOut2(i+1,5);
        xx(j,2) = MatrixCondOut2(i,6)-MatrixCondOut2(i+1,6);
        xx(j,3) = MatrixCondOut2(i,7)-MatrixCondOut2(i+1,7);
        xx(j,4) = MatrixCondOut2(i,8)-MatrixCondOut2(i+1,8);
        xx(j,5) = MatrixCondOut2(i,9)-MatrixCondOut2(i+1,9);
        i = i+1;
        j = j+1;
    end
end

%flipping the same but reversed number pairs (e.g. 3 versus 1 to 1 versus 3)
cg1=0;
ig1=0;
for j = 1:5;
    for i = 1:length(ya);
        if ya(i,1) < 0 ;
            if xx(i,j) < 0;
                regr1(i,j) = abs(xx(i,j));
                regr2(i,j) = abs(ya(i,1));
                cg(i,j) = cg1 +1;
            else
                regr1(i,j) = -(xx(i,j));
                regr2(i,j) = abs(ya(i,1));
                ig(i,j) = ig1+1;
            end
        else
            if xx(i,j) < 0;
                regr1(i,j) = (xx(i,j));
                regr2(i,j) = (ya(i,1));
                ig(i,j) = ig1+1;

            else
                regr1(i,j) = (xx(i,j));
                regr2(i,j) = (ya(i,1));
                cg(i,j) = cg1+1;
            end
        end
    end
end

regr1(:,6) = regr2(:,1);
regr1(:,7) = yw(:,1);
clear j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputfiles:
%
% all.txt comprises the following: trialnumber, codenr(refers to first
% standard, last standard and number deviant stimulus), numerosity presented, weber
% fraction, area subtended, density, total surface of the dots, average diameter,
% total circumference.
%
% regr.txt comprises the difference of the visual properties of the stimuli
% as well as the abolute and relative number distance: area subtended,
% density, total surface of the dots, average diameter, total
% circumference, number distance, weber fraction
%
% Thus the "all" file specifies all values for each stimulus separately
% whereas the "regr" file specifies the differences between the two stimuli
% of each pair. The "regr" file should be used to verify whether no
% relation exists between the difference in visual properties and number
% distance.



dlmwrite('all.txt',MatrixCondOut);
dlmwrite('regr.txt',regr1);
if images == 2;
    namefile = sprintf('%s\n', picturename{:});
    dlmwrite('picturenames.txt',namefile,'');
end

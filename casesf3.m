%%%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%
%%%%%%%%% Masimula Jabulani - Matlab code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Code calculates the load point relibilty and system indices %%%%
%%%%%%%%% Four load points (LP14-LP17)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Ten failure components taken into consideraration %%%%%%%%%%%%%%
%%%%%%%%% Date - 29 August 2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%
% Reset Matlab
%function ECOST = ecostfdr3()
clear all
% Start a stopwatch timer
tic
% User input either SAIDI or ECOST
prompt = 'Do you want to place switches in terms of SAIDI (press 1) or ECOST (press 2)';
x = input(prompt)

if x==1
    % Number of simulations to be perfomed
    N = 100000;
    % Begin a for loop
    for i = 1:N
        % Initialize (Number of failures of each component)
        tf1 = 0; tf2 = 0; tf3 = 0; tf4 = 0; tf5 = 0; tf6 = 0;
        tf7 = 0; tf8 = 0; tf9 = 0; tf10 = 0;
        % Set a new random state
        rand( 'state', sum(100*clock));
        % Repair rates of each load point
        trLp1=0;trLp2=0;trLp3=0;trLp4=0;
        tr2Lp1=0;
        tr3Lp1=0;tr3Lp2=0;tr3Lp3=0;tr3Lp4=0;
        tr4Lp1=0;tr4Lp2=0;tr4Lp3=0;tr4Lp4=0;
        tr5Lp1=0;tr5Lp2=0;tr5Lp3=0;tr5Lp4=0;
        tr6Lp1=0;tr6Lp2=0;tr6Lp3=0;tr6Lp4=0;
        tr7Lp1=0;tr7Lp2=0;tr7Lp3=0;tr7Lp4=0;
        tr8Lp4=0;
        tr9 =0;
        tr10=0;
        % Number of customers at each load point
        d=22;d1=10;d2=1;d3=1;d4=10;a='29';
        % Mission time of 1 year and Failure rate of each component
        t = 1;
        a1 = 0.04875; a2 = 0.039; a3 = 0.052; a4 = 0.04875; a5 = 0.039;
        a6 = 0.052; a7 = 0.04875; a8 = 0.039; a9 = 0.015; a10 = 0.006;
        % Define c = failure rate * number of n years
        c1 = 48.75; c2 = 39; c3 = 52; c4 = 48.75; c5 = 39; c6 = 52; c7 = 48.75;
        c8 = 39; c9 = 15; c10 = 6;
        % Simulation period n years
        n = 1000;
        % start a for loop for n simulation years
        for j = 1:n
            % Generate random number of each component and Convert them to
            b1 = rand; b2 = rand; b3 = rand; b4 = rand; b5 = rand; b6 = rand;
            b7 = rand; b8 = rand; b9 = rand; b10 = rand;
            % Convert random numbers to times-to-failure (TTF)
            T1 = (-1/a1)*reallog(b1); T2 = (-1/a2)*reallog(b2);
            T3 = (-1/a3)*reallog(b3); T4 = (-1/a4)*reallog(b4);
            T5 = (-1/a5)*reallog(b5); T6 = (-1/a6)*reallog(b6);
            T7 = (-1/a7)*reallog(b7); T8 = (-1/a8)*reallog(b8);
            T9 = (-1/a9)*reallog(b9); T10 = (-1/a10)*reallog(b10);
            % Capture the number of failures (tf) and the number of repair rate (tr)
            if T1 < t
                tf1 = tf1 + 1;
                trLp1 = trLp1 + 5/c1;
                trLp2 = trLp2 + 5/c1;
                trLp3 = trLp3 + 5/c1;
                trLp4 = trLp4 + 5/c1;
            end
            if T2 < t
                tf2 = tf2 + 1;
                tr2Lp1 = tr2Lp1 + 5/c2;
            end
            if T3 < t
                tf3 = tf3 +1;
                tr3Lp1 = tr3Lp1 + 1/c3;
                tr3Lp2 = tr3Lp2 + 5/c3;
                tr3Lp3 = tr3Lp3 + 5/c3;
                tr3Lp4 = tr3Lp4 + 5/c3;
            end
            if T4 < t
                tf4 = tf4 + 1;
                %tr4Lp1 = tr4Lp1 + 1/c4;
                tr4Lp2 = tr4Lp2 + 5/c4;
                %tr4Lp3 = tr4Lp3 + 5/c4;
                %tr4Lp4 = tr4Lp4 + 5/c4;
            end
            if T5 < t
                tf5 = tf5 + 1;
                tr5Lp1 = tr5Lp1 + 1/c5;
                tr5Lp2 = tr5Lp2 + 5/c5;
                tr5Lp3 = tr5Lp3 + 5/c5;
                tr5Lp4 = tr5Lp4 + 5/c5;
            end
            if T6 < t
                tf6 = tf6 + 1;
                %tr6Lp1 = tr6Lp1 + 1/c6;
                %tr6Lp2 = tr6Lp2 + 1/c6;
                tr6Lp3 = tr6Lp3 + 5/c6;
                %tr6Lp4 = tr6Lp4 + 5/c6;
            end
            if T7 < t
                tf7 = tf7 + 1;
                tr7Lp1 = tr7Lp1 + 1/c7;
                tr7Lp2 = tr7Lp2 + 5/c7;
                tr7Lp3 = tr7Lp3 + 5/c7;
                tr7Lp4 = tr7Lp4 + 5/c7;
            end
            if T8 < t
                tf8 = tf8 + 1;
                tr8Lp4 = tr8Lp4 + 5/c8;
            end
            if T9 < t
                tf9 = tf9 + 1;
                tr9 = tr9 + 200/c9;
            end
            if T10 < t
                tf10 = tf10 + 1;
                tr10 = tr10 + 4/c10;
            end
        end
        % Finds the failure rate, lamda, for each load point
        lamda1 = (tf1+tf2+tf3+tf4+tf5+tf6+tf7+tf9+tf10)/n;
        lamda2 = (tf1+tf3+tf4+tf5+tf6+tf7+tf10)/n;
        lamda3 = (tf1+tf3+tf4+tf5+tf6+tf7+tf10)/n;
        lamda4 = (tf1+tf3+tf4+tf5+tf6+tf7+tf8+tf9+tf10)/n;
        lamda1=lamda1-0.093; lamda4=lamda4-0.093;
        % Finds the Unavailability, U, for each load point
        ULp1=(tf1/n).*trLp1+(tf2/n).*tr2Lp1+(tf3/n).*tr3Lp1+(tf4/n).*tr4Lp1+(tf5/n).*tr5Lp1+(tf6/n).*tr6Lp1+(tf7/n).*tr7Lp1+(tf9/n).*tr9+(tf10/n).*tr10;
        ULp2=(tf1/n).*trLp2+(tf3/n).*tr3Lp2+(tf4/n).*tr4Lp2+(tf5/n).*tr5Lp2+(tf6/n).*tr6Lp2+(tf7/n).*tr7Lp2+(tf10/n).*tr10;
        ULp3=(tf1/n).*trLp3+(tf3/n).*tr3Lp3+(tf4/n).*tr4Lp3+(tf5/n).*tr5Lp3+(tf6/n).*tr6Lp3+(tf7/n).*tr7Lp3+(tf10/n).*tr10;
        ULp4=(tf1/n).*trLp4+(tf3/n).*tr3Lp4+(tf4/n).*tr4Lp4+(tf5/n).*tr5Lp4+(tf6/n).*tr6Lp4+(tf7/n).*tr7Lp4+(tf8/n).*tr8Lp4+(tf9/n).*tr9+(tf10/n).*tr10;
        % Compute Outage duration, r, for each load point
        r1 = ULp1/lamda1;
        r2 = ULp2/lamda2;
        r3 = ULp3/lamda3;
        r4 = ULp4/lamda4;
        % Determine System indices using equation discussed in section 2 of report
        % Find SAIFI
        SAIFI = ((lamda1*d1)+(lamda2*d2)+(lamda3*d3)+(lamda4*d4))/d;
        % Find SAIDI
        SAIDI = ((ULp1*d1)+(ULp2*d2)+(ULp3*d3)+(ULp4*d4))/d;
        % Find CAIDI
        num = ((ULp1*d1)+(ULp2*d2)+(ULp3*d3)+(ULp4*d4));
        den = ((lamda1*d1)+(lamda2*d2)+(lamda3*d3)+(lamda4*d4));
        CAIDI = num/den;
        % Compute ASAI
        ASAI = ((d*8760) - num)/(d*8760);
        % Find ASUI
        ASUI = 1 - ASAI;
        % Calculates the ENS
        ENS = ((ULp1*0.4697)+(ULp2*1.6391)+(ULp3*0.9025)+(ULp4*0.4697));
        % Determine the  AENS
        AENS = ENS/d;
        % Determining Reliability cost/worth indices LP14-LP17
        COST1 = ( ((83.01-31.32)/(8-4))*r1 + (83.01 -  ((83.01-31.32)/(8-4))*8))*0.4697;
        COST2 = ( ((83.01-31.32)/(8-4))*r2 + (83.01 -  ((83.01-31.32)/(8-4))*8))*1.6391;
        COST3 = ( ((83.01-31.32)/(8-4))*r3 + (83.01 -  ((83.01-31.32)/(8-4))*8))*0.9025;
        COST4 = ( ((83.01-31.32)/(8-4))*r4 + (83.01 -  ((83.01-31.32)/(8-4))*8))*0.4697;
        %Total ECOST
        ECOST = (COST1*lamda1)+(COST2*lamda2)+(COST3*lamda3)+(COST4*lamda4);
        % IAER calculation
        IAER=(COST1*lamda1)/(ULp1*0.4697)+(COST2*lamda2)/(ULp2*1.6391)+(COST3*lamda3)/(ULp3*0.9025)+(COST4*lamda4)/(ULp4*0.4697);
        
        % Display Matrix of results
        results(i,:) = [lamda1 lamda2 lamda3 lamda4 r1 r2 r3 r4 ULp1 ULp2 ULp3 ULp4];
        y(i,:) = [SAIDI  ENS  ECOST];
        % Matrix of SAIFI and SAIDI
        s1(i,:) = [SAIFI];
        s2(i,:) = [SAIDI];
       
    end
  disp(['The Swicth is in position :', a]);
else
  
    % Number of simulations to be perfomed
    N = 100000;
    % Begin a for loop
    for i = 1:N
        % Initialize (Number of failures of each component)
        tf1 = 0; tf2 = 0; tf3 = 0; tf4 = 0; tf5 = 0; tf6 = 0;
        tf7 = 0; tf8 = 0; tf9 = 0; tf10 = 0;
        % Set a new random state
        rand( 'state', sum(100*clock));
        % Repair rates of each load point
        trLp1=0;trLp2=0;trLp3=0;trLp4=0;
        tr2Lp1=0;
        tr3Lp1=0;tr3Lp2=0;tr3Lp3=0;tr3Lp4=0;
        tr4Lp1=0;tr4Lp2=0;tr4Lp3=0;tr4Lp4=0;
        tr5Lp1=0;tr5Lp2=0;tr5Lp3=0;tr5Lp4=0;
        tr6Lp1=0;tr6Lp2=0;tr6Lp3=0;tr6Lp4=0;
        tr7Lp1=0;tr7Lp2=0;tr7Lp3=0;tr7Lp4=0;
        tr8Lp4=0;
        tr9 =0;
        tr10=0;
        % Number of customers at each load point
        d=22;d1=10;d2=1;d3=1;d4=10; e='31';
        % Mission time of 1 year and Failure rate of each component
        t = 1;
        a1 = 0.04875; a2 = 0.039; a3 = 0.052; a4 = 0.04875; a5 = 0.039;
        a6 = 0.052; a7 = 0.04875; a8 = 0.039; a9 = 0.015; a10 = 0.006;
        % Define c = failure rate * number of n years
        c1 = 48.75; c2 = 39; c3 = 52; c4 = 48.75; c5 = 39; c6 = 52; c7 = 48.75;
        c8 = 39; c9 = 15; c10 = 6; 
        % Simulation period n years
        n = 1000;
        % start a for loop for n simulation years
        for j = 1:n
            % Generate random number of each component and Convert them to
            b1 = rand; b2 = rand; b3 = rand; b4 = rand; b5 = rand; b6 = rand;
            b7 = rand; b8 = rand; b9 = rand; b10 = rand;
            % Convert random numbers to times-to-failure (TTF)
            T1 = (-1/a1)*reallog(b1); T2 = (-1/a2)*reallog(b2);
            T3 = (-1/a3)*reallog(b3); T4 = (-1/a4)*reallog(b4);
            T5 = (-1/a5)*reallog(b5); T6 = (-1/a6)*reallog(b6);
            T7 = (-1/a7)*reallog(b7); T8 = (-1/a8)*reallog(b8);
            T9 = (-1/a9)*reallog(b9); T10 = (-1/a10)*reallog(b10);
            % Capture the number of failures (tf) and the number of repair rate (tr)
            if T1 < t
                tf1 = tf1 + 1;
                trLp1 = trLp1 + 5/c1;
                trLp2 = trLp2 + 5/c1;
                trLp3 = trLp3 + 5/c1;
                trLp4 = trLp4 + 5/c1;
            end
            if T2 < t
                tf2 = tf2 + 1;
                tr2Lp1 = tr2Lp1 + 5/c2;
            end
            if T3 < t
                tf3 = tf3 +1;
                tr3Lp1 = tr3Lp1 + 5/c3;
                tr3Lp2 = tr3Lp2 + 5/c3;
                tr3Lp3 = tr3Lp3 + 5/c3;
                tr3Lp4 = tr3Lp4 + 5/c3;
            end
            if T4 < t
                tf4 = tf4 + 1;
                %tr4Lp1 = tr4Lp1 + 1/c4;
                tr4Lp2 = tr4Lp2 + 5/c4;
                %tr4Lp3 = tr4Lp3 + 5/c4;
                %tr4Lp4 = tr4Lp4 + 5/c4;
            end
            if T5 < t
                tf5 = tf5 + 1;
                tr5Lp1 = tr5Lp1 + 1/c5;
                tr5Lp2 = tr5Lp2 + 1/c5;
                tr5Lp3 = tr5Lp3 + 5/c5;
                tr5Lp4 = tr5Lp4 + 5/c5;
            end
            if T6 < t
                tf6 = tf6 + 1;
                %tr6Lp1 = tr6Lp1 + 1/c6;
                %tr6Lp2 = tr6Lp2 + 1/c6;
                tr6Lp3 = tr6Lp3 + 5/c6;
                %tr6Lp4 = tr6Lp4 + 5/c6;
            end
            if T7 < t
                tf7 = tf7 + 1;
                tr7Lp1 = tr7Lp1 + 1/c7;
                tr7Lp2 = tr7Lp2 + 1/c7;
                tr7Lp3 = tr7Lp3 + 5/c7;
                tr7Lp4 = tr7Lp4 + 5/c7;
            end
            if T8 < t
                tf8 = tf8 + 1;
                tr8Lp4 = tr8Lp4 + 5/c8;
            end
            if T9 < t
                tf9 = tf9 + 1;
                tr9 = tr9 + 200/c9;
            end
            if T10 < t
                tf10 = tf10 + 1;
                tr10 = tr10 + 4/c10;
            end
        end
        b=2;
        % Finds the failure rate, lamda, for each load point
        lamda1 = (tf1+tf2+tf3+tf4+tf5+tf6+tf7+tf9+tf10)/n;
        lamda2 = (tf1+tf3+tf4+tf5+tf6+tf7+tf10)/n;
        lamda3 = (tf1+tf3+tf4+tf5+tf6+tf7+tf10)/n;
        lamda4 = (tf1+tf3+tf4+tf5+tf6+tf7+tf8+tf9+tf10)/n;
        lamda1=lamda1-0.093; lamda4=lamda4-0.093;
        % Finds the Unavailability, U, for each load point
        ULp1=(tf1/n).*trLp1+(tf2/n).*tr2Lp1+(tf3/n).*tr3Lp1+(tf4/n).*tr4Lp1+(tf5/n).*tr5Lp1+(tf6/n).*tr6Lp1+(tf7/n).*tr7Lp1+(tf9/n).*tr9+(tf10/n).*tr10;
        ULp2=(tf1/n).*trLp2+(tf3/n).*tr3Lp2+(tf4/n).*tr4Lp2+(tf5/n).*tr5Lp2+(tf6/n).*tr6Lp2+(tf7/n).*tr7Lp2+(tf10/n).*tr10;
        ULp3=(tf1/n).*trLp3+(tf3/n).*tr3Lp3+(tf4/n).*tr4Lp3+(tf5/n).*tr5Lp3+(tf6/n).*tr6Lp3+(tf7/n).*tr7Lp3+(tf10/n).*tr10;
        ULp4=(tf1/n).*trLp4+(tf3/n).*tr3Lp4+(tf4/n).*tr4Lp4+(tf5/n).*tr5Lp4+(tf6/n).*tr6Lp4+(tf7/n).*tr7Lp4+(tf8/n).*tr8Lp4+(tf9/n).*tr9+(tf10/n).*tr10;
        % Compute Outage duration, r, for each load point
        r1 = ULp1/lamda1;
        r2 = ULp2/lamda2;
        r3 = ULp3/lamda3;
        r4 = ULp4/lamda4;
        % Determine System indices using equation discussed in section 2 of report
        % Find SAIFI
        SAIFI = ((lamda1*d1)+(lamda2*d2)+(lamda3*d3)+(lamda4*d4))/d;
        % Find SAIDI
        SAIDI = ((ULp1*d1)+(ULp2*d2)+(ULp3*d3)+(ULp4*d4))/d;
        % Find CAIDI
        num = ((ULp1*d1)+(ULp2*d2)+(ULp3*d3)+(ULp4*d4));
        den = ((lamda1*d1)+(lamda2*d2)+(lamda3*d3)+(lamda4*d4));
        CAIDI = num/den;
        % Compute ASAI
        ASAI = ((d*8760) - num)/(d*8760);
        % Find ASUI
        ASUI = 1 - ASAI;
        % Calculates the ENS
        ENS = ((ULp1*0.4697)+(ULp2*1.6391)+(ULp3*0.9025)+(ULp4*0.4697));
        % Determine the  AENS
        AENS = ENS/d;
        % Determining Reliability cost/worth indices LP14-LP17
        COST1 = ( ((83.01-31.32)/(8-4))*r1 + (83.01 -  ((83.01-31.32)/(8-4))*8))*0.4697;
        COST2 = ( ((83.01-31.32)/(8-4))*r2 + (83.01 -  ((83.01-31.32)/(8-4))*8))*1.6391;
        COST3 = ( ((83.01-31.32)/(8-4))*r3 + (83.01 -  ((83.01-31.32)/(8-4))*8))*0.9025;
        COST4 = ( ((83.01-31.32)/(8-4))*r4 + (83.01 -  ((83.01-31.32)/(8-4))*8))*0.4697;
        %Total ECOST
        ECOST = (COST1*lamda1)+(COST2*lamda2)+(COST3*lamda3)+(COST4*lamda4);
        % IAER calculation
        IAER=(COST1*lamda1)/(ULp1*0.4697)+(COST2*lamda2)/(ULp2*1.6391)+(COST3*lamda3)/(ULp3*0.9025)+(COST4*lamda4)/(ULp4*0.4697);
        
        % Display Matrix of results
        results(i,:) = [lamda1 lamda2 lamda3 lamda4 r1 r2 r3 r4 ULp1 ULp2 ULp3 ULp4];
        y(i,:) = [ SAIDI  ENS  ECOST ];
        % Matrix of SAIFI and SAIDI
        s1(i,:) = [SAIFI];
        s2(i,:) = [SAIDI];
  %disp(['The Swicth is in position :', e]);
    end
    disp(['The Swicth is in position :', e]);
end

% Outputs
results;
disp('Respective Load points indices are :')
x = (sum(results))/N
disp('The SAIDI, ENS & ECOST are: ')
yy = (sum(y))/N
% Mean and Standard deviation of SAIFI and SAIDI
% msf = mean(s1)
% msd = mean(s2)
% stdsf = std(s1)
% stdsd = std(s2)
% Read stop watch
toc
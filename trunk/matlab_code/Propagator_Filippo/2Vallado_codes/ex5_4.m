%     -----------------------------------------------------------------
%
%                              Ex5_4.m
%
%  this file demonstrates example 5-4.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            22 jan 11  david vallado
%                         original
%  changes :
%            22 jan 11  david vallado
%                         original baseline
%
%     *****************************************************************

constmath;

% --------  moon         - moon rise set
jd = jday( 1998, 8, 21, 0, 0, 0.00 );
latgd = 40.0/rad;
lon = 0.00 / rad;

[utmoonrise,utmoonset,moonphaseang,error] = moonrise( jd,latgd,lon,'y' )
fprintf(1,'moon moonrise %14.4f    moonset %14.4f hrs \n',utmoonrise,utmoonset );
fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );


jd = jday( 1990, 3, 5, 0, 0, 0.00 );
latgd =  40.94 /rad;
lon   = -73.97 / rad;

%        [utmoonrise,utmoonset,moonphaseang,error] = moonrise1( jd,latgd,lon,'y' )
%        fprintf(1,'moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n',utmoonrise, (utmoonrise-floor(utmoonrise))*60, utmoonset, (utmoonset-floor(utmoonset))*60 );
%        fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );

[utmoonrise,utmoonset,moonphaseang,error] = moonrise2( jd,latgd,lon,'y' )
fprintf(1,'moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n',utmoonrise, (utmoonrise-floor(utmoonrise))*60, utmoonset, (utmoonset-floor(utmoonset))*60 );
fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );



jd = jday( 2006, 6, 28, 0, 0, 0.00 );
latgd = 40.0/rad;
lon = 0.00 / rad;

[utmoonrise,utmoonset,moonphaseang,error] = moonrise2( jd,latgd,lon,'y' )
fprintf(1,'moon moonrise %14.4f    moonset %14.4f hrs \n',utmoonrise,utmoonset );
fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );

pause;

fprintf(1,'     40    42    44    46    48    50    52    54    56    58    60    62    64    66  \n' );
for i = 8:30  %8:30
    jd = jday( 2006, 6, i, 0, 0, 0.00 );
    for j = 0:13  %0:13
        latgd = (40.0 + j * 2.0)/rad;
        lon = 0.00 / rad;
        
        [utmoonrise,utmoonset,moonphaseang,error] = moonrise2( jd,latgd,lon,'n' );
        %                if strcmp(error,'ok') == 0 % 1 if true, 0 if false
        %                    fprintf(1,'error');
        %                end;
        
        [hr,min,sec] = rad2hms( utmoonrise*15.0 * pi/180.0 ); % be sure to take hrs and convert to rad befor separating
        [hr1,min1,sec1] = rad2hms( utmoonset*15.0 * pi/180.0 ); % be sure to take hrs and convert to rad befor separating
        %  print out header date for each section of results
        if j == 0
            fprintf(1,'%2i ',i );
        end;
        
%         if utmoonrise > 9998.0 utmoonset > 9998.0 
%              fprintf(1,' none ');
%         else
%         if utmoonrise > 9998.0 && utmoonset < 24.0 
%              fprintf(1,' nors ');
%         else 
% %        if utmoonset - utmoonrise > 14.0 
% %             fprintf(1,' none ');
% %        else
%             fprintf(1,'%2i:%2i ',hr, min);            
% %        end;
%         end;
%         end;
            
%            if hr >= 24
                jd = jday( 2006, 6, i, 0, 0, 0.00 );
                [el1] = moonel( jd,latgd,lon );
                jd = jday( 2006, 6, i+1, 0, 0, 0.00 );
                [el2] = moonel( jd,latgd,lon );

%                if utmoonrise > 9998.0
%                    fprintf(1,' nors ');
%                else
%                    if utmoonset > 9998.0
%                        fprintf(1,' nost ');
%                    else
%                        if (el1 < 0.0 && el2 < 0.0 && utmoonrise > 9998.0 && utmoonset > 9998.0)
%                            fprintf(1,' nors ');
%                        else
%                        if (el1 > 0.0 && el2 > 0.0 && utmoonrise > 9998.0 && utmoonset > 9998.0)
%                            fprintf(1,' nost ');
%                        else
%                            fprintf(1,'%2i:%2i ',hr, min);
% %                            fprintf(1,'%3.0f|%3.0f ',el1*rad, el2*rad);
%                        end    
%                    end;
%                end;
%             end;
                  % print all out
                fprintf(1,'| %2i:%2i %3.0f:%3.0f %3.0f:%3.0f ',hr, min, el1*rad, el2*rad, utmoonrise, utmoonset);
         end; % if j = 0
        fprintf(1,'\n');
    end;
        
    
    
    
%
% ------------------------------------------------------------------------------
%
%                           function coe2rv
%
%  this function finds the position and velocity vectors in geocentric
%    equatorial (ijk) system given the classical orbit elements.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%    vallado     - add constant file use                         29 jun 2003
%
%  inputs          description                    range / units
%    p           - semilatus rectum               km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%
%  outputs       :
%    r           - ijk position vector            km
%    v           - ijk velocity vector            km / s
%
%  locals        :
%    temp        - temporary real*8 value
%    rpqw        - pqw position vector            km
%    vpqw        - pqw velocity vector            km / s
%    sinnu       - sine of nu
%    cosnu       - cosine of nu
%    tempvec     - pqw velocity vector
%
%  coupling      :
%    mag         - magnitude of a vector
%    rot3        - rotation about the 3rd axis
%    rot1        - rotation about the 1st axis
%
%  references    :
%    vallado       2007, 126, alg 10, ex 2-5
%
% [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper );
% ------------------------------------------------------------------------------

function [r,v] = coe2rv ( p0 );

    mu=3.986004e5;
	ecc=p0(2);
	m=p0(6)*pi/180;

	[e0,nu] = newtonm ( ecc,m );	
                
        p=p0(1)*0.001*(1-ecc^2);
        incl=p0(3)*pi/180;
        omega=p0(4)*pi/180;
        argp=p0(5)*pi/180;

        % ----------  form pqw position and velocity vectors ----------
        cosnu= cos(nu);
        sinnu= sin(nu);
        temp = p / (1.0  + ecc*cosnu);
        rpqw(1)= temp*cosnu;
        rpqw(2)= temp*sinnu;
        rpqw(3)=     0.0;
        if ( abs(p) < 0.0001)
            p= 0.0001;
        end
        vpqw(1)=    -sinnu*sqrt(mu)  / sqrt(p);
        vpqw(2)=  (ecc + cosnu)*sqrt(mu) / sqrt(p);
        vpqw(3)=      0.0;

        % ----------------  perform transformation to ijk  ------------
        [tempvec] = rot3( rpqw   , -argp );
        [tempvec] = rot1( tempvec, -incl );
        [r] = rot3( tempvec, -omega );

        [tempvec] =rot3( vpqw   , -argp );
        [tempvec] =rot1( tempvec, -incl );
        [v] = rot3( tempvec, -omega );

        r=1000*r'
        v=1000*v'


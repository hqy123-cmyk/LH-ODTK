function [L, LP, F, D, OM] = FUNDARG(T)
%+
%  - - - - - - - - - - -
%   F U N D A R G 
%  - - - - - - - - - - -
%
%  This routine is part of the International Earth Rotation and
%  Reference Systems Service (IERS) Conventions software collection.
%
%  This subroutine computes the lunisolar fundamental arguments.
%  The model used is from Simon et al. (1994) as recommended by the IERS
%  Conventions (2010).  Refer to IERS Conventions (2010) Chapter 5 
%  Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
%
%  In general, Class 1, 2, and 3 models represent physical effects that
%  act on geodetic parameters while canonical models provide lower-level
%  representations or basic computations that are used by Class 1, 2, or
%  3 models.
% 
%  Status: Canonical model
%
%     Class 1 models are those recommended to be used a priori in the
%     reduction of raw space geodetic data in order to determine
%     geodetic parameter estimates.
%     Class 2 models are those that eliminate an observational
%     singularity and are purely conventional in nature.
%     Class 3 models are those that are not required as either Class
%     1 or 2.
%     Canonical models are accepted as is and cannot be classified as a
%     Class 1, 2, or 3 model.
%
%  Given:
%     T           d      TT, Julian centuries since J2000 (Note 1)
%
%  Returned:
%     L           d      Mean anomaly of the Moon (Note 2)
%     LP          d      Mean anomaly of the Sun (Note 2)
%     F           d      L - OM (Notes 2 and 3)
%     D           d      Mean elongation of the Moon from the Sun
%                                                         (Note 2)
%     OM          d      Mean longitude of the ascending node of
%                                                the Moon (Note 2)
%
%  Notes:
%
%  1) Though T is strictly TDB, it is usually more convenient to use
%     TT, which makes no significant difference.  Julian centuries since
%     J2000 is (JD - 2451545.0)/36525.
%
%  2) The expression used is as adopted in IERS Conventions (2010) and
%     is from Simon et al. (1994).  Arguments are in radians.
%
%  3) L in this instance is the Mean Longitude of the Moon. OM is the 
%     Mean longitude of the ascending node of the Moon.
%
%  Test case:
%     given input: T = 0.07995893223819302 Julian centuries since J2000
%                  (MJD = 54465)
%     expected output:  L = 2.291187512612069099 radians
%                       LP = 6.212931111003726414 radians
%                       F = 3.658025792050572989 radians
%                       D = 4.554139562402433228 radians
%                       OM = -0.5167379217231804489 radians
%
%  References:
%
%     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
%     Francou, G., Laskar, J., 1994, Astron.Astrophys. 282, 663-683
%
%     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
%     IERS Technical Note No. 36, BKG (2010)
%
%  Revisions:
%  2008 January 18 B.E.Stetzler  Initial changes to header
%               and used 2PI instead of PI as parameter
%  2008 January 25 B.E. Stetzler Additional changes to header
%               and defined fundamental arguments
%  2008 January 28 B.E. Stetzler Additional changes to header
%  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
%  2008 April   03 B.E. Stetzler Provided example test case.
%  2009 February 11 B.E. Stetzler Corrected term in OM from 6962890.2665
%                                 to 6962890.5431 and updated test case
%  2009 May     07 B.E. Stetzler Code formatting changes based on 
%                                client recommendations
%  2009 May     07 B.E. Stetzler Updated test case due to above changes
%  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
%-----------------------------------------------------------------------

   %  Arcseconds to radians
   DAS2R = 4.848136811095359935899141D-6;

   %  Arcseconds in a full circle
      TURNAS = 1296000D0;

   %  2Pi
   D2PI = 6.283185307179586476925287D0;

   % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   %  Compute the fundamental argument L.
      L = mod( 485868.249036D0 + ...
                        T*( 1717915923.2178D0 + ...
                        T*(  31.8792D0 + ...
                        T*(   0.051635D0 + ...
                        T*( - 0.00024470D0 )))), TURNAS ) * DAS2R;

   %  Compute the fundamental argument LP.
   LP = mod( 1287104.79305D0 + ...
                  T*( 129596581.0481D0 + ...
                  T*(       - 0.5532D0 + ...
                  T*(         0.000136D0 + ...
                  T*(       - 0.00001149D0 )))), TURNAS ) * DAS2R;

   %  Compute the fundamental argument F.
   F  = mod( 335779.526232D0 + ...
               T*( 1739527262.8478D0 + ...
               T*(       - 12.7512D0 + ...
               T*(       -  0.001037D0 + ...
               T*(          0.00000417D0 )))), TURNAS ) * DAS2R;

   %  Compute the fundamental argument D.
   D = mod( 1072260.70369D0 + ...
             T*( 1602961601.2090D0 + ...
             T*(        - 6.3706D0 + ...
             T*(          0.006593D0 + ...
             T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R;

   %  Compute the fundamental argument OM.
   OM = mod( 450160.398036D0 + ...
             T*( - 6962890.5431D0 + ...
             T*(         7.4722D0 + ...
             T*(         0.007702D0 + ...
             T*(       - 0.00005939D0 )))), TURNAS ) * DAS2R;

   %  Finished.

end
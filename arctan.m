function lamda = arctan(y, x)
% ! ----------------------------------------------------------------------
% ! Purpose:
% !  Angle computation based on the x,y Cartesian coordinates.
% !  Test for identifying the angle quadrant.
% ! ----------------------------------------------------------------------
% ! Remarks:
% !  Angles are considered counter-clockwise
% !  The result is computed in radians. 
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - x,y:			2D Cartesian coordinates
% ! Output arguments:
% ! - angle:			Orientation angle (counter-clockwise) starting from x axis (radians)
% ! ----------------------------------------------------------------------

  if (x ~= 0.d0)  % ~=等价于!=
      a = atan( abs( y/x ) );
      if (x > 0.0d0)
          if (y > 0.0d0)
              angle = a;
          elseif (y < 0.0d0)
              angle = 2.0d0 * pi - a;
          else
              angle = 0.0d0;
          end
      elseif (x < 0.0d0)
          if (y > 0.0d0)
              angle = pi - a;
          elseif (y < 0.0d0)
              angle = pi + a;
          else
              angle = pi;
          end
      end
  else
      if (y > 0.0d0)
          angle = pi / 2.0d0;
      elseif (y < 0.0d0)
          angle = 3.0d0 * pi / 2.0d0;
      else
          angle = 0.d0;
      end
  end

end
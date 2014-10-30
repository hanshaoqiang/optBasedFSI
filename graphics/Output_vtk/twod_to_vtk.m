function twod_to_vtk ( xy, e_conn, uvp, output_filename, title )

%*****************************************************************************80
%
%% TWOD_TO_VTK writes out a TWOD dataset to a legacy VTK file.
%
%  Discussion:
%
%    The VTK file can be read and displayed by the Paraview program.
%
%    Thanks to Mike Sussman for suggesting that real data should be
%    written with the "double" attribute rather than the "float",
%    JVB, 20 December 2010.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 December 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real XY(NODE_NUM,2), the node coordinates.
%
%    Input, integer E_CONN(ELEMENT_NUM,ELEMENT_ORDER), the nodes that
%    form each element.  Node indices are 1-based.
%
%    Input, real UVP(NODE_NUM,3), U and V velocity components and 
%    pressure at each node.
%
%    Input, string OUTPUT_FILENAME, the name of the output file.
%    By convention, this file should have the extension ".vtk".
%
%    Input, string TITLE, a title for the data.
%

%
%  Determine the sizes of things.
%
  [ node_num, dim_num ] = size ( xy );
  [ element_num,  element_order ] = size ( e_conn );
%
%  Open the output file.
%
  if ( isempty ( output_filename ) )
    output_filename = 'ns2d_fem.vtk';
  end

  output_unit = fopen ( output_filename, 'w' );
%
%  Transpose or otherwise modify the data.
%
  xyz = zeros ( 3, node_num );
  xyz(1:2,1:node_num) = xy(1:node_num,1:2)';

  if ( element_order == 6 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TWO_TO_VTK:\n' );
    fprintf ( 1, '  The input data uses quadratic elements.\n' );
    fprintf ( 1, '  The output data will use linear elements.\n' );
  end

  element_order2 = 3;

  element_node = zeros ( element_order2, element_num );
  element_node(1:element_order2,1:element_num) = ...
    e_conn(1:element_num,1:element_order2)' - 1;

  p = zeros ( 1, node_num );
  p(1,1:node_num) = uvp(1:node_num,3)';

  uvw = zeros ( 3, node_num );
  uvw(1:2,1:node_num) = uvp(1:node_num,1:2)';
%
%  Write the data.
%
  vtk_puv_write ( output_unit, title, node_num, element_num, ...
    element_order2, xyz, element_node, p, uvw );

  fclose ( output_unit );

%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  The data was written to "%s"\n', output_filename );

  return
end

function vtk_puv_write ( output_unit, title, node_num, element_num, ...
  element_order, xyz, element_node, p, uvw )

%*****************************************************************************80
%
%% VTK_PUV_WRITE writes pressure and velocity data to a VTK file.
%
%  Discussion:
%
%    The data is assumed to have been computed by a finite element program
%    for a 2D geometry which has been meshed using triangular elements 
%    of 3 or 6 nodes.
%
%    The solution data includes the pressure and velocity vector at each node.
%
%    Note that the VTK format used here is known as the "legacy" or "old style"
%    format.  It has been superseded by a family of XML based formats.  The
%    appropriate replacement for the VTK format used here is known as "VTU",
%    which is the Visual Toolkit format for unstructured grid data.
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
%    Input, integer OUTPUT_UNIT, the output unit.
%
%    Input, string TITLE, a title for the data.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, real XYZ(3,NODE_NUM), the node coordinates.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
%    nodes that make up each element.  Node indices are zero-based.
%
%    Input, real P(1,NODE_NUM), the pressure at each node.
%
%    Input, real UVW(3,NODE_NUM), the velocity at each node.
%
  fprintf ( output_unit, '# vtk DataFile Version 2.0\n' );
  fprintf ( output_unit, '%s\n', title );
  fprintf ( output_unit, 'ASCII\n' );
  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'DATASET UNSTRUCTURED_GRID\n' );
  fprintf ( output_unit, 'POINTS %d double\n', node_num );

  for node = 1 : node_num
    fprintf ( output_unit, '  %f  %f  %f\n', xyz(1:3,node) );
  end
%
%  Note that CELL_SIZE uses ELEMENT_ORDER+1 because the order of each element
%  is included as a data item.
%
  cell_size = element_num * ( element_order + 1 );

  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
  for element = 1 : element_num
    fprintf ( output_unit, '  %d', element_order );
    for order = 1 : element_order
      fprintf ( output_unit, '  %d', element_node(order,element) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  VTK has a cell type 22 for quadratic triangles.  However, we
%  are going to strip the data down to linear triangles for now,
%  which is cell type 5.
%
  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'CELL_TYPES %d\n', element_num );

  if ( element_order == 3 )
    for element = 1 : element_num
      fprintf ( output_unit, '5\n' );
    end
  elseif ( element_order == 6 )
    for element = 1 : element_num
      fprintf ( output_unit, '22\n' );
    end
  end

  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'POINT_DATA %d\n', node_num );
  fprintf ( output_unit, 'SCALARS pressure double\n' );
  fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '  %f\n', p(1,node) );
  end
  fprintf ( output_unit, 'VECTORS velocity double\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '  %f  %f  %f\n', uvw(1:3,node) );
  end

  return
end

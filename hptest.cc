/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2011 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#define USE_FACE_COUPLING
// #define USE_JIMMY_MODEL //BC2
#define USE_LINEAR_ELASTICITY
#define USE_TRIANGLE
#define IS_AXISYMMETRIC
// #define USE_TERM13_FOR_AXISYMMETRY
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/table.h>

#include <deal.II/differentiation/sd/symengine_math.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

// remeshing stuff
#include <deal.II/gmsh/utilities.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>


using namespace dealii;
void test(const std::string &file_name)
{
  enum
  {
    fluid_domain_id,
    hydrogel_domain_id
  };

  const unsigned int dim = 2;
  const unsigned int degree = 2;
  FE_SimplexP<dim> volume_fe(degree);
  Triangulation<dim> triangulation(Triangulation<dim>::maximum_smoothing);
  std::ifstream grid_file(file_name);

  GridIn<2, 2> gridin;
  gridin.attach_triangulation(triangulation);
  gridin.read_msh(grid_file);

  GridOut grid_out;
  std::ofstream out(file_name + "mesh.out.vtk");
  grid_out.write_vtk(triangulation, out);

  hp::FECollection<dim> volume_fe_collection;
  volume_fe_collection.push_back(FE_Nothing<dim>(ReferenceCells::Triangle));
  volume_fe_collection.push_back(volume_fe);

  DoFHandler<dim>       volume_dof_handler(triangulation);
  volume_dof_handler.distribute_dofs(volume_fe_collection);


  hp::QCollection<dim> quadrature_collection;
  const QGaussSimplex<dim> quadrature_formula(2 + degree);
  quadrature_collection.push_back(quadrature_formula);
  quadrature_collection.push_back(quadrature_formula);


  std::unique_ptr<MappingFE<dim>> mapping_pointer;
  mapping_pointer =
    std::make_unique<MappingFE<dim>>(FE_SimplexP<dim>(1));
  hp::MappingCollection<dim> mapping_collection;
  mapping_collection.push_back(*mapping_pointer);
  mapping_collection.push_back(*mapping_pointer);

  // set material id:
  const double inner_radius = 1;
  Point<dim>         quarter_circle_center; // initialize as 0
  for (const auto &cell : volume_dof_handler.active_cell_iterators())
    if (quarter_circle_center.distance(cell->center()) < inner_radius)
        cell->set_material_id(hydrogel_domain_id);
     else
      cell->set_material_id(fluid_domain_id);

  const double hmin = GridTools::minimal_cell_diameter(triangulation,
                                          *mapping_pointer);

  const double hmax =
    GridTools::maximal_cell_diameter(triangulation,
                                     *mapping_pointer);

  std::cout << " hmin = " << hmin << " hmax = " << hmax << std::endl
            << std::endl;

  // set active fe
  const unsigned int stokes_active_fe_index   = 0;
  const unsigned int hydrogel_active_fe_index = 1;
  for (const auto &cell : volume_dof_handler.active_cell_iterators())
    {
      if (cell->material_id() == hydrogel_domain_id)
        {
          cell->set_active_fe_index(hydrogel_active_fe_index);
        }
      else if (cell->material_id() == fluid_domain_id)
        {
          cell->set_active_fe_index(stokes_active_fe_index);
        }
      else
        Assert(
              false,
              ExcNotImplemented(
                " the cell is neither hydrogel nor stokes! Please set the correct material ids"));
    }

  hp::FEValues<dim>    volume_fe_values(mapping_collection,
                                        volume_fe_collection,
                                        quadrature_collection,
                                        update_values | update_JxW_values | update_quadrature_points);


  for(const auto &cell : volume_dof_handler.active_cell_iterators())
    if(cell->material_id() == hydrogel_domain_id)
    {
      volume_fe_values.reinit(cell);
      const FEValues<dim> &volume_fe_values =
        volume_fe_values.get_present_fe_values();
    }


}

int
main()
{
  try
    {
      test("extension_triangle_mesh_1.msh");
      test("extension_triangle_mesh.msh");
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}

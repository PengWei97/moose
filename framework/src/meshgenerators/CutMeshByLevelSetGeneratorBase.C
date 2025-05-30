//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CutMeshByLevelSetGeneratorBase.h"
#include "MooseMeshElementConversionUtils.h"
#include "MooseMeshUtils.h"

#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/cell_tet4.h"

// C++ includes
#include <cmath>

InputParameters
CutMeshByLevelSetGeneratorBase::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params += FunctionParserUtils<false>::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The input mesh that needs to be trimmed.");

  params.addParam<boundary_id_type>("cut_face_id",
                                    "The boundary id of the face generated by the cut. An "
                                    "id will be automatically assigned if not provided.");
  params.addParam<BoundaryName>(
      "cut_face_name", BoundaryName(), "The boundary name of the face generated by the cut.");

  params.addClassDescription(
      "This CutMeshByLevelSetGeneratorBase object is designed to be the base class of mesh "
      "generator that cuts a 3D mesh based on an analytic level set function. The level set "
      "function could be provided explicitly or indirectly.");

  return params;
}

CutMeshByLevelSetGeneratorBase::CutMeshByLevelSetGeneratorBase(const InputParameters & parameters)
  : MeshGenerator(parameters),
    FunctionParserUtils<false>(parameters),
    _input_name(getParam<MeshGeneratorName>("input")),
    _cut_face_name(getParam<BoundaryName>("cut_face_name")),
    _input(getMeshByName(_input_name))
{
  _cut_face_id = isParamValid("cut_face_id") ? getParam<boundary_id_type>("cut_face_id") : -1;
}

std::unique_ptr<MeshBase>
CutMeshByLevelSetGeneratorBase::generate()
{
  auto replicated_mesh_ptr = dynamic_cast<ReplicatedMesh *>(_input.get());
  if (!replicated_mesh_ptr)
    paramError("input", "Input is not a replicated mesh, which is required");
  if (*(replicated_mesh_ptr->elem_dimensions().begin()) != 3 ||
      *(replicated_mesh_ptr->elem_dimensions().rbegin()) != 3)
    paramError(
        "input",
        "Only 3D meshes are supported. Use XYMeshLineCutter for 2D meshes if applicable. Mixed "
        "dimensional meshes are not supported at the moment.");

  ReplicatedMesh & mesh = *replicated_mesh_ptr;

  if (!_cut_face_name.empty())
  {
    if (MooseMeshUtils::hasBoundaryName(mesh, _cut_face_name))
    {
      const boundary_id_type exist_cut_face_id =
          MooseMeshUtils::getBoundaryID(_cut_face_name, mesh);
      if (_cut_face_id != -1 && _cut_face_id != exist_cut_face_id)
        paramError("cut_face_id",
                   "The cut face boundary name and id are both provided, but they are inconsistent "
                   "with an existing boundary name which has a different id.");
      else
        _cut_face_id = exist_cut_face_id;
    }
    else
    {
      if (_cut_face_id == -1)
        _cut_face_id = MooseMeshUtils::getNextFreeBoundaryID(mesh);
      mesh.get_boundary_info().sideset_name(_cut_face_id) = _cut_face_name;
    }
  }
  else if (_cut_face_id == -1)
  {
    _cut_face_id = MooseMeshUtils::getNextFreeBoundaryID(mesh);
  }

  // Subdomain ID for new utility blocks must be new
  std::set<subdomain_id_type> subdomain_ids_set;
  mesh.subdomain_ids(subdomain_ids_set);
  const subdomain_id_type max_subdomain_id = *subdomain_ids_set.rbegin();
  const subdomain_id_type block_id_to_remove = max_subdomain_id + 1;

  // For the boolean value in the pair, true means the element is crossed by the cutting plane
  // false means the element is on the remaining side
  std::vector<std::pair<dof_id_type, bool>> cross_and_remained_elems_pre_convert;

  for (auto elem_it = mesh.active_elements_begin(); elem_it != mesh.active_elements_end();
       elem_it++)
  {
    const unsigned int & n_vertices = (*elem_it)->n_vertices();
    unsigned short elem_vertices_counter(0);
    for (unsigned int i = 0; i < n_vertices; i++)
    {
      // We define elem_vertices_counter in this way so that those elements with one face on the
      // plane are also processed to have the cut face boundary assigned.
      if (pointLevelSetRelation((*(*elem_it)->node_ptr(i))) !=
          PointLevelSetRelationIndex::level_set_in_side)
        elem_vertices_counter++;
    }
    if (elem_vertices_counter == n_vertices)
      (*elem_it)->subdomain_id() = block_id_to_remove;
    else
    {
      // Check if any elements to be processed are not first order
      if ((*elem_it)->default_order() != Order::FIRST)
        mooseError("Only first order elements are supported for cutting.");
      cross_and_remained_elems_pre_convert.push_back(
          std::make_pair((*elem_it)->id(), elem_vertices_counter > 0));
    }
  }

  std::vector<dof_id_type> converted_elems_ids_to_cut;
  // Then convert these non-TET4 elements into TET4 elements
  MooseMeshElementConversionUtils::convert3DMeshToAllTet4(mesh,
                                                          cross_and_remained_elems_pre_convert,
                                                          converted_elems_ids_to_cut,
                                                          block_id_to_remove,
                                                          false);

  std::vector<const Node *> new_on_plane_nodes;
  // We build the sideset information now, as we only need the information of the elements before
  // cutting
  BoundaryInfo & boundary_info = mesh.get_boundary_info();
  const auto bdry_side_list = boundary_info.build_side_list();
  // Cut the TET4 Elements
  for (const auto & converted_elems_id_to_cut : converted_elems_ids_to_cut)
  {
    tet4ElemCutter(
        mesh, bdry_side_list, converted_elems_id_to_cut, block_id_to_remove, new_on_plane_nodes);
  }

  // Delete the block to remove
  for (auto elem_it = mesh.active_subdomain_elements_begin(block_id_to_remove);
       elem_it != mesh.active_subdomain_elements_end(block_id_to_remove);
       elem_it++)
    mesh.delete_elem(*elem_it);

  mesh.contract();
  mesh.set_isnt_prepared();
  return std::move(_input);
}

CutMeshByLevelSetGeneratorBase::PointLevelSetRelationIndex
CutMeshByLevelSetGeneratorBase::pointLevelSetRelation(const Point & point)
{
  const Real level_set_value = levelSetEvaluator(point);
  if (MooseUtils::absoluteFuzzyLessThan(level_set_value, 0.0))
    return PointLevelSetRelationIndex::level_set_in_side;
  else if (MooseUtils::absoluteFuzzyGreaterThan(level_set_value, 0.0))
    return PointLevelSetRelationIndex::level_set_out_side;
  else
    return PointLevelSetRelationIndex::on_level_set;
}

Point
CutMeshByLevelSetGeneratorBase::pointPairLevelSetInterception(const Point & point1,
                                                              const Point & point2)
{
  Real dist1 = levelSetEvaluator(point1);
  Real dist2 = levelSetEvaluator(point2);

  if (MooseUtils::absoluteFuzzyEqual(dist1, 0.0) || MooseUtils::absoluteFuzzyEqual(dist2, 0.0))
    mooseError("At least one of the two points are on the plane.");
  if (MooseUtils::absoluteFuzzyGreaterThan(dist1 * dist2, 0.0))
    mooseError("The two points are on the same side of the plane.");

  Point p1(point1);
  Point p2(point2);
  Real dist = abs(dist1) + abs(dist2);
  Point mid_point;

  // Bisection method to find midpoint
  while (MooseUtils::absoluteFuzzyGreaterThan(dist, 0.0))
  {
    mid_point = 0.5 * (p1 + p2);
    const Real dist_mid = levelSetEvaluator(mid_point);
    // We do not need Fuzzy here as it will be covered by the while loop
    if (dist_mid == 0.0)
      dist = 0.0;
    else if (dist_mid * dist1 < 0.0)
    {
      p2 = mid_point;
      dist2 = levelSetEvaluator(p2);
      dist = abs(dist1) + abs(dist2);
    }
    else
    {
      p1 = mid_point;
      dist1 = levelSetEvaluator(p1);
      dist = abs(dist1) + abs(dist2);
    }
  }
  return mid_point;
}

const Node *
CutMeshByLevelSetGeneratorBase::nonDuplicateNodeCreator(
    ReplicatedMesh & mesh,
    std::vector<const Node *> & new_on_plane_nodes,
    const Point & new_point) const
{
  for (const auto & new_on_plane_node : new_on_plane_nodes)
  {
    if (MooseUtils::absoluteFuzzyEqual((*new_on_plane_node - new_point).norm(), 0.0))
      return new_on_plane_node;
  }
  new_on_plane_nodes.push_back(mesh.add_point(new_point));
  return new_on_plane_nodes.back();
}

void
CutMeshByLevelSetGeneratorBase::tet4ElemCutter(
    ReplicatedMesh & mesh,
    const std::vector<libMesh::BoundaryInfo::BCTuple> & bdry_side_list,
    const dof_id_type elem_id,
    const subdomain_id_type & block_id_to_remove,
    std::vector<const Node *> & new_on_plane_nodes)
{
  // Retrieve boundary information for the mesh
  BoundaryInfo & boundary_info = mesh.get_boundary_info();
  // Create a list of sidesets involving the element to be split
  // It might be complex to assign the boundary id to the new elements
  // In TET4, none of the four faces have the same normal vector
  // So we are using the normal vector to identify the faces that
  // need to be assigned the same boundary id
  std::vector<Point> elem_side_normal_list;
  elem_side_normal_list.resize(4);
  for (const auto i : make_range(4))
  {
    auto elem_side = mesh.elem_ptr(elem_id)->side_ptr(i);
    elem_side_normal_list[i] = (*elem_side->node_ptr(1) - *elem_side->node_ptr(0))
                                   .cross(*elem_side->node_ptr(2) - *elem_side->node_ptr(1))
                                   .unit();
  }

  std::vector<std::vector<boundary_id_type>> elem_side_list;
  MooseMeshElementConversionUtils::elementBoundaryInfoCollector(
      bdry_side_list, elem_id, 4, elem_side_list);

  std::vector<PointLevelSetRelationIndex> node_plane_relation(4);
  std::vector<const Node *> tet4_nodes(4);
  std::vector<const Node *> tet4_nodes_on_plane;
  std::vector<const Node *> tet4_nodes_outside_plane;
  std::vector<const Node *> tet4_nodes_inside_plane;
  // Sort tetrahedral nodes based on their positioning wrt the plane
  for (unsigned int i = 0; i < 4; i++)
  {
    tet4_nodes[i] = mesh.elem_ptr(elem_id)->node_ptr(i);
    node_plane_relation[i] = pointLevelSetRelation(*tet4_nodes[i]);
    if (node_plane_relation[i] == PointLevelSetRelationIndex::on_level_set)
      tet4_nodes_on_plane.push_back(tet4_nodes[i]);
    else if (node_plane_relation[i] == PointLevelSetRelationIndex::level_set_out_side)
      tet4_nodes_outside_plane.push_back(tet4_nodes[i]);
    else
      tet4_nodes_inside_plane.push_back(tet4_nodes[i]);
  }
  std::vector<Elem *> elems_tet4;
  bool new_elements_created(false);
  // No cutting needed if no nodes are outside the plane
  if (tet4_nodes_outside_plane.size() == 0)
  {
    if (tet4_nodes_on_plane.size() == 3)
    {
      // Record the element for future cross section boundary assignment
      elems_tet4.push_back(mesh.elem_ptr(elem_id));
    }
  }
  // Remove the element if all the nodes are outside the plane
  else if (tet4_nodes_inside_plane.size() == 0)
  {
    mesh.elem_ptr(elem_id)->subdomain_id() = block_id_to_remove;
    if (tet4_nodes_on_plane.size() == 3)
    {
      // I think the neighboring element will be handled,
      // So we do not need to assign the cross section boundary here
    }
  }
  // As we have nodes on both sides, six different scenarios are possible
  else
  {
    new_elements_created = true;
    if (tet4_nodes_inside_plane.size() == 1 && tet4_nodes_outside_plane.size() == 3)
    {
      std::vector<const Node *> new_plane_nodes;
      // A smaller TET4 element is created, this solution is unique
      for (const auto & tet4_node_outside_plane : tet4_nodes_outside_plane)
      {
        new_plane_nodes.push_back(nonDuplicateNodeCreator(
            mesh,
            new_on_plane_nodes,
            pointPairLevelSetInterception(*tet4_node_outside_plane, *tet4_nodes_inside_plane[0])));
      }
      auto new_elem_tet4 = std::make_unique<Tet4>();
      new_elem_tet4->set_node(0, const_cast<Node *>(tet4_nodes_inside_plane[0]));
      new_elem_tet4->set_node(1, const_cast<Node *>(new_plane_nodes[0]));
      new_elem_tet4->set_node(2, const_cast<Node *>(new_plane_nodes[1]));
      new_elem_tet4->set_node(3, const_cast<Node *>(new_plane_nodes[2]));
      new_elem_tet4->subdomain_id() = mesh.elem_ptr(elem_id)->subdomain_id();
      elems_tet4.push_back(mesh.add_elem(std::move(new_elem_tet4)));
    }
    else if (tet4_nodes_inside_plane.size() == 2 && tet4_nodes_outside_plane.size() == 2)
    {
      std::vector<const Node *> new_plane_nodes;
      // 3 smaller TET3 elements are created
      for (const auto & tet4_node_outside_plane : tet4_nodes_outside_plane)
      {
        for (const auto & tet4_node_inside_plane : tet4_nodes_inside_plane)
        {
          new_plane_nodes.push_back(nonDuplicateNodeCreator(
              mesh,
              new_on_plane_nodes,
              pointPairLevelSetInterception(*tet4_node_outside_plane, *tet4_node_inside_plane)));
        }
      }
      std::vector<const Node *> new_elems_nodes = {tet4_nodes_inside_plane[1],
                                                   new_plane_nodes[3],
                                                   new_plane_nodes[1],
                                                   tet4_nodes_inside_plane[0],
                                                   new_plane_nodes[2],
                                                   new_plane_nodes[0]};
      std::vector<std::vector<unsigned int>> rotated_tet_face_indices;
      std::vector<std::vector<const Node *>> optimized_node_list;
      MooseMeshElementConversionUtils::prismNodesToTetNodesDeterminer(
          new_elems_nodes, rotated_tet_face_indices, optimized_node_list);

      for (unsigned int i = 0; i < optimized_node_list.size(); i++)
      {
        auto new_elem_tet4 = std::make_unique<Tet4>();
        new_elem_tet4->set_node(0, const_cast<Node *>(optimized_node_list[i][0]));
        new_elem_tet4->set_node(1, const_cast<Node *>(optimized_node_list[i][1]));
        new_elem_tet4->set_node(2, const_cast<Node *>(optimized_node_list[i][2]));
        new_elem_tet4->set_node(3, const_cast<Node *>(optimized_node_list[i][3]));
        new_elem_tet4->subdomain_id() = mesh.elem_ptr(elem_id)->subdomain_id();
        elems_tet4.push_back(mesh.add_elem(std::move(new_elem_tet4)));
      }
    }
    else if (tet4_nodes_inside_plane.size() == 3 && tet4_nodes_outside_plane.size() == 1)
    {
      std::vector<const Node *> new_plane_nodes;
      // 3 smaller Tet4 elements are created
      for (const auto & tet4_node_inside_plane : tet4_nodes_inside_plane)
      {
        new_plane_nodes.push_back(nonDuplicateNodeCreator(
            mesh,
            new_on_plane_nodes,
            pointPairLevelSetInterception(*tet4_node_inside_plane, *tet4_nodes_outside_plane[0])));
      }
      std::vector<const Node *> new_elems_nodes = {tet4_nodes_inside_plane[0],
                                                   tet4_nodes_inside_plane[1],
                                                   tet4_nodes_inside_plane[2],
                                                   new_plane_nodes[0],
                                                   new_plane_nodes[1],
                                                   new_plane_nodes[2]};
      std::vector<std::vector<unsigned int>> rotated_tet_face_indices;
      std::vector<std::vector<const Node *>> optimized_node_list;
      MooseMeshElementConversionUtils::prismNodesToTetNodesDeterminer(
          new_elems_nodes, rotated_tet_face_indices, optimized_node_list);

      for (unsigned int i = 0; i < optimized_node_list.size(); i++)
      {
        auto new_elem_tet4 = std::make_unique<Tet4>();
        new_elem_tet4->set_node(0, const_cast<Node *>(optimized_node_list[i][0]));
        new_elem_tet4->set_node(1, const_cast<Node *>(optimized_node_list[i][1]));
        new_elem_tet4->set_node(2, const_cast<Node *>(optimized_node_list[i][2]));
        new_elem_tet4->set_node(3, const_cast<Node *>(optimized_node_list[i][3]));
        new_elem_tet4->subdomain_id() = mesh.elem_ptr(elem_id)->subdomain_id();
        elems_tet4.push_back(mesh.add_elem(std::move(new_elem_tet4)));
      }
    }
    else if (tet4_nodes_inside_plane.size() == 1 && tet4_nodes_outside_plane.size() == 1)
    {
      auto new_plane_node = nonDuplicateNodeCreator(
          mesh,
          new_on_plane_nodes,
          pointPairLevelSetInterception(*tet4_nodes_inside_plane[0], *tet4_nodes_outside_plane[0]));
      // A smaller Tet4 is created, this solution is unique
      auto new_elem_tet4 = std::make_unique<Tet4>();
      new_elem_tet4->set_node(0, const_cast<Node *>(new_plane_node));
      new_elem_tet4->set_node(1, const_cast<Node *>(tet4_nodes_on_plane[0]));
      new_elem_tet4->set_node(2, const_cast<Node *>(tet4_nodes_on_plane[1]));
      new_elem_tet4->set_node(3, const_cast<Node *>(tet4_nodes_inside_plane[0]));
      new_elem_tet4->subdomain_id() = mesh.elem_ptr(elem_id)->subdomain_id();
      elems_tet4.push_back(mesh.add_elem(std::move(new_elem_tet4)));
    }
    else if (tet4_nodes_inside_plane.size() == 1 && tet4_nodes_outside_plane.size() == 2)
    {
      std::vector<const Node *> new_plane_nodes;
      // A smaller Tet4 element is created, this solution is unique
      for (const auto & tet4_node_outside_plane : tet4_nodes_outside_plane)
      {
        new_plane_nodes.push_back(nonDuplicateNodeCreator(
            mesh,
            new_on_plane_nodes,
            pointPairLevelSetInterception(*tet4_node_outside_plane, *tet4_nodes_inside_plane[0])));
      }
      auto new_elem_tet4 = std::make_unique<Tet4>();
      new_elem_tet4->set_node(0, const_cast<Node *>(new_plane_nodes[0]));
      new_elem_tet4->set_node(1, const_cast<Node *>(new_plane_nodes[1]));
      new_elem_tet4->set_node(2, const_cast<Node *>(tet4_nodes_on_plane[0]));
      new_elem_tet4->set_node(3, const_cast<Node *>(tet4_nodes_inside_plane[0]));
      new_elem_tet4->subdomain_id() = mesh.elem_ptr(elem_id)->subdomain_id();
      elems_tet4.push_back(mesh.add_elem(std::move(new_elem_tet4)));
    }
    else if (tet4_nodes_inside_plane.size() == 2 && tet4_nodes_outside_plane.size() == 1)
    {
      std::vector<const Node *> new_plane_nodes;
      // 2 smaller TET4 elements are created
      for (const auto & tet4_node_inside_plane : tet4_nodes_inside_plane)
      {
        new_plane_nodes.push_back(nonDuplicateNodeCreator(
            mesh,
            new_on_plane_nodes,
            pointPairLevelSetInterception(*tet4_node_inside_plane, *tet4_nodes_outside_plane[0])));
      }
      std::vector<const Node *> new_elems_nodes = {tet4_nodes_inside_plane[0],
                                                   tet4_nodes_inside_plane[1],
                                                   new_plane_nodes[1],
                                                   new_plane_nodes[0],
                                                   tet4_nodes_on_plane[0]};
      std::vector<std::vector<unsigned int>> rotated_tet_face_indices;
      std::vector<std::vector<const Node *>> optimized_node_list;
      MooseMeshElementConversionUtils::pyramidNodesToTetNodesDeterminer(
          new_elems_nodes, rotated_tet_face_indices, optimized_node_list);

      for (unsigned int i = 0; i < optimized_node_list.size(); i++)
      {
        auto new_elem_tet4 = std::make_unique<Tet4>();
        new_elem_tet4->set_node(0, const_cast<Node *>(optimized_node_list[i][0]));
        new_elem_tet4->set_node(1, const_cast<Node *>(optimized_node_list[i][1]));
        new_elem_tet4->set_node(2, const_cast<Node *>(optimized_node_list[i][2]));
        new_elem_tet4->set_node(3, const_cast<Node *>(optimized_node_list[i][3]));
        new_elem_tet4->subdomain_id() = mesh.elem_ptr(elem_id)->subdomain_id();
        elems_tet4.push_back(mesh.add_elem(std::move(new_elem_tet4)));
      }
    }
    else
      mooseError("Unexpected scenario.");

    mesh.elem_ptr(elem_id)->subdomain_id() = block_id_to_remove;
  }

  for (auto & elem_tet4 : elems_tet4)
  {
    if (new_elements_created)
    {
      if (elem_tet4->volume() < 0.0)
      {
        Node * temp = elem_tet4->node_ptr(0);
        elem_tet4->set_node(0, elem_tet4->node_ptr(1));
        elem_tet4->set_node(1, temp);
      }
    }
    // Find the boundary id of the new element
    for (unsigned int i = 0; i < 4; i++)
    {
      const Point & side_pt_0 = *elem_tet4->side_ptr(i)->node_ptr(0);
      const Point & side_pt_1 = *elem_tet4->side_ptr(i)->node_ptr(1);
      const Point & side_pt_2 = *elem_tet4->side_ptr(i)->node_ptr(2);

      const Point side_normal = (side_pt_1 - side_pt_0).cross(side_pt_2 - side_pt_1).unit();
      for (unsigned int j = 0; j < 4; j++)
      {
        if (new_elements_created)
        {
          if (MooseUtils::absoluteFuzzyEqual(side_normal * elem_side_normal_list[j], 1.0))
          {
            for (const auto & side_info : elem_side_list[j])
            {
              boundary_info.add_side(elem_tet4, i, side_info);
            }
          }
        }
      }
      if (MooseUtils::absoluteFuzzyEqual(levelSetEvaluator(side_pt_0), 0.0) &&
          MooseUtils::absoluteFuzzyEqual(levelSetEvaluator(side_pt_1), 0.0) &&
          MooseUtils::absoluteFuzzyEqual(levelSetEvaluator(side_pt_2), 0.0))
      {

        boundary_info.add_side(elem_tet4, i, _cut_face_id);
      }
    }
  }
}

Real
CutMeshByLevelSetGeneratorBase::levelSetEvaluator(const Point & point)
{
  _func_params[0] = point(0);
  _func_params[1] = point(1);
  _func_params[2] = point(2);
  return evaluate(_func_level_set);
}

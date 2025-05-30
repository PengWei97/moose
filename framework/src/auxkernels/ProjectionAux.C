//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ProjectionAux.h"
#include "SystemBase.h"
#include "libmesh/system.h"

registerMooseObjectRenamed("MooseApp", SelfAux, "01/30/2024 24:00", ProjectionAux);
registerMooseObject("MooseApp", ProjectionAux);

InputParameters
ProjectionAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Returns the specified variable as an auxiliary variable with a projection of the source "
      "variable. If they are the same type, this amounts to a simple copy.");
  params.addRequiredCoupledVar("v", "Variable to take the value of.");

  params.addParam<bool>("use_block_restriction_for_source",
                        false,
                        "Whether to use the auxkernel block restriction to also restrict the "
                        "locations selected for source variable values");

  // Technically possible to project from nodal to elemental and back
  params.set<bool>("_allow_nodal_to_elemental_coupling") = true;

  // We need some ghosting for all elemental to nodal projections
  params.addParam<unsigned short>("ghost_layers", 1, "The number of layers of elements to ghost.");
  params.addRelationshipManager("ElementPointNeighborLayers",
                                Moose::RelationshipManagerType::ALGEBRAIC,
                                [](const InputParameters & obj_params, InputParameters & rm_params)
                                {
                                  rm_params.set<unsigned short>("layers") =
                                      obj_params.get<unsigned short>("ghost_layers");
                                  rm_params.set<bool>("use_displaced_mesh") =
                                      obj_params.get<bool>("use_displaced_mesh");
                                });
  return params;
}

ProjectionAux::ProjectionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _v(coupledValue("v")),
    _source_variable(*getFieldVar("v", 0)),
    _source_sys(_c_fe_problem.getSystem(coupledName("v"))),
    _use_block_restriction_for_source(getParam<bool>("use_block_restriction_for_source"))
{
  // Output some messages to user
  if (_source_variable.order() > _var.order())
    mooseInfo("Projection lowers order, please expect a loss of accuracy");
  // Check the dimension of the block restriction
  for (const auto & sub_id : blockIDs())
    if (_mesh.isLowerD(sub_id))
      paramError("block",
                 "ProjectionAux's block restriction must not include lower dimensional blocks");
}

Real
ProjectionAux::computeValue()
{
  if (!isNodal() || (_source_variable.isNodal() && _source_variable.order() >= _var.order()))
    return _v[_qp];
  // projecting continuous elemental variable onto a nodal one
  // AND nodal low order -> nodal higher order
  else if (isNodal() && _source_variable.getContinuity() != DISCONTINUOUS &&
           _source_variable.getContinuity() != SIDE_DISCONTINUOUS)
  {
    return _source_sys.point_value(
        _source_variable.number(), *_current_node, elemOnNodeVariableIsDefinedOn());
  }
  // Handle discontinuous elemental variable projection into a nodal variable
  else
  {
    // Custom projection rule : use nodal values, if discontinuous the one coming from each element,
    // weighted by element volumes
    // First, find all the elements that this node is part of
    auto elem_ids = _mesh.nodeToElemMap().find(_current_node->id());
    mooseAssert(elem_ids != _mesh.nodeToElemMap().end(),
                "Should have found an element around node " + std::to_string(_current_node->id()));

    // Get the neighbor element centroid values & element volumes
    Real sum_weighted_values = 0;
    Real sum_volumes = 0;
    for (auto & id : elem_ids->second)
    {
      const auto & elem = _mesh.elemPtr(id);
      const auto block_id = elem->subdomain_id();
      // Only use higher D elements
      // We allow full-dimensional elements in a higher dimension mesh
      if (_source_variable.hasBlocks(block_id) && (!_mesh.isLowerD(block_id)) &&
          (!_use_block_restriction_for_source || hasBlocks(block_id)))
      {
        const auto elem_volume = elem->volume();
        sum_weighted_values +=
            _source_sys.point_value(_source_variable.number(), *_current_node, elem) * elem_volume;
        sum_volumes += elem_volume;
      }
    }
    if (sum_volumes == 0)
      mooseError("Did not find a valid source variable value for node: ", *_current_node);
    return sum_weighted_values / sum_volumes;
  }
}

const Elem *
ProjectionAux::elemOnNodeVariableIsDefinedOn() const
{
  for (const auto & elem_id : _mesh.nodeToElemMap().find(_current_node->id())->second)
  {
    const auto & elem = _mesh.elemPtr(elem_id);
    const auto block_id = elem->subdomain_id();
    if (_source_variable.hasBlocks(block_id) &&
        (!_mesh.isLowerD(block_id) && elem->dim() == _mesh.dimension()) &&
        (!_use_block_restriction_for_source || hasBlocks(block_id)))
      return elem;
  }
  mooseError("Source variable is not defined everywhere the target variable is");
}

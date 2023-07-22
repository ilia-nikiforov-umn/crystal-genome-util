"""
Function for adding one instance of the ``binding-energy-crystal`` and ``crystal-structure-npt`` properties each to the
"""

from kim_property import kim_property_create, kim_property_modify
from typing import Dict, List, Union
import numpy as np
import kim_edn
from .aflow_util import get_stoich_reduced_list_from_prototype, AFLOW, read_shortnames
from .physics_validator import validate_binding_energy_crystal, validate_crystal_structure_npt
import os
from kim_property import kim_property_dump

ENERGY_PROPERTY_ID = "tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal"
STRUCTURE_PROPERTY_ID = "tag:staff@noreply.openkim.org,2023-02-21:property/crystal-structure-npt"

###############################################################################
def add_property_inst(
    energy_atom:float,stoichiometric_species:List[str],proto_des: Dict, libproto: str, shortname: str, property_instances: Union[str,None] = None
) -> str:
    """
    Add a pair of property instances to the property_instances object.

    Create the KIM property instances for this Test or Reference data. See the documentation for the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_  and the definition of the ``binding-energy-crystal`` and ``crystal-structure-npt`` properties
    
    Args:
        energy_atom:
            Energy per atom
        stoichiometric_species:
            Element symbols corresponding to the atom types in the stoichiometric formula which appears at the start of the prototype label (e.g. ['Mo','S'] for the AB2 stoichiometric formula, means that the 'A' atom is 'Mo' and the 'B' atom is 'S' for the MoS_2 structure).
        proto_des:
            AFLOW prototype designation
        libproto:
            AFLOW library prototype
        shortname:
            Material shortname
        property_instances:
            Existing property_instances object, if any
    Returns:        
        Updated property_instances object
    """
    
    # figure out required values
    prototype_label = proto_des["aflow_prototype_label"]
    a = proto_des["aflow_prototype_params_values"][0]
    n_formula = sum(get_stoich_reduced_list_from_prototype(prototype_label))

    # create energy property    
    if property_instances is None:
        energy_index = 1
    else:
        energy_index=len(kim_edn.loads(property_instances))+1
    property_instances = kim_property_create(
        energy_index,
        ENERGY_PROPERTY_ID,property_instances
    )

    structure_index=len(kim_edn.loads(property_instances))+1
    # create structure property
    property_instances = kim_property_create(
        structure_index,
        STRUCTURE_PROPERTY_ID,property_instances
    )    

    # Set all required key-value pairs for energy
    property_instances = kim_property_modify(
        property_instances,
        energy_index,
        #
        "key",
        "prototype-label",
        "source-value",
        prototype_label,
        #
        "key",
        "stoichiometric-species",
        "source-value",
        "1:{}".format(len(stoichiometric_species)),
        *stoichiometric_species,
        #
        "key",
        "a",
        "source-value",
        a,
        "source-unit",
        "angstrom",
        #
        "key",
        "binding-potential-energy-per-atom",
        "source-value",
        energy_atom,
        "source-unit",
        "eV",
        #
        "key",
        "binding-potential-energy-per-formula",
        "source-value",
        energy_atom*n_formula,
        "source-unit",
        "eV",
    )

    # Set all required key-value pairs for structure
    property_instances = kim_property_modify(
        property_instances,
        structure_index,
        #
        "key",
        "prototype-label",
        "source-value",
        prototype_label,
        #
        "key",
        "stoichiometric-species",
        "source-value",
        "1:{}".format(len(stoichiometric_species)),
        *stoichiometric_species,
        #
        "key",
        "a",
        "source-value",
        a,
        "source-unit",
        "angstrom",
        #
        "key",
        "cell-cauchy-stress",
        "source-value",
        "1:6",
        *np.zeros(6),
        "source-unit",
        "eV/angstrom^3",
        #
        "key",
        "temperature",
        "source-value",
        0.,
        "source-unit",
        "K",
    )

    # write non-a parameters if present
    if len(proto_des["aflow_prototype_params_values"])>1:
        parameter_names=proto_des["aflow_prototype_params_list"][1:]
        parameter_values=proto_des["aflow_prototype_params_values"][1:]
        for index in (energy_index,structure_index):
            property_instances = kim_property_modify(
                property_instances,
                index,
                "key",
                "parameter-names",
                "source-value",
                "1:{}".format(len(parameter_names)),
                *parameter_names,
                #
                "key",
                "parameter-values",
                "source-value",
                "1:{}".format(len(parameter_values)),
                *parameter_values,
            )

    if libproto is not None:
        for index in (energy_index,structure_index):        
            property_instances = kim_property_modify(
                property_instances,
                index,                         
                #   
                "key",
                "library-prototype-label",
                "source-value",
                libproto
            )

    if shortname is not None:
        for index in (energy_index,structure_index):        
            property_instances = kim_property_modify(
                property_instances,
                index,                         
                #   
                "key",
                "short-name",
                "source-value",
                shortname
            )                            
    
    return property_instances

def find_unique_materials_and_write_properties(compare_dir: str, prototype_label: str, energy_per_atom: List[float], species):
    """
    Given a directory of relaxed structures, find unique ones and write their properties.

    1) Use :func:`util.aflow_util.AFLOW.compare_materials_dir` to identify groups of duplicates among relaxed structures
    2) For each group of duplicates:
        a) Loop over the inidividual files, try to add a property instance until one is successfully added:
            i) Use :func:`.aflow_util.AFLOW.get_prototype` to check that the prototype label didn't change
            ii) Get the library prototype and shortname using :func:`.aflow_util.AFLOW.get_library_prototype_label_and_shortname`
            iii) Try to add the properties using :func:`add_property_inst`
            iv) Validate the properties using :func:`.physics_validator.validate_binding_energy_crystal` and :func:`.physics_validator.validate_crystal_structure_npt`
    3) If at least one property instance was added, write the ``output/results.edn`` file

    Args:
        compare_dir:
            Directory of structures to compare. This should contain filenames that are just integers corresponding to indices
        prototype_label:
            Prototype label of structures. This must be unchanged for property to write.
        energy_per_atom:
            List of relaxed energies per atom corresponding to the relaxed structures
        species:
            Chemical species
            
    """
    aflow = AFLOW()
    shortnames = read_shortnames()

    # compare relaxed structures for duplicates
    comparison=aflow.compare_materials_dir(compare_dir)

    # Try to write at most one property per group of duplicates
    property_inst = None
    for materials_group in comparison:
        # get all "indices" belonging to this group, i.e. numerical filenames
        indices = [int(materials_group['structure_representative']['name'].split('/')[-1])]
        for structure in materials_group['structures_duplicate']:
            indices.append(int(structure['name'].split('/')[-1]))
        print ("Parameter sets %s relaxed to duplicate structures, attempting to write only one of them."%str(indices))

        # loop over all materials in this group until we find one that successfully adds a property instance
        for i in indices:
            relax_poscar_path=os.path.join(compare_dir,str(i))
            relax_proto_des=aflow.get_prototype(relax_poscar_path)
            relax_proto_label=relax_proto_des['aflow_prototype_label']
            if relax_proto_label != prototype_label:
                print("Prototype label changed during relaxation: test template prototype is %s, while relaxed is %s. Skipping parameter set %d."
                %(prototype_label,relax_proto_label,i))
                if i==indices[-1]:
                    print("No parameter sets in this group successfully added a property instance. Skipping this group.")
                continue
            (libproto,shortname)=aflow.get_library_prototype_label_and_shortname(relax_poscar_path,shortnames)
            # Try to add property instances. Back up original so we can keep going even if one parameter set fails
            property_inst_new = property_inst # these are just serialized edn strings, so no need for deepcopy or anything
            try:
                property_inst_new = add_property_inst(energy_per_atom[i],species,relax_proto_des,libproto,shortname,property_inst_new)
                validate_crystal_structure_npt(property_inst_new,STRUCTURE_PROPERTY_ID)
                validate_binding_energy_crystal(property_inst_new,ENERGY_PROPERTY_ID)
                property_inst = property_inst_new
                print("Successfully added property instance for parameter set %d"%i)
                break
            except Exception as e:
                print("Skipping parameter set %d because of error while adding or validating property:\n %s"%(i,e))
                if i==indices[-1]:
                    print("No parameter sets in this group successfully added a property instance. Skipping this group.")
            
    # Dump the results in a file
    if not (property_inst is None or \
            property_inst in ('None', '', '[]')):
        with open("output/results.edn", "w") as fp:
            kim_property_dump(property_inst, fp)
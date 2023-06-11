"""
Function for adding one instance of the ``binding-energy-crystal`` and ``crystal-structure-npt`` properties each to the
"""

from kim_property import kim_property_create, kim_property_modify
from typing import Dict, List, Union
import numpy as np
import kim_edn
from .aflow_util import get_stoich_reduced_list_from_prototype

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
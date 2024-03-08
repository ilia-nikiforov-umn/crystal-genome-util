"""
Functions for writing and validating common fields of Crystal Genome properties
"""

from kim_property import kim_property_modify
from typing import Dict, List, Union
import numpy as np
from ..aflow_util import AFLOW
from curses.ascii import isdigit
import subprocess


ALLOWED_FULL_PARAMETER_NAMES = ['b/a','c/a','alpha','beta','gamma']
XYZ = ['x','y','z']
ELEMENTS = [
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
        "Cn",
        "Uut",
        "Fl",
        "Uup",
        "Lv",
        "Uus",
        "Uuo",
    ]

def convert_units(value: float, unit_from: str, unit_to: str)->float:
    """
    Invoke the GNU units utility

    """
    return float(subprocess.check_output("units -t \"%s(%f)\" %s"%(unit_from,value,unit_to),shell=True).strip())

def validate_common_fields(property_inst: Dict):
    """
    Validate common fields of Crystal Genome properties

    Args:
        property_inst: A python dictionary containing a single property instance
    Raises:
        RuntimeError: When validation fails
    """
    for species in property_inst['stoichiometric-species']['source-value']:
        if species not in ELEMENTS:
            raise RuntimeError("Unphysical species %s reported"%species)
    
    if property_inst["a"]["source-value"] < 0:
        raise RuntimeError("Negative lattice constant")

    free_params = [property_inst["a"]["source-value"]]

    if "parameter-names" in property_inst:
        if len(property_inst["parameter-names"]["source-value"]) != len(property_inst["parameter-values"]["source-value"]):
            raise RuntimeError("Length of parameter names does not match length of parameter values")
        for i,name in enumerate(property_inst["parameter-names"]["source-value"]):
            parameter_value = property_inst["parameter-values"]["source-value"][i]
            if name not in ALLOWED_FULL_PARAMETER_NAMES:
                if (name[0] in XYZ) and (len(name) > 1):
                    for char in name[1:]:
                        if not isdigit(char):
                            raise RuntimeError("Illegal parameter name %s"%name)
                    # This parameter is a fractional coordinate
                    if not (0<=parameter_value<=1):
                        raise RuntimeError("Parameter %s has illegal value %f" %(name,parameter_value))        
                else:
                    raise RuntimeError("Illegal parameter name %s"%name)
            # This parameter is a ratio or angle. Just needs to be a positive float
            if parameter_value < 0:
                raise RuntimeError("Parameter %s has illegal value %f" %(name,parameter_value))
            free_params.append(parameter_value)

    aflow = AFLOW()
    aflow.write_poscar(property_inst["prototype-label"]["source-value"],free_params=free_params)

    if "library-prototype-label" in property_inst:
        aflow.write_poscar(property_inst["library-prototype-label"]["source-value"])

    if "temperature" in property_inst:
        # GNU Units will crash if the absolute temperature is negative
        convert_units(property_inst["temperature"]["source-value"],property_inst["temperature"]["source-unit"],"K")
    
    if "cell-cauchy-stress" in property_inst:
        stress={
            "xx":property_inst["cell-cauchy-stress"]["source-value"][0],
            "yy":property_inst["cell-cauchy-stress"]["source-value"][1],
            "zz":property_inst["cell-cauchy-stress"]["source-value"][2],
            "yz":property_inst["cell-cauchy-stress"]["source-value"][3],
            "xz":property_inst["cell-cauchy-stress"]["source-value"][4],
            "xy":property_inst["cell-cauchy-stress"]["source-value"][5]
        }
        if "parameter-names" not in property_inst:
            # Cubic crystal, stress must be isotropic
            if not((stress["xx"]==stress["yy"]==stress["zz"]) and (stress["yz"]==stress["xz"]==stress["xy"]==0.)):
                raise RuntimeError ("Non-isotropic stress state "+str(stress)+" in cubic crystal")
        else:
            if "b/a" not in property_inst["parameter-names"]["source-value"]:
                if stress["xx"] != stress["yy"]:
                    raise RuntimeError(
                        "Stress component xx %f is not equal to stress component yy %f. This is incompatible with a crystal prototype lacking parameter b/a"%
                        (stress["xx"], stress["yy"])
                        )
            if "c/a" not in property_inst["parameter-names"]["source-value"]:
                if stress["xx"] != stress["zz"]:
                    raise RuntimeError(
                        "Stress component xx %f is not equal to stress component zz %f. This is incompatible with a crystal prototype lacking parameter c/a"%
                        (stress["xx"], stress["zz"])
                        )
            if "alpha" not in property_inst["parameter-names"]["source-value"]:
                if stress["yz"] != 0.:
                    raise RuntimeError(
                        "Stress component yz %f is nonzero. This is incompatible with a crystal prototype lacking parameter alpha"%
                        stress["yz"]
                        )
            if "beta" not in property_inst["parameter-names"]["source-value"]:
                if stress["xz"] != 0.:
                    raise RuntimeError(
                        "Stress component xz %f is nonzero. This is incompatible with a crystal prototype lacking parameter beta"%
                        stress["xz"]
                        )
            if "gamma" not in property_inst["parameter-names"]["source-value"]:
                if stress["xy"] != 0.:
                    raise RuntimeError(
                        "Stress component xy %f is nonzero. This is incompatible with a crystal prototype lacking parameter gamma"%
                        stress["xy"]
                        )                                                        

def add_common_fields(property_instances: str, property_index: int, stoichiometric_species:List[str], proto_des: Dict, libproto: str, shortname: str, 
                      stress: Union[np.ndarray,None] = None, stress_unit: Union[str,None] = None, temperature: Union[float,None] = None) -> str:
    """
    Add common crystal-genome fields to an initialized property

    Args:
        property_instances: Serialized kim-edn property instances object
        property_index: index of the property we are writing to in the instances object
        stoichiometric_species:
            Element symbols corresponding to the atom types in the stoichiometric formula which appears at the start of the prototype label (e.g. ['Mo','S'] for the AB2 stoichiometric formula, 
            means that the 'A' atom is 'Mo' and the 'B' atom is 'S' for the MoS_2 structure).
        proto_des: AFLOW prototype designation from crystal_genome_util.aflow_util.AFLOW.get_prototype
        libproto: AFLOW library prototype
        shortname: Material shortname
        stress: Cauchy stress
        stress_unit: unit of stress
        temperature: temperature in K
    Returns:
        Updated property_instances object
    """
    
    prototype_label = proto_des["aflow_prototype_label"]
    a = proto_des["aflow_prototype_params_values"][0]
    
    property_instances = kim_property_modify(
        property_instances,
        property_index,
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
    )

    # write non-a parameters if present
    if len(proto_des["aflow_prototype_params_values"])>1:
        parameter_names=proto_des["aflow_prototype_params_list"][1:]
        parameter_values=proto_des["aflow_prototype_params_values"][1:]
        property_instances = kim_property_modify(
            property_instances,
            property_index,
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
        property_instances = kim_property_modify(
            property_instances,
            property_index,                         
            #   
            "key",
            "library-prototype-label",
            "source-value",
            libproto
        )

    if shortname is not None:
        property_instances = kim_property_modify(
            property_instances,
            property_index,                         
            #   
            "key",
            "short-name",
            "source-value",
            "1:1",
            shortname
        )                            

    if stress is not None:
        property_instances = kim_property_modify(
            property_instances,
            property_index,             
            #
            "key",
            "cell-cauchy-stress",
            "source-value",
            "1:6",
            *stress,
            "source-unit",
            stress_unit,
        )                            

    if temperature is not None:
        property_instances = kim_property_modify(
            property_instances,
            property_index,                         
            #
            "key",
            "temperature",
            "source-value",
            temperature,
            "source-unit",
            "K",
        )

    return property_instances
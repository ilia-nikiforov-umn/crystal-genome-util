"""
Python functions for validating binding-energy-crystal and crystal-structure-npt properties
"""
from typing import Any, Callable, Dict, List, Optional, Set
from curses.ascii import isdigit
from .aflow_util import AFLOW, get_stoich_reduced_list_from_prototype
import kim_edn
import subprocess

EV_ENERGY_EPSILON = 1e-3 #For checking whether energy per formula is a correct multiple of energy per atom. Don't want to be strict here, so just use 1 meV/atom
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
    The two properties this package generates have many common fields. Validate them in this function.

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
                    
def validate_binding_energy_crystal(property_instances: str, energy_property_id: str):
    """
    Validate "binding-energy-crystal" property

    Args:
        property_instances: Serialized kim-edn property instances object
        energy_property_id: name of the property
    Raises:
        RuntimeError: When validation fails
    """
    property_instances_list = kim_edn.loads(property_instances)
    if not any(property_instance["property-id"]==energy_property_id for property_instance in property_instances_list):
        raise RuntimeError("Could not find any instances of property %s in the indicated results file" %energy_property_id)
    for property_instance in property_instances_list:
        if property_instance["property-id"]==energy_property_id:
            validate_common_fields(property_instance)
            energy_per_atom_ev = convert_units(
                property_instance["binding-potential-energy-per-atom"]["source-value"],
                property_instance["binding-potential-energy-per-atom"]["source-unit"],
                "eV"
                )
            energy_per_formula_ev = convert_units(
                property_instance["binding-potential-energy-per-formula"]["source-value"],
                property_instance["binding-potential-energy-per-formula"]["source-unit"],
                "eV"
                )            
            num_atoms_per_formula = sum(get_stoich_reduced_list_from_prototype(property_instance["prototype-label"]["source-value"]))
            if abs(num_atoms_per_formula*energy_per_atom_ev-energy_per_formula_ev)>3*EV_ENERGY_EPSILON:
                raise RuntimeError(
                    "Energy per atom %f eV multiplied by %d atoms per formula differs from energy per formula %s eV by more than %f eV" %
                    (energy_per_atom_ev,num_atoms_per_formula,energy_per_formula_ev,EV_ENERGY_EPSILON)
                    )
        
def validate_crystal_structure_npt(property_instances: str, structure_property_id: str):
    """
    Validate "crystal-structure-npt" property

    Args:
        property_instances: Serialized kim-edn property instances object
        structure_property_id: name of the structure property
    Raises:
        RuntimeError: When validation fails        
    """
    property_instances_list = kim_edn.loads(property_instances)
    if not any(property_instance["property-id"]==structure_property_id for property_instance in property_instances_list):
        raise RuntimeError("Could not find any instances of property %s in the indicated results file" %structure_property_id)
    for property_instance in property_instances_list:
        if property_instance["property-id"]==structure_property_id:
            validate_common_fields(property_instance)
            # GNU Units will crash if the absolute temperature is negative
            convert_units(property_instance["temperature"]["source-value"],property_instance["temperature"]["source-unit"],"K")
            stress={
                "xx":property_instance["cell-cauchy-stress"]["source-value"][0],
                "yy":property_instance["cell-cauchy-stress"]["source-value"][1],
                "zz":property_instance["cell-cauchy-stress"]["source-value"][2],
                "yz":property_instance["cell-cauchy-stress"]["source-value"][3],
                "xz":property_instance["cell-cauchy-stress"]["source-value"][4],
                "xy":property_instance["cell-cauchy-stress"]["source-value"][5]
            }
            if "parameter-names" not in property_instance:
                # Cubic crystal, stress must be isotropic
                if not((stress["xx"]==stress["yy"]==stress["zz"]) and (stress["yz"]==stress["xz"]==stress["xy"]==0.)):
                    raise RuntimeError ("Non-isotropic stress state "+str(stress)+" in cubic crystal")
            else:
                if "b/a" not in property_instance["parameter-names"]["source-value"]:
                    if stress["xx"] != stress["yy"]:
                        raise RuntimeError(
                            "Stress component xx %f is not equal to stress component yy %f. This is incompatible with a crystal prototype lacking parameter b/a"%
                            (stress["xx"], stress["yy"])
                            )
                if "c/a" not in property_instance["parameter-names"]["source-value"]:
                    if stress["xx"] != stress["zz"]:
                        raise RuntimeError(
                            "Stress component xx %f is not equal to stress component zz %f. This is incompatible with a crystal prototype lacking parameter c/a"%
                            (stress["xx"], stress["zz"])
                            )
                if "alpha" not in property_instance["parameter-names"]["source-value"]:
                    if stress["yz"] != 0.:
                        raise RuntimeError(
                            "Stress component yz %f is nonzero. This is incompatible with a crystal prototype lacking parameter alpha"%
                            stress["yz"]
                            )
                if "beta" not in property_instance["parameter-names"]["source-value"]:
                    if stress["xz"] != 0.:
                        raise RuntimeError(
                            "Stress component xz %f is nonzero. This is incompatible with a crystal prototype lacking parameter beta"%
                            stress["xz"]
                            )
                if "gamma" not in property_instance["parameter-names"]["source-value"]:
                    if stress["xy"] != 0.:
                        raise RuntimeError(
                            "Stress component xy %f is nonzero. This is incompatible with a crystal prototype lacking parameter gamma"%
                            stress["xy"]
                            )                                                        

    
        


"""
Functions for writing and validating ``binding-energy-crystal`` and ``crystal-structure-npt`` properties
"""

from kim_property import kim_property_create, kim_property_modify
from typing import Dict, List, Union
import numpy as np
import kim_edn
from ..aflow_util import get_stoich_reduced_list_from_prototype, AFLOW, read_shortnames
import os
from kim_property import kim_property_dump
from curses.ascii import isdigit
import subprocess

ENERGY_PROPERTY_ID = "tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal"
STRUCTURE_PROPERTY_ID = "tag:staff@noreply.openkim.org,2023-02-21:property/crystal-structure-npt"

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
                    
def validate_binding_energy_crystal(property_instances: str, energy_property_id: str = ENERGY_PROPERTY_ID):
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
        
def validate_crystal_structure_npt(property_instances: str, structure_property_id: str = STRUCTURE_PROPERTY_ID):
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
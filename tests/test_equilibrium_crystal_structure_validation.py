import os
from crystal_genome_util.property_util.equilibrium_crystal_structure import validate_binding_energy_crystal,validate_crystal_structure_npt

BINDING_ENERGY_DATA_DIR="data/binding_energy_crystal"
CRYSTAL_STRUCTURE_DATA_DIR="data/crystal_structure_npt"


def validation_test(datadir,validate_func):
    for filename in os.listdir(datadir):
        if filename.split(".")[0]=="correct":
            correct=True
        else:
            correct = False
        filepath = os.path.join(datadir,filename)
        with open(filepath) as f:
            prop_inst = f.read()
            if correct:
                validate_func(prop_inst)
                print ("\nFile %s correctly passed validation\n"%filepath)
            else:
                try:
                    validate_func(prop_inst)
                    assert False, "\nFile %s should have failed validation but it passed\n"%filepath
                except RuntimeError as e:
                    print ("File %s correctly failed validation with the following error:"%filepath)
                    print (e)
                    pass

def test_main():
    validation_test(BINDING_ENERGY_DATA_DIR,validate_binding_energy_crystal)
    validation_test(CRYSTAL_STRUCTURE_DATA_DIR,validate_crystal_structure_npt)

if __name__ == "__main__":
    test_main()


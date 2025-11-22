######################################################
# Created: August, 2024
# By: Md Kamruzzaman (@mkamruz)
# Objective:
#          1. Helper methods
######################################################
import os
import logging
from padelpy import padeldescriptor
from collections import OrderedDict
# from datetime import datetime
# import random
from csv import DictReader
import warnings
from time import sleep
from os import remove
# from utils import helper
import tempfile

# Configure logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

PKGDIR = os.getcwd()

class Descriptor:
    def __init__(self):
        self.temp_dir = 'tmp'
        # Create a unique temporary file path that the OS guarantees is writable.
        # delete=False means we handle the deletion explicitly.
        # mode='w' ensures we can write text content.
        with tempfile.NamedTemporaryFile(suffix=".smi", delete=False, mode='w') as temp_smi_file:
            self.smi_file = temp_smi_file.name
            
        #print(f"Temporary SMILES file created at: {self.smi_file}")

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode='w') as temp_smi_file:
            self.output_csv = temp_smi_file.name
            
        #print(f"Temporary CSV file created at: {self.output_csv}")



    # This method is copied from https://github.com/ecrl/padelpy/
    # We have modified a few lines to suit our needs.
    def from_smiles(self, smiles) -> OrderedDict:
        """ from_smiles: converts SMILES string to QSPR descriptors/fingerprints.

        Args:
            smiles (str, list): SMILES string for a given molecule, or a list of
                SMILES strings

        Returns:
            list or OrderedDict: if multiple SMILES strings provided, returns a
                list of OrderedDicts, else single OrderedDict; each OrderedDict
                contains labels and values for each descriptor generated for each
                supplied molecule
        """
        #output_csv = None
        descriptors =True # if `True`, calculates descriptors
        fingerprints = False # if `True`, calculates fingerprints
        timeout = 10 # maximum time, in seconds, for conversion, default=60s
        maxruntime = -1 # maximum running time per molecule in seconds. default=-1.
        threads = -1 # number of threads to use; defaults to -1 for max available
        save_csv = False
        max_run = 1 # Number of trials
                
        # unit conversion for maximum running time per molecule
        # seconds -> milliseconds
        if maxruntime != -1:
            maxruntime = maxruntime * 1000

        #timestamp = datetime.now().strftime("%Y%m%d%H%M%S%f")#[:-3]
        #filename = timestamp + str(random.randint(int(1e8),int(1e9)))
        #smi_file_path = helper.resource_path(os.path.join(self.temp_dir, "{}.smi".format(filename)))

        with open(self.smi_file, "w") as smi_file:
            if type(smiles) == str:
                smi_file.write(smiles)
            elif type(smiles) == list:
                smi_file.write("\n".join(smiles))
            else:
                raise RuntimeError("Unknown input format for `smiles`: {}".format(
                    type(smiles)
                ))
        smi_file.close()

        # if output_csv is None:
        #     save_csv = False
        #     output_csv = helper.resource_path(os.path.join(self.temp_dir, "{}.csv".format(timestamp)))


        for attempt in range(max_run):
            try:
                padeldescriptor(
                    mol_dir=self.smi_file,
                    d_file=self.output_csv,
                    convert3d=True,
                    retain3d=True,
                    d_2d=descriptors,
                    d_3d=descriptors,
                    fingerprints=fingerprints,
                    sp_timeout=timeout,
                    retainorder=True,
                    maxruntime=maxruntime,
                    threads=threads
                )
                break
            except RuntimeError as exception:
                if attempt == max_run-1:
                    remove(self.smi_file)
                    if not save_csv:
                        sleep(0.5)
                        try:
                            remove(self.output_csv)
                        except FileNotFoundError as e:
                            warnings.warn(e, RuntimeWarning)
                    raise RuntimeError(exception)
                else:
                    continue
            except KeyboardInterrupt as kb_exception:
                remove(self.smi_file)
                if not save_csv:
                    try:
                        remove(self.output_csv)
                    except FileNotFoundError as e:
                        warnings.warn(e, RuntimeWarning)
                raise kb_exception

        with open(self.output_csv, "r", encoding="utf-8") as desc_file:
            reader = DictReader(desc_file)
            rows = [row for row in reader]
        desc_file.close()

        remove(self.smi_file)
        if not save_csv:
            remove(self.output_csv)

        if type(smiles) == list and len(rows) != len(smiles):
            raise RuntimeError("PaDEL-Descriptor failed on one or more mols." +
                            " Ensure the input structures are correct.")
        elif type(smiles) == str and len(rows) == 0:
            raise RuntimeError(
                "PaDEL-Descriptor failed on {}.".format(smiles) +
                " Ensure input structure is correct."
            )

        for idx, r in enumerate(rows):
            if len(r) == 0:
                raise RuntimeError(
                    "PaDEL-Descriptor failed on {}.".format(smiles[idx]) +
                    " Ensure input structure is correct."
                )

        for idx in range(len(rows)):
            del rows[idx]["Name"]

        if type(smiles) == str:
            return rows[0]
        
        return rows
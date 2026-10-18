"""
Pseudo dinámica molecular
A partir de un código PDB genera un gif con el pseudo movimiento de la proteína plegada

Autor: Adrián Berenguer Agustí
version 5
"""

import argparse
import os
import subprocess
import sys
from typing import TextIO

from Bio import Blast, Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import sys

#####################################
#                                   #
#             FUNCTIONS             #
#                                   #
#####################################


def open_file(path: str, mode: str = "r"):
    """Opens file"""
    try:
        file = open(path, mode)
    except Exception as e:
        print(f"Opening file {path} failed: {e}")
        sys.exit(1)

    return file


def close_file(file: TextIO):
    """Closes file"""
    try:
        file.close()
    except Exception as e:
        print(f"Closing file failed: {e}")
        sys.exit(1)


def get_sequence(id_pdb: str) -> SeqIO:
    """Using an id_pdb gets the sequence of the protein"""

    Entrez.email = "adrian9420agusti@gmail.com"
    query = id_pdb + "[All Fields] AND pdb[filter]"

    handle = Entrez.esearch(db="protein", term=query)
    records = Entrez.read(handle)
    id = records["IdList"][0]
    handle.close()

    fasta = Entrez.efetch(db="protein", id=id, rettype="fasta")
    fasta_record = SeqIO.read(fasta, "fasta")

    return fasta_record


def get_pdb_variants(fasta: SeqIO, n: int = 5) -> str:
    """Using a sequence it does a BLAST"""

    Blast.email = "adrian9420agusti@gmail.com"
    result_blast = Blast.qblast(
        "blastp", "pdb", fasta.seq, expect=0.05, perc_ident=99, hitlist_size=n
    )
    data = result_blast.read()
    result_blast.close()

    return data


def save_blast(data, output_dir: str):
    """creates an XML with the blast results"""

    path = os.path.join(output_dir, "alignment.txt")
    file = open_file(path, "wb")
    file.write(data)
    close_file(file)

    return path


def list_pdb(path: str) -> dict:
    """With one input file extracts pdb id and sequence"""

    file = open_file(path)
    blast = NCBIXML.read(file)
    close_file(file)
    pdb_id_dict = {}

    for alignment in blast.alignments:
        pdb_id = alignment.accession[:4].upper()
        if pdb_id not in pdb_id_dict:
            pdb_id_dict[pdb_id] = alignment.hsps[0].sbjct

    return pdb_id_dict


def read_pdb_list(path: str) -> dict:
    """Reads a plain text file containing one PDB ID per line"""
    file = open_file(path)
    pdb_id_dict = {}

    for line in file:
        line_clean = line.strip().upper()

        if line_clean and not line_clean.startswith("#"):

            pdb_id = line_clean[:4]
            if len(pdb_id) == 4 and pdb_id.isalnum():
                # Not a PDB
                pdb_id_dict[pdb_id] = ""

    close_file(file)

    if not pdb_id_dict:
        print(f"id PDBs not found at {path}")
        sys.exit(1)

    return pdb_id_dict


def already_down(pdb_id_dict, path) -> dict:
    """Verifies which pdb have been already downloaded
    and returns not downloaded pdb ids"""
    already = os.listdir(path)

    if already:
        pdb_files = [
            file.upper()
            for file in already
            if file.upper().endswith((".PDB", ".ENT", ".CIF"))
        ]

        new_pdb_id_dict = {
            pdb_id: seq
            for pdb_id, seq in pdb_id_dict.items()
            if not any(pdb_id in filename for filename in pdb_files)
        }

    else:
        new_pdb_id_dict = pdb_id_dict

    return new_pdb_id_dict


def download_3D(pdb_id_dict: dict, output_dir: str) -> str:
    """Obtains 3D structure files from the dict with pdb_ids"""

    pdbs_dir = os.path.join(output_dir, "pdbs")
    os.makedirs(pdbs_dir, exist_ok=True)

    pdb_id_dict = already_down(pdb_id_dict, pdbs_dir)

    pdb_list = PDBList()
    for pdb_id in pdb_id_dict.keys():
        pdb_list.retrieve_pdb_file(pdb_id, pdir=pdbs_dir, file_format="mmCif")

    return pdbs_dir


def _count_atoms(filepath: str) -> int:
    """Atom counter"""
    if filepath.lower().endswith(".cif"):
        try:
            mmcif_dict = MMCIF2Dict(filepath)
            return len(mmcif_dict.get("_atom_site.id", []))
        except Exception:
            n_atom = 0
    else:
        with open(filepath, "r") as f:
            n_atom = sum(1 for line in f if line.startswith(("ATOM", "HETATM")))

    return n_atom


def resolve_chains_to_auth(file_path: str, requested_chains: str) -> str:
    """
    Author - original author chain translator. Most of PDB chains are wrong annotated
    """
    if file_path.lower().endswith(".cif"):

        try:
            mmcif_dict = MMCIF2Dict(file_path)
            labels = mmcif_dict.get("_atom_site.label_asym_id", [])
            auths = mmcif_dict.get("_atom_site.auth_asym_id", [])

            if labels and auths:   

                unique_pairs = set(zip(labels, auths))
                auth_set = set(auths)
                label_to_auth = {label: author for label, author in unique_pairs}

                sub_chains = requested_chains.replace(",", "+").split("+")
                resolved_list = []

                for chain in sub_chains:
                    chain = chain.strip()
                    if chain in auth_set:
                        resolved_list.append(chain)
                    elif chain in label_to_auth:
                        resolved_list.append(label_to_auth[chain])
                    else:
                        resolved_list.append(chain)

                requested_chains = "+".join(sorted(set(resolved_list)))

        except:
            requested_chains = requested_chains
            
        
    return requested_chains


#####################################
#                                   #
#   Pymol complementary functions   #
#                                   #
#####################################


def _filter_pdbs(files_path: str, output_dir: str) -> tuple[list[str], int]:
    """PDBs filtered respect to number of atoms compared to original sequence.
    Creates a file for knowing whcih PDB have been selected"""

    pdb_files = [
        file
        for file in os.listdir(files_path)
        if file.lower().endswith((".pdb", ".ent", ".cif"))
    ]
    # Order PDB of txt
    list_arg = next((sys.argv[i+1] for i, arg in enumerate(sys.argv) if arg in ("-l", "--list") and i+1 < len(sys.argv)), None)

    if list_arg and os.path.exists(list_arg):
        f = open_file(list_arg, "r")
        orden_txt = [line.strip()[:4].lower() for line in f if line.strip() and not line.strip().startswith('#')]
        close_file(f)
        pdb_files.sort(key=lambda x: orden_txt.index(x[:4].lower()) if x[:4].lower() in orden_txt else 9999)
        
    pdb_files_clean = []

    if pdb_files:

        max_size = _count_atoms(os.path.join(files_path, pdb_files[0]))
        
        for file in pdb_files:
            atom_count = _count_atoms(os.path.join(files_path, file))
            # Range limit for filtering PDB with excess of too low nº of atoms
            # YOU CAN CHANGE THIS RANGE. default 0.5 and 2. For accepting al 0 and 2000
            if (max_size * 0.5) < atom_count < (max_size * 2):
                pdb_files_clean.append(file)
        # Creates a file for knowing whcih PDB have been selected
        filtered_file = open_file(os.path.join(output_dir, "filtered_pdb.txt"), "w")
        
        for pdb in pdb_files_clean:
            filtered_file.write(pdb + "\n")

        close_file(filtered_file)

    return pdb_files_clean, len(pdb_files)


def _write_load_and_extract(pymol_file: TextIO, files_path: str, pdb_list: list,
    chains: str, only_backbone: bool = False):
    """Says pymol which pdb must use and which chain(s)
    Only backbone used for discriminating between superposition (looks full protein)
    and structural alignment (look only alfa carbons)"""

    pymol_file.write("set cif_use_auth, 1\n")

    for idx, file in enumerate(pdb_list, start=1):
        # Absolute paths
        abs_pdb_path = os.path.abspath(os.path.join(files_path, file))
        resolved_chains = resolve_chains_to_auth(abs_pdb_path, chains)
        # Says Pymol which PDB has to load
        pymol_file.write(f"load {abs_pdb_path}, raw_{idx}\n")

        if only_backbone:
            # For Structural alignment
            pymol_file.write(
                f"create c_{idx}, raw_{idx} and polymer and chain {resolved_chains} and name CA+C+N+O and (not alt B)\n"
            )
        else:
            # For superposition
            pymol_file.write(
                f"create c_{idx}, raw_{idx} and polymer and chain {resolved_chains} and (not alt B)\n"
            )

        pymol_file.write(f"delete raw_{idx}\n")


def _write_trajectory_build(pymol_file: TextIO, num_structures: int):
    """Says pymol trajectory of each pdb"""

    for idx in range(1, num_structures + 1):
        pymol_file.write(f"create trajectory, c_{idx}, 1, {idx}\n")
    pymol_file.write("delete c_*\n")


def _write_render_and_export(
    pymol_file: TextIO, obj_name: str, mset_cmd: str, output_dir: str, color_mode: str
):
    """Says pymol visualizing configuration"""
    # Same orientation for all PDBs
    pymol_file.write(f"orient {obj_name}\n")
    # Smoothing. Troubleshooting
    #pymol_file.write(f"smooth {obj_name}, 5, 3\n")
    pymol_file.write("hide all\n")
    pymol_file.write(f"show cartoon, {obj_name}\n")
    # Color options
    if color_mode == "1":
        pymol_file.write(f"spectrum count, green cyan deepblue, {obj_name}\n")
    elif color_mode == "2":
        pymol_file.write(f"spectrum resi, green cyan deepblue, {obj_name}\n")
    else:
        pymol_file.write(f"color {color_mode}, {obj_name}\n")
    # RAM liberation
    pymol_file.write("set cache_frames, 0\n")
    pymol_file.write("set ray_trace_frames, 0\n")
    pymol_file.write(f"mset {mset_cmd}\n")
    pymol_file.write("set movie_loop, 1\n")

    gif_path = os.path.abspath(os.path.join(output_dir, "animated.gif"))
    pymol_file.write(f"movie.produce {gif_path}\n")


#####################################
#                                   #
#          Pymol functions          #
#                                   #
#####################################


def pymol_script_superposition(files_path: str, output_dir: str, 
                               chains: str, color_mode: str) -> str:
    """Superposition approach"""

    pdb_clean, n = _filter_pdbs(files_path, output_dir)
    print(f"\n{len(pdb_clean)} filtered PDBS of ", n)

    pml_path = os.path.join(output_dir, "molecular_dynamics_superposition.pml")
    pymol_file = open_file(pml_path, "w")

    _write_load_and_extract(
        pymol_file, files_path, pdb_clean, chains, only_backbone=False
    )
    # Superposition
    for idx in range(2, len(pdb_clean) + 1):
        pymol_file.write(f"super c_{idx}, c_1\n")

    _write_trajectory_build(pymol_file, len(pdb_clean))

    mset_str = " ".join(f"{i}x10" for i in range(1, len(pdb_clean) + 1))
    _write_render_and_export(pymol_file, "trajectory", mset_str, output_dir, color_mode)

    close_file(pymol_file)

    return pml_path


def pymol_script_morph(files_path: str, output_dir: str, chains: str, 
                       color_mode: str = "2") -> str:

    """Structural alignment approach"""
    pdb_clean, n = _filter_pdbs(files_path, output_dir)
    print(f"\n{len(pdb_clean)} PDBs filtered of", n)

    pml_path = os.path.join(output_dir, "molecular_dynamics_morph.pml")
    pymol_file = open_file(pml_path, "w")

    # PDB loading to pymol
    _write_load_and_extract(
        pymol_file, files_path, pdb_clean, chains, only_backbone=True
    )

    # Structural alingment
    for idx in range(2, len(pdb_clean) + 1):
        pymol_file.write(f"align c_{idx}, c_1\n")

    # trajectories
    for idx in range(1, len(pdb_clean) + 1):
        pymol_file.write(f"create trajectory, c_{idx}, 1, {idx}\n")
    pymol_file.write("delete c_*\n")

    # Interporlation with morph and linear + super position calculus
    n_structs = len(pdb_clean)
    #pymol_file.write(f"morph pseudo_md, trajectory, trajectory, 1, 2, method=linear, match=super\n")
    pymol_file.write(f"morph pseudo_md, trajectory, trajectory, 1, 2, method=linear, match=align\n")

    for i in range(2, n_structs):
        # pymol_file.write(
        #     f"morph tramo_{i}, trajectory, trajectory, {i}, {i+1}, method=linear, match=super\n"
        # )
        pymol_file.write(
            f"morph tramo_{i}, trajectory, trajectory, {i}, {i+1}, method=linear, match=align\n"
        )
        pymol_file.write(f"create pseudo_md, tramo_{i}, 0, -1\n")
        pymol_file.write(f"delete tramo_{i}\n")

    pymol_file.write("delete trajectory\n")
    total_frames = (n_structs - 1) * 30

    # Gif generation
    _write_render_and_export(pymol_file, "pseudo_md", f"1 -{total_frames}", output_dir, color_mode)

    close_file(pymol_file)

    return pml_path


#####################################
#                                   #
#               MAIN                #
#                                   #
#####################################


def main():
    print(
        """
#####################################
#                                   #
#         DYNAMIC SIMULATOR         #
#                                   #
#####################################
"""
    )
    parser = argparse.ArgumentParser(description="GIF molecular dynamics simulator")
    parser.add_argument(
        "-id",
        "--id_pdb",
        type=str,
        help="id PDB",
    )
    parser.add_argument(
        "-al", "--alignment", type=str, help="BLAST directory as optional input. Accurate results not guaranted" \
        "Accurate results guaranted if percentage identity is 99 and number of hits is 5"
    )
    parser.add_argument(
        "-l",
        "--list",
        type=str,
        help="For selecting only some PDBs. The file must contain one PDB per line",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default=".",
        help="Output directory, optional",
    )
    parser.add_argument(
        "-ht", "--hits", type=int, default=5, help="BLAST hits, default 5"
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=["superposition", "structural_alignment"],
        default="superposition",
        help="Algorithmic approach: superposition (fails less) "
        "or structural alignment (more accurate) with interpolation. "
        "Default: superposition",
    )
    parser.add_argument(
        "-c",
        "--chains",
        type=str,
        default="A",
        help="Which chains you want to send to PyMOL (eg: Am A+B, A,B). Default A. " \
        "Take care, most PDB chains are wrong annotated",
    )
    parser.add_argument(
        "-col",
        "--color",
        type=str,
        default="2",
        help="Color: '1' atom model, '2' residue model, or choose (eg: 'deepblue'). Default: 2",
    )

    args = parser.parse_args()

    # 
    if not (args.id_pdb or args.alignment or args.list):
        parser.error(
            "Lack of arguments. Add at least: -id, -al or o -l"
        )

    out_dir = os.path.abspath(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)
    print("\nOutput directory: ", out_dir)

    pymol_chains_syntax = args.chains.replace(",", "+").replace(" ", "+")
    print(f"Requested Target Chains: {pymol_chains_syntax}")

    # PDB list option
    if args.list:
        list_file = os.path.abspath(args.list)
        if not os.path.exists(list_file):
            print("\nPDB list not found at: ", list_file)
            sys.exit(1)
        print("\nLoading PDBs from: ", list_file)
        pdb_id_dict = read_pdb_list(list_file)

    # Blast added as argument
    elif args.alignment:
        alignment_file = os.path.abspath(args.alignment)
        if not os.path.exists(alignment_file):
            print("\nBLAST not found at: ", alignment_file)
            sys.exit(1)
        print("\nLoading BLAST from: ", alignment_file)
        pdb_id_dict = list_pdb(alignment_file)

    # Do BLAST
    else:
        fasta = get_sequence(args.id_pdb)
        print("\nExecuting BLAST...")
        data = get_pdb_variants(fasta, args.hits)

        alignment_file = save_blast(data, out_dir)
        print("\nBLAST results saved to: ", alignment_file)
        pdb_id_dict = list_pdb(alignment_file)

    files_path = download_3D(pdb_id_dict, out_dir)
    print("\n3D structures downloaded/verified (mmCIF format).")

    # Structural alignment mode
    if args.method == "structural_alignment":
        pml_script_path = pymol_script_morph(
            files_path, out_dir, pymol_chains_syntax, args.color
        )
        print(
            "\nPYMOL license is required for interpolation execution. "
            f"Please, run the command --pymol {os.path.basename(pml_script_path)}--"
        )
        print("\nPYMOL script generated located at:", pml_script_path)

    # Superposition mode
    else:
        pml_script_path = pymol_script_superposition(
            files_path, out_dir, pymol_chains_syntax, args.color
        )
        print("\nPYMOL script generated.")
        print("\nExecuting PyMOL in background...")

        log_file_path = os.path.join(out_dir, "pymol_process.log")
        log_file = open_file(log_file_path, "w")

        process = subprocess.run(
            ["pymol", "-c", "-q", pml_script_path],
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )

        if process.returncode != 0:
            print(f"\nPyMOL failed. More details at: {log_file_path}")
        else:
            print(f"\nGIF successfully located at: {out_dir}")


if __name__ == "__main__":
    main()
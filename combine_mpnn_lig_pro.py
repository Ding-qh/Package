#!/data/PRG/tools/Biomolecules/miniconda3/bin/python -B
import os
import sys
import argparse
import subprocess
import re
from subprocess import Popen, PIPE

# Chinese punctuation -> English (for normalizing user input in string args)
_CN_TO_EN_PUNCT = str.maketrans({
    '，': ',',
    '。': '.',
    '：': ':',
    '；': ';',
    '（': '(',
    '）': ')',
    '【': '[',
    '】': ']',
    '“': '"',
    '”': '"',
    '‘': "'",
    '’': "'",
})


def norm_arg_punct(s):
    """Convert Chinese punctuation in string to English. Return unchanged if not str or empty."""
    if not s or not isinstance(s, str):
        return s
    return s.translate(_CN_TO_EN_PUNCT)


def norm_file_punct_save(input_path, output_path=None):
    """
    Read file, convert Chinese punctuation to English in content, save as new file.
    If output_path is None, save to same dir with suffix '_en' before extension (e.g. bias.txt -> bias_en.txt).
    Returns path to the new file.
    """
    with open(input_path, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    normalized = [norm_arg_punct(line) for line in lines]
    if output_path is None:
        base, ext = os.path.splitext(input_path)
        output_path = base + '_en' + ext
    with open(output_path, 'w', encoding='utf-8') as f:
        f.writelines(normalized)
    return output_path


def run_ext_cmder(cmds: list, query:str=None) -> str:
    ''' Run external command.
    - cmds: list of [cmder, option, value...]
    '''
    process = Popen(cmds, stdout=sys.stdout, stderr=sys.stderr)
    #process = Popen(cmds, stdout=PIPE, stderr=sys.stderr)
    process.communicate() if query is None else process.communicate(query)
    retcode = process.wait()
    if retcode:
        raise RuntimeError('Run Ext Cmd Failed')


def parse_bias_AA_file_for_ligand_mpnn(bias_file_path):
    """
    Read residue bias from file, convert to LigandMPNN --bias_AA string format.
    File: one line per AA, "letter,value" (e.g. H,1.39). Returns e.g. "H:1.39,D:1.0".
    """
    pairs = []
    with open(bias_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',', 1)]
            if len(parts) != 2:
                continue
            aa, val = parts[0], parts[1]
            if len(aa) != 1:
                continue
            try:
                pairs.append(f"{aa.upper()}:{float(val)}")
            except ValueError:
                continue
    return ",".join(pairs) if pairs else ""


def convert_positions_chain_range_to_default(positions_raw, chains_raw=None):
    """
    Convert positions from "A1-10,A30,B12-25,B30" to default format:
    "1 2 ... 10 30, 12 13 ... 25 30" (grouped by --chains order).
    """
    if not positions_raw:
        return positions_raw

    pattern = re.compile(r"^([A-Za-z])\s*(\d+)(?:\s*-\s*(\d+))?$")
    chain_order = [c.strip() for c in chains_raw.split(",") if c.strip()] if chains_raw else []
    chain_to_positions = {ch: [] for ch in chain_order}
    chain_to_seen = {ch: set() for ch in chain_order}

    entries = [e.strip() for e in positions_raw.split(",") if e.strip()]
    if not entries:
        raise ValueError("Invalid --positions: empty value.")

    for entry in entries:
        m = pattern.match(entry)
        if not m:
            raise ValueError(
                f"Invalid --positions entry '{entry}'. Expected format like A1-10 or B30."
            )
        chain_id, start_str, end_str = m.group(1), m.group(2), m.group(3)
        if chain_id not in chain_to_positions:
            if chains_raw:
                raise ValueError(
                    f"Chain '{chain_id}' in --positions is not in --chains '{chains_raw}'."
                )
            # Infer chain order from positions by first appearance.
            chain_order.append(chain_id)
            chain_to_positions[chain_id] = []
            chain_to_seen[chain_id] = set()
        start = int(start_str)
        end = int(end_str) if end_str is not None else start
        if start <= 0 or end <= 0:
            raise ValueError(f"Invalid residue index in '{entry}': index must be positive.")
        if start > end:
            raise ValueError(f"Invalid range in '{entry}': start cannot be larger than end.")

        for pos in range(start, end + 1):
            if pos not in chain_to_seen[chain_id]:
                chain_to_seen[chain_id].add(pos)
                chain_to_positions[chain_id].append(pos)

    # Keep compatibility with current script behavior: one position group per chain.
    missing = [ch for ch in chain_order if not chain_to_positions[ch]]
    if missing:
        raise ValueError(
            "Missing positions for chain(s): "
            + ",".join(missing)
            + ". Please provide entries like "
            + ",".join([f"{ch}1-2" for ch in missing])
            + "."
        )

    groups = []
    for ch in chain_order:
        groups.append(" ".join(str(p) for p in chain_to_positions[ch]))
    return ",".join(chain_order), ", ".join(groups)

if __name__ == '__main__':
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_path', required=True, type=str, help="The PDB file to be designed.")
    parser.add_argument('--chains', required=False, type=str, default="",  help="Define which chains need to be designed , e.g. 'A,B', use comma to separate. Optional when --positions uses chain-prefix format like A1-10,B20.")

    parser.add_argument('--num_seq_per_target', type=int, default=5, help="Number of sequences to generate per target.")
    parser.add_argument('--sampling_temp', type=float, default=0.1, help="Sampling temperature for amino acids, T=0.0 means taking argmax, T>>1.0 means sample randomly. Suggested values 0.1, 0.15, 0.2, 0.25, 0.3. Higher values will lead to more diversity.")
    
    parser.add_argument('--positions', type=str, default=None,  help="When the option pdb_path_chains='A C', the positions='1 2 3 4 5 6 7 8 23 25, 10 11 12 13 14 15 16 17 18 19 20 40' means fixing or designing residues 1 2 3...25 in chain A and residues 10 11 12...40 in chain C. Note that The first amino acid in the chain corresponds to 1 and not PDB residues index for now.")
    parser.add_argument('--pos_type', type=str, default=None,  choices=["Fix","Design"], help="To Fix or Design the residues in positions.")
    parser.add_argument('--omit_AAs', type=str, default="X",  help="Specify which amino acids should be omitted in the generated sequence, e.g. 'AC' would omit alanine and cystine.")
    
    # proteinMPNN
    parser.add_argument('--homomer', action="store_true", default=False, help="Design homomer")
    #parser.add_argument('--unconditional_probs_only', type=int, default=0, help="0 for False, 1 for True; output unconditional probabilities p(s_i given backbone) in one forward pass")
    parser.add_argument("--save_probs", type=str, default="False", help="0 for False, 1 for True; save MPNN predicted probabilites per position")
    parser.add_argument('--soluble', action="store_true", default=False, help="Flag to load ProteinMPNN weights trained on soluble proteins only.")
    parser.add_argument('--antibody', action="store_true", default=False, help="Flag to design antibody.")

    # ligandMPNN
    parser.add_argument('--ligandMPNN', action="store_true", default=False, help="Flag to design ligandMPNN.")
    parser.add_argument('--pack_side_chains', action="store_true", default=False, help="to run side chain packer or not")

    # Amino acid bias file (one line per AA: residue,value e.g. H,1.39; ProteinMPNN gets path, LigandMPNN gets AA:val string)
    parser.add_argument('--bias_AA_file', type=str, default='', help="Path to file with per-residue bias (one line per AA: single letter, comma, value, e.g. 'H,1.39'). Passed as file to ProteinMPNN, as --bias_AA string to LigandMPNN.")

    args = parser.parse_args()
    # Normalize string args: Chinese punctuation -> English
    args.chains = norm_arg_punct(args.chains) if args.chains else args.chains
    args.positions = norm_arg_punct(args.positions) if args.positions else args.positions
    args.omit_AAs = norm_arg_punct(args.omit_AAs) if args.omit_AAs else args.omit_AAs
    args.bias_AA_file = norm_arg_punct(args.bias_AA_file) if args.bias_AA_file else args.bias_AA_file
    if getattr(args, 'save_probs', None):
        args.save_probs = norm_arg_punct(args.save_probs)
    if getattr(args, 'pos_type', None):
        args.pos_type = norm_arg_punct(args.pos_type)
    args.pdb_path = norm_arg_punct(args.pdb_path)
    if args.positions:
        args.chains, args.positions = convert_positions_chain_range_to_default(args.positions, args.chains)
    elif not args.chains:
        raise ValueError("Either provide --positions (with chain-prefix format) or provide --chains.")

    #print(args.chains, args.positions)
    # Convert input file content (e.g. bias_AA_file): Chinese punct -> English, save as new file and use it
    if args.bias_AA_file and os.path.isfile(args.bias_AA_file):
        args.bias_AA_file = norm_file_punct_save(args.bias_AA_file)

    cmd = []
    if args.ligandMPNN:        
        cmd.append("/data/PRG/tools/Biomolecules/apps/LigandMPNN/my_run.py")
        
        if args.pack_side_chains:
            cmd.append("--pack_side_chains")
        if args.bias_AA_file:
            bias_str = parse_bias_AA_file_for_ligand_mpnn(args.bias_AA_file)
            if bias_str:
                cmd.append("--bias_AA")
                cmd.append(bias_str)
    else:
        cmd.append("/data/PRG/tools/Biomolecules/apps/ProteinMPNN/run.py")
        if args.homomer:
            cmd.append("--homomer")
        if args.soluble:
            cmd.append("--soluble")
        if args.antibody:
            cmd.append("--antibody")
        if args.save_probs:
            cmd.append("--save_probs")
            cmd.append(args.save_probs)
        if args.bias_AA_file:
            cmd.append("--bias_AA_file")
            cmd.append(args.bias_AA_file)

    cmd.append("--pdb_path")
    cmd.append(args.pdb_path)
    cmd.append("--chains")
    cmd.append(args.chains)
    cmd.append("--num_seq_per_target")
    cmd.append(str(args.num_seq_per_target))
    cmd.append("--sampling_temp")
    cmd.append(str(args.sampling_temp))
    
    if args.positions:
        cmd.append("--positions")
        cmd.append(args.positions)
    if args.pos_type:
        cmd.append("--pos_type")
        cmd.append(args.pos_type)
    if args.omit_AAs:
        cmd.append("--omit_AAs")
        cmd.append(args.omit_AAs)

    print(cmd)
    run_ext_cmder(cmd)
import argparse
import ast
import pathlib

def create_parser():
    parser = argparse.ArgumentParser(
        description="Simulate Sequences by inferred parameters"
    )

    parser.add_argument(
        "--abspath",
        type=pathlib.Path,
        help="The absolute path of this project"
    )

    parser.add_argument(
        "--dataset_location",
        type=pathlib.Path,
        help="The path which contain all the data"
    )

    parser.add_argument(
        "--rawAAdir",
        type=pathlib.Path,
        help="Specify the raw amino acid sequences' location"
    )

    parser.add_argument(
        "--rawCDSdir",
        type=pathlib.Path,
        help="Specify the raw Codon sequences' location"
    )

    parser.add_argument(
        "--dat_dir",
        type=str,
        help="Substituition Matrix Dat file located"
    )

    parser.add_argument(
        "--dat_file",
        type=str,
        help="Substitution Matrix Dat file name"
    )

    parser.add_argument(
        "--combine_dir",
        type=str,
        help="Combine the sequences into multiple huge file"
    )

    parser.add_argument(
        "--freq_mode",
        type=str,
        help="Specify the frequency mode when generate sequences"
    )

    parser.add_argument(
        "--model_location",
        type=str,
        help="PyTorch model file OR name of pretrained model to download (see README for models)",
    )

    parser.add_argument(
        "--fasta_file",
        type=pathlib.Path,
        help="FASTA file on which to extract representations",
    )

    parser.add_argument(
        "--repr_layers",
        type=int,
        default=[12],
        nargs="+",
        help="layers indices from which to extract representations (0 to num_layers, inclusive)",
    )

    parser.add_argument(
        "--bottleneck_weight",
        type=pathlib.Path,
        help="Extract the Bottleneck weight From Path"
    )

    parser.add_argument(
        "--include",
        type=str,
        nargs="+",
        choices=["mean", "per_tok", "bottleneck"],
        help="specify which representations to return"
    )

    parser.add_argument(
        "--device",
        type=str,
        help="specify the device to load model and data"
    )

    parser.add_argument(
        "--case_sps1",
        type=str,
        nargs="+",
        help="specify the species which are involved in convergence detection group1"
    )

    parser.add_argument(
        "--case_sps2",
        type=str,
        nargs="+",
        help="specify the species which are involved in convergence detection group2"
    )

    parser.add_argument(
        "--embed_single",
        type=ast.literal_eval,
        choices=[True, False],
        help="If this parameter is specified as True, program will only embed the case species instead of all mammals"
    )

    parser.add_argument(
        "--simtime",
        type=int,
        default=100,
        help="specify the null simulation times"
    )

    parser.add_argument(
        "--seqid",
        type=str,
        help="sequence id"
    )

    parser.add_argument(
        "--codeml_str",
        type=pathlib.Path,
        help="specify the path of codeml control file template"
    )

    parser.add_argument(
        "--tree_file",
        type=pathlib.Path,
        help="specify the path of speices tree topology file"
    )

    return parser
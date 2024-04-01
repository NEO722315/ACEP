import os
from ete3 import Tree
import sys
from Bio.Phylo.PAML import codeml
sys.path.append("../../")
import args_prepare

def parse_prepare():
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    return args

class Fileprepare:

    def __init__(self, args):
        self.args = args
        self.dataset_location = self.args.dataset_location
        self.codeml_str = self.args.codeml_str
        self.tree_file = self.args.tree_file
        self.dat_file = self.args.dat_file

    def _judge_filename(self, seqid):
        # make sure the cleaned fasta files exist
        assert os.path.exists(f"{self.dataset_location}/{seqid}/{seqid}.aa.cleaned.fasta")

    def _exract_species(self, seqid):
        species_set = [species[1:-1] for species in open(f"{self.dataset_location}/{seqid}/{seqid}.aa.cleaned.fasta", "r").readlines()[::2]]
        return species_set

    def _make_pamlparas(self, seqid):
        os.makedirs(f"{self.dataset_location}/{seqid}/{seqid}_pamlparas", exist_ok=True)
        os.system(f"cp -f {self.dataset_location}/{seqid}/{seqid}.aa.cleaned.fasta {self.dataset_location}/{seqid}/{seqid}_pamlparas")
        codeml_str = open(self.codeml_str, "r").read()
        codeml_str = codeml_str.replace("SEQ", f"./{seqid}.aa.cleaned.fasta")
        codeml_str = codeml_str.replace("NWK", f"./{seqid}.nwk")
        os.system(f"cp -f {self.dat_file} {self.dataset_location}/{seqid}/{seqid}_pamlparas")
        file = open(f"{self.dataset_location}/{seqid}/{seqid}_pamlparas/{seqid}.ctl", "w")
        print(codeml_str, file=file)
        file.close()

    def _make_tree(self, seqid):
        tree_nwk = open(self.args.tree_file, "r").read()
        tree = Tree(tree_nwk, format=9)
        species_set = self._exract_species(seqid)
        tree.prune(species_set)
        tree.write(outfile=f"{self.dataset_location}/{seqid}/{seqid}_pamlparas/{seqid}.nwk", format=9)

    def pamlinfer(self, seqid):
        working_dir = f"{self.args.abspath}/{self.dataset_location}/{seqid}/{seqid}_pamlparas"
        os.chdir(working_dir)
        os.system(f"codeml ./{seqid}.ctl")

    def main(self):
        for seqid in os.listdir(self.dataset_location):
            self._judge_filename(seqid)
            self._make_pamlparas(seqid)
            self._make_tree(seqid)
            self.pamlinfer(seqid)

if __name__ == '__main__':
    args = parse_prepare()
    fp = Fileprepare(args)
    fp.main()

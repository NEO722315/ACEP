from ete3 import Tree
import numpy
import sys
import os

sys.path.append("./Simulations/scripts")
import sequence

sys.path.append("../../")
import args_prepare

class CalcExpConv:

    def __init__(self, args):
        self.args = args
        self.freq_mode = self.args.freq_mode
        self.aas = 'ARNDCQEGHILKMFPSTWYV'

    def main(self, SeqID):
        os.makedirs(f"{self.args.dataset_location}/{SeqID}", exist_ok=True)
        fasta_file_path = f"{self.args.dataset_location}/{SeqID}/{SeqID}.aa.cleaned.fasta"
        content = open(fasta_file_path, 'r').readlines()
        ancestral_seq = self.get_ancestral_sequence(f"{self.args.dataset_location}/{SeqID}/{SeqID}_pamlparas/rst")
        rates = self.read_rate(f"{self.args.dataset_location}/{SeqID}/{SeqID}_pamlparas/rates")
        smat = self.read_mat(f"{self.args.abspath}/{self.args.dat_file}")
        if not os.path.exists(f"{self.args.dataset_location}/{SeqID}/{SeqID}_inferred_tree/tree_inferred_branch.nwk"):
            self.get_inferred_tree(f"{self.args.dataset_location}/{SeqID}/{SeqID}_pamlparas", f"{self.args.dataset_location}/{SeqID}/{SeqID}_inferred_tree")
        inferred_tree_path = f"{self.args.dataset_location}/{SeqID}/{SeqID}_inferred_tree/tree_inferred_branch.nwk"
        seqs = [content[i][:-1] for i in range(1, len(content), 2)]
        freq = self.get_freqs(seqs, self.freq_mode)
        return ancestral_seq, rates, smat, freq, inferred_tree_path

    def get_ancestral_sequence(self, rst_file_name, site_indices=None):
        ancestral_seq = ""
        ancestral_num = -1
        rst_file = open(rst_file_name, "r")
        flag = False
        for line in rst_file.readlines():
            if flag == False and line.startswith("Nodes"):
                elements = line.split()
                ancestral_num = elements[1]
                flag = True
            elif flag == True and line.startswith(f"node #{ancestral_num}"):
                ancestral_seq = "".join(line.split()[2:])
        if site_indices is None:
            return ancestral_seq
        else:
            return "".join(list(numpy.array(list(ancestral_seq))[site_indices]))

    def read_rate(self, file_name, site_indices=None):
        rates = []
        flag = False
        with open(file_name) as f:
            for line in f:
                cols = line.strip().split()
                if (flag == True) and (len(cols) == 5):
                    rates.append(float(cols[3]))
                if 'posterior' in line:
                    flag = True
        if site_indices is None:
            return rates
        else:
            return list(numpy.array(rates)[site_indices])

    def read_mat(self, mat_file):
        mat = numpy.zeros((len(self.aas), len(self.aas)))
        string_list = open(mat_file).read().split()
        i = 0
        for j in range(len(self.aas)):
            for k in range(j):
                mat[j, k] = float(string_list[i])
                i += 1
        mat += mat.T
        return mat

    def get_freqs(self, seqs, freq_mode='gene'):
        if freq_mode == 'gene':
            freq = []
            seq = ''.join(seqs)
            for aa in self.aas:
                freq.append(seq.count(aa))
            freq = numpy.array(freq)
            return freq / freq.sum()
        elif freq_mode == 'site':
            d = dict(zip(self.aas, list(range(len(self.aas)))))
            freq = numpy.zeros((len(seqs[0]), len(self.aas)))
            for seq in seqs:
                for i in range(len(seq)):
                    if seq[i] in list(self.aas):
                        freq[i, d[seq[i]]] += 1.
            return freq / freq.sum(axis=1, keepdims=True)
        else:
            return numpy.loadtxt(freq_mode)


    def get_inferred_tree(self, gene_in_dir, gene_out_dir):
        with open(f'{gene_in_dir}/rst') as rst_f:
            content = rst_f.readlines()
            rst_f.close()
        os.makedirs(gene_out_dir, exist_ok=True)
        with open(f'{gene_out_dir}/tree_inferred_branch.nwk', 'w') as f:
            f.write(content[7][:-1])
            f.close()

class Simulators:

    def __init__(self, args):
        self.args = args

    def calc_expconv(self, seqID):
        cec = CalcExpConv(self.args)
        ancestral_seq, rates, mat, freq, tree_result_path = cec.main(seqID)
        return ancestral_seq, rates, mat, freq, tree_result_path

    # Use the key info from paml result file to generate sequences without convergence
    def simulate(self, seqID):
        ancestral_seq, rates, mat, freq, tree_result_path = self.calc_expconv(seqID)
        simulator = sequence.HeterogeneousProteinSequenceSimulator(
            rates=rates,
            mat=mat,
            ancestral_seq=ancestral_seq,
            freq=freq,
            substitution_model=None,
            profile=None,
            profile_resampler=None,
            heterogeneous_branch_ratio=-1,
            rate_swap_ratio='random',
            profile_swap_model=None,
        )
        # Sequence Simulations
        tree = Tree(tree_result_path)
        tree = simulator.generate_sequences(tree)
        return tree

    def write_seq(self, seqID, cycle_time):
        tree = self.simulate(seqID)
        simulate_sequence_path = f"{self.args.dataset_location}/{seqID}/{seqID}_SimSequences"
        os.makedirs(simulate_sequence_path, exist_ok=True)
        f = open(f"{simulate_sequence_path}/{seqID}_{cycle_time}.fasta", "w")
        for node in tree.iter_leaves():
            print(">" + node.name, file=f)
            print(node.sequence, file=f)
        f.close()

if __name__ == '__main__':
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    Sim = Simulators(args)
    for cycle_time in range(1, args.simtime+1):
        Sim.write_seq(args.seqid, cycle_time)
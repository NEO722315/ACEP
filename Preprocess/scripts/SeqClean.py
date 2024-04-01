import os
import sys
import numpy
sys.path.append("../../")
import args_prepare

def parse_prepare():
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    return args

class PreProcessor:

    def __init__(self, args):
        self.args = args
        self.dataset_location = self.args.dataset_location
        self.rawAAdir = self.args.rawAAdir
        self.rawCDSdir = self.args.rawCDSdir
        self.case_sps1 = self.args.case_sps1
        self.case_sps2 = self.args.case_sps2
        self.aas = list('ARNDCQEGHILKMFPSTWYV')

    def Clean(self, aa_fasta_file, cds_fasta_file, aa_output_file, cds_output_file):
        content = open(aa_fasta_file, "r").readlines()
        seqid = os.path.basename(aa_fasta_file).split('.')[0]
        spnames = [spname[1:-1] for spname in content[::2]]
        seqs = [seq[:-1] for seq in content[1::2]]
        # judge if there are focus species in both groups
        if not (len(set(spnames).intersection(self.case_sps1)) > 1 and
            len(set(spnames).intersection(self.case_sps2)) > 1):
            return
        orignal_length = numpy.mean([len(seq) - seq.count('-') for seq in seqs])
        new_seqs, new_spnames = [], []
        for i, spname in enumerate(spnames):
            # exclude the sequence whose proportion of gaps exceed 5%
            if seqs[i].count('-') > 0.05 * orignal_length:
                continue
            new_spnames.append(spname)
            new_seqs.append(seqs[i])
        if not (len(set(new_spnames).intersection(self.case_sps1)) > 1 and
            len(set(new_spnames).intersection(self.case_sps2)) > 1):
            return
        original_shape, clean_seq_array, gapcol_indices = self.clean_alignment_aa(new_seqs)
        if clean_seq_array == "null":
            return
        clean_cds_seqs = self.clean_alignment_cds(cds_fasta_file, original_shape, gapcol_indices, new_spnames)
        self.save_file(seqid, clean_seq_array, new_spnames, aa_output_file)
        self.save_file(seqid, clean_cds_seqs, new_spnames, cds_output_file)

    def clean_alignment_aa(self, seqs):
        seq_array = numpy.array([list(seq) for seq in seqs])
        original_shape = seq_array.shape
        ins = numpy.isin(seq_array, self.aas)
        gap_indices = numpy.where(ins == False)
        gapcol_indices = list(set(gap_indices[1]))
        clean_seq_array = numpy.delete(seq_array, gapcol_indices, axis=1)
        if clean_seq_array.shape[1] > 1024:
            pos = [(i+1)*1024 for i in range(clean_seq_array.shape[1] // 1024)]
            split_seqs = numpy.split(clean_seq_array, (*pos,), axis=1)
        elif clean_seq_array.shape[1] < 100:
            split_seqs = "null"
        else:
            split_seqs = [clean_seq_array]
        return original_shape, split_seqs, gapcol_indices

    def clean_alignment_cds(self, fasta_file, original_shape, gapcol_indices, new_spnames):
        cds_content = open(fasta_file, "r").readlines()
        cds_seqs = numpy.array([[seq[:-1][i:i+3] for i in range(0, len(seq[:-1]), 3)] for seq in cds_content[1::2]])
        spnames_indices = [i for i, spname in enumerate(cds_content[::2]) if spname[1:-1] in new_spnames]
        keepcol_indices = [i for i in range(original_shape[1]) if not i in gapcol_indices]
        clean_cds_seqs = cds_seqs[:, keepcol_indices]
        clean_cds_seqs = clean_cds_seqs[spnames_indices, :]
        return clean_cds_seqs

    def save_file(self, seqid, seq_array, spnames, outputfile):
        if outputfile.endswith("aa.cleaned.fasta"):
            for i in range(len(seq_array)):
                dir_path = f"{self.dataset_location}/01_AA_data/{seqid}_{i}"
                os.makedirs(dir_path, exist_ok=True)
                outputpath = f"{dir_path}/{outputfile}"
                file = open(outputpath.replace('.aa.cleaned', f'_{i}.aa.cleaned'), "w")
                for m in range(seq_array[i].shape[0]):
                    seq = "".join(list(seq_array[i][m]))
                    print(f">{spnames[m]}", file=file)
                    print(seq, file=file)
                file.close()
        else:
            dir_path = f"{self.dataset_location}/02_CDS_data/{seqid}"
            os.makedirs(dir_path, exist_ok=True)
            file = open(f"{dir_path}/{outputfile}", "w")
            for i in range(seq_array.shape[0]):
                seq = "".join(list(seq_array[i]))
                print(f">{spnames[i]}", file=file)
                print(seq, file=file)
            file.close()


    def main(self):
        for file in os.listdir(self.rawAAdir):
            seqid = file.split('.')[0]
            print(seqid)
            os.makedirs(self.dataset_location, exist_ok=True)
            self.Clean(f"{self.rawAAdir}/{file}", f"{self.rawCDSdir}/{file}",
                       f"{seqid}.aa.cleaned.fasta", f"{seqid}.cds.cleaned.fasta")

if __name__ == '__main__':
    args = parse_prepare()
    prep = PreProcessor(args)
    prep.main()

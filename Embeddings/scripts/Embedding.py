import os
import torch
from torch.nn import functional as F
import sys
from esm import pretrained
sys.path.append("../../")

import args_prepare

def parse_prepare():
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    return args

class Embedding:

    def __init__(self, args):
        self.args = args
        self.model, self.alphabet = self._load_model()

    def _load_model(self):
        # load model and the alphabet encoder
        model, alphabet = pretrained.load_model_and_alphabet(f"{self.args.abspath}/{self.args.model_location}")
        model.eval()
        model = model.cuda(device=self.args.device)
        return model, alphabet

    def SingleCombine(self):
        os.makedirs(self.args.combine_dir, exist_ok=True)
        combinefile = open(f"{self.args.combine_dir}/{self.args.seqid}.fasta", "a")
        self.args.fasta_file = f"{self.args.combine_dir}/{self.args.seqid}.fasta"
        for seq_file in os.listdir(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_SimSequences"):
            content = open(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_SimSequences/{seq_file}", "r").readlines()
            spnames = content[::2]
            seqcontent = content[1::2]
            for i, spname in enumerate(spnames):
                file_chuck = seq_file.split('.')[0].rsplit("_", 1)
                print(f">{file_chuck[0]}|{file_chuck[-1]}|{spname[1:-1]}", file=combinefile)
                print(seqcontent[i][:-1], file=combinefile)
        ori_content = open(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}.aa.cleaned.fasta", "r").readlines()
        spnames = ori_content[::2]
        seqcontent = ori_content[1::2]
        for i, spname in enumerate(spnames):
            print(f">{self.args.seqid}|ORIGINAL|{spname[1:-1]}", file=combinefile)
            print(seqcontent[i][:-1], file=combinefile)
        combinefile.close()

    def MultiEmbed(self):
        print("Transferred model to GPU")
        if "bottleneck" in self.args.include:
            linear1_weight = torch.load(f"{self.args.abspath}/{self.args.bottleneck_weight}", map_location=torch.device(self.args.device)).T
            print("linear1 weight", linear1_weight.shape)
        content = open(self.args.fasta_file, "r").readlines()
        names, seqs = [name[1:-1] for name in content[::2]], [seq[:-1] for seq in content[1::2]]
        data = zip(names, seqs)
        print(f"Read {self.args.fasta_file} and simulated sequences")
        with torch.no_grad():
            for idx, (name, seq) in enumerate(data):
                pure_name = name.split('|')[-1]
                if self.args.embed_single:
                    if not pure_name in self.args.case_sps1 and not pure_name in self.args.case_sps2:
                        continue
                toki = self.alphabet.encode('<cls>' + seq + '<eos>')
                toki = torch.tensor(toki).unsqueeze(0).unsqueeze(0)
                toki_c = toki.to(self.args.device)
                out = self.model(toki_c, repr_layers=self.args.repr_layers)
                emb = out['representations'][12]
                emb = emb.reshape((1, *emb.shape[-2:]))
                if "raw" in self.args.include:
                    result = {"raw_representations": {12: emb.cpu().detach()}}
                if "mean" in self.args.include:
                    result = {"mean_representations": {12: torch.mean(emb, dim=1).cpu().detach()}}
                if "bottleneck" in self.args.include:
                    emb = F.pad(emb, (0, 0, 0, 1024 - emb.shape[-2]))
                    emb = emb.reshape(emb.shape[0], -1)
                    result = {"bottleneck_representations": {12: torch.matmul(emb, linear1_weight).cpu().detach().squeeze()}}
                self.Distribute(result, name)

    def Distribute_func(self, label, value, label_elements):
        if isinstance(label_elements, list):
            if "ORIGINAL" == label_elements[1]:
                hypoth = "alter"
            else:
                hypoth = "null"
        emb_path = f"{self.args.dataset_location}/{label_elements[0]}/{label_elements[0]}_{hypoth}_Embedding"
        os.makedirs(emb_path, exist_ok=True)
        torch.save({label: value}, f"{emb_path}/{label}|{self.args.include[0]}.pt")

    def Distribute(self, value, label):
        # Every Pt File contains the whole simulation sequences which belong to multiple genes
        if '|' in label:
            label_elements = label.split("|")
        else:
            label_elements = [label]
        if self.args.embed_single:
            if label_elements[-1] in self.args.case_sps1 or label_elements[-1] in self.args.case_sps2:
                self.Distribute_func(label, value, label_elements)
        else:
            self.Distribute_func(label, value, label_elements)


if __name__ == '__main__':
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    Embedding = Embedding(args)
    for seqid in os.listdir(args.dataset_location):
        Embedding.args.seqid = seqid
        Embedding.SingleCombine()
        Embedding.MultiEmbed()

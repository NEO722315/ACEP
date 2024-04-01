import os
import pickle as pkl
import torch
import numpy as np
from numpy.linalg import norm
import sys

sys.path.append("../../")
import args_prepare

class ConvCalculation:

    def __init__(self, args):
        self.args = args
        self.allsps = [sp[1:-1] for sp in open(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}.aa.cleaned.fasta", "r").readlines()[::2]]

    def DistanceCal(self, vector_a, vector_b):
        vector_a, vector_b = vector_a.reshape(-1), vector_b.reshape(-1)
        cos_dis = 1.0 - self.cosine_distance(vector_a, vector_b)
        eu_dis = np.sqrt(np.sum((vector_a - vector_b) ** 2))
        corr_dis = 1.0 - np.corrcoef(vector_a, vector_b)[1, 0]
        hyper_dis = self.hyperbolic_dist(vector_a, vector_b)
        dis_res = {"cosine distance": cos_dis, "euclidean distance": eu_dis, "correlation distance": corr_dis, "hyperbolic distance": hyper_dis}
        return dis_res

    def cosine_distance(self, Vec_a, Vec_b):
        x1 = Vec_a / norm(Vec_a, axis=-1, keepdims=True)
        x2 = Vec_b / norm(Vec_b, axis=-1, keepdims=True)
        cos = np.dot(x1, x2.T)
        return cos

    def hyperbolic_dist(self, u, v, c=-1):
        return np.arccosh(
            1.0 - 2.0 * c * (np.sum((u - v) ** 2.0)) / (1.0 + c * np.sum(u ** 2.0)) / (1.0 + c * np.sum(v ** 2.0)))

    def FocusSpecies(self):
        focus_case_sps1, focus_case_sps2 = [], []
        spnames = [spname[1:-1] for spname in open(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}.aa.cleaned.fasta").readlines()[::2]]
        for spname in spnames:
            if spname in self.args.case_sps1:
                focus_case_sps1.append(spname)
            if spname in self.args.case_sps2:
                focus_case_sps2.append(spname)
        return focus_case_sps1, focus_case_sps2

    # extract the focus species are included in this MSA
    def SpeciesGroupPair(self, group1, group2):
        sps_pairs = []
        for spname1 in group1:
            for spname2 in group2:
                sps_pairs.append((spname1, spname2))
        return sps_pairs

    def NullDistribution(self):
        focus_case_sps1, focus_case_sps2 = self.FocusSpecies()
        case_sps_num = len(focus_case_sps1) + len(focus_case_sps2)
        nullseqs_num = len([file for file in os.listdir(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_SimSequences")])
        if self.args.embed_single:
            embedding_num = len([file for file in os.listdir(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_Embedding")]) // case_sps_num
        else:
            sequence_num = len(open(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}.aa.cleaned.fasta", "r").readlines()) // 2
            embedding_num = len([file for file in os.listdir(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_Embedding")]) // sequence_num
        assert nullseqs_num == self.args.simtime
        assert embedding_num == nullseqs_num
        sps_pairs = self.SpeciesGroupPair(focus_case_sps1, focus_case_sps2)
        all_results = []
        for cycle_time in range(1, nullseqs_num + 1):
            result = {}
            for pair in sps_pairs:
                emb_a = torch.load(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_Embedding/{self.args.seqid}|{cycle_time}|{pair[0]}|{self.args.include[0]}.pt")
                emb_b = torch.load(f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_Embedding/{self.args.seqid}|{cycle_time}|{pair[1]}|{self.args.include[0]}.pt")
                emb_a = emb_a[f"{self.args.seqid}|{cycle_time}|{pair[0]}"][f"{self.args.include[0]}_representations"][self.args.repr_layers[0]]
                emb_b = emb_b[f"{self.args.seqid}|{cycle_time}|{pair[1]}"][f"{self.args.include[0]}_representations"][self.args.repr_layers[0]]
                emb_a = emb_a.numpy()
                emb_b = emb_b.numpy()
                distances = self.DistanceCal(emb_a, emb_b)
                result[f"{pair[0]}&{pair[1]}&{cycle_time}"] = distances
            all_results.append(result)
        pkl_file = open(
            f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_distribution.pickle", "wb")
        pkl.dump(all_results, pkl_file)
        pkl_file.close()

    def AlterDistribution(self):
        alter_dir = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_alter_Embedding"
        pkl_file = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_alter_distribution.pickle"
        focus_case_sps1, focus_case_sps2 = self.FocusSpecies()
        sps_pairs = self.SpeciesGroupPair(focus_case_sps1, focus_case_sps2)
        result = {}
        for pair in sps_pairs:
            emb_a = torch.load(f"{alter_dir}/{self.args.seqid}|ORIGINAL|{pair[0]}|{self.args.include[0]}.pt")
            emb_b = torch.load(f"{alter_dir}/{self.args.seqid}|ORIGINAL|{pair[1]}|{self.args.include[0]}.pt")
            emb_a = emb_a[f"{self.args.seqid}|ORIGINAL|{pair[0]}"][f"{self.args.include[0]}_representations"][
                self.args.repr_layers[0]]
            emb_b = emb_b[f"{self.args.seqid}|ORIGINAL|{pair[1]}"][f"{self.args.include[0]}_representations"][
                self.args.repr_layers[0]]
            emb_a = emb_a.numpy()
            emb_b = emb_b.numpy()
            distances = self.DistanceCal(emb_a, emb_b)
            result[f"{pair[0]}&{pair[1]}&ORIGINAL"] = distances
        pkl_file = open(pkl_file, "wb")
        pkl.dump([result], pkl_file)
        pkl_file.close()

class ConvStat:

    def __init__(self, args):
        self.args = args
        self.metrics = ["cosine distance", "euclidean distance", "correlation distance", "hyperbolic distance"]

    def load_data(self, NA):
        path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_{NA}_distribution.pickle"
        distances = pkl.load(open(path, "rb"))
        mat = np.zeros((len(distances), len(distances[0]), len(self.metrics)))
        for j, cycle in enumerate(distances):
            for i, (pair_label, pair) in enumerate(cycle.items()):
                for m, (metric, value) in enumerate(pair.items()):
                    mat[j, i, m] = value
        np.save(f"{path.rsplit('.', 1)[0]}.npy", mat)

    def cal_mean_pvalue(self):
        null_path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_distribution.npy"
        alter_path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_alter_distribution.npy"
        null_mat, alter_mat = np.load(null_path), np.load(alter_path)
        length = null_mat.shape[0]
        null_mean_mat, alter_mean_mat = np.mean(null_mat, axis=1, keepdims=True), np.mean(alter_mat, axis=1, keepdims=True)
        mat = null_mean_mat - alter_mean_mat
        mat = mat < 0
        prop_mat = np.sum(mat, axis=0, keepdims=True)
        return prop_mat / length

    # def cal_meadian_pvalue(self):
    #     null_path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_null_distribution.npy"
    #     alter_path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_alter_distribution.npy"
    #     null_mat, alter_mat = np.load(null_path), np.load(alter_path)
    #     length = null_mat.shape[0]
    #     null_median_mat, alter_median_mat = np.median(null_mat, axis=1, keepdims=True), np.median(alter_mat, axis=1, keepdims=True)
    #     mat = null_median_mat - alter_median_mat
    #     mat = mat < 0
    #     prop_mat = np.sum(mat, axis=0, keepdims=True)
    #     return prop_mat / length

    def output_results(self):
        self.load_data('null')
        self.load_data("alter")
        mean_result = self.cal_mean_pvalue()
        # median_result = self.cal_meadian_pvalue()
        mean_path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_mean_convresult.pickle"
        # median_path = f"{self.args.dataset_location}/{self.args.seqid}/{self.args.seqid}_median_convresult.pickle"
        mean_file = open(mean_path, "wb")
        # median_file = open(median_path, "wb")
        pkl.dump(['mean', self.metrics, mean_result.squeeze().tolist()], mean_file)
        # pkl.dump(['median', self.metrics, median_result.squeeze().tolist()], median_file)
        mean_file.close()
        # median_file.close()

if __name__ == '__main__':
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    for seqid in os.listdir(args.dataset_location):
        args.seqid = seqid
        ConvCal = ConvCalculation(args)
        ConvCal.NullDistribution()
        ConvCal.AlterDistribution()
        ConvS = ConvStat(args)
        ConvS.output_results()

import os
import args_prepare
from Preprocess.scripts.SeqClean import PreProcessor
from PhyloInfer.scripts.pamlparas import Fileprepare
from Simulations.scripts import Simulation,sequence
from Embeddings.scripts import Embedding
from Calculation.scripts import ConvStatistic

def SingleConv(args):
    # Single gene Convergence Detection
    sim_path = f"{args.dataset_location}/{args.seqid}/{args.seqid}_SimSequences"
    comb_path = f"{args.combine_dir}/{args.seqid}.fasta"
    if not os.path.exists(sim_path):
        print("Simulating ... ...")
        Sim = Simulation.Simulators(args)
        for cycle_time in range(1, args.simtime + 1):
            Sim.write_seq(args.seqid, cycle_time)
        print("Simulation Done")
    if os.path.exists(comb_path):
        os.remove(comb_path)
    print("Embedding sequences ... ...")
    Embedder = Embedding.Embedding(args)
    Embedder.SingleCombine()
    Embedder.MultiEmbed()
    print("Embedding Done")
    print("Statistical inference of convergent signal ... ...")
    ConvCal = ConvStatistic.ConvCalculation(args)
    ConvCal.NullDistribution()
    ConvCal.AlterDistribution()
    ConvS = ConvStatistic.ConvStat(args)
    ConvS.output_results()
    print("All Done")

def Parapre(args):
    # preprocess the raw sequences dataset
    print("Preprocessing of raw sequences ... ...")
    prep = PreProcessor(args)
    prep.main()
    print("Preprocessing of raw sequences Done.")
    # prepare the paml inference files
    args.dataset_location = f"{args.dataset_location}/01_AA_data"
    print("Paml parameters inferencing ... ...")
    fp = Fileprepare(args)
    fp.main()
    print("Paml parameters inference Done")


if __name__ == '__main__':
    parser = args_prepare.create_parser()
    args = parser.parse_args()
    Parapre(args)
    args.dataset_location = f"{args.abspath}/{args.dataset_location}"
    for seqid in os.listdir(args.dataset_location):
        args.seqid = seqid
        SingleConv(args)

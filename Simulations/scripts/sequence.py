import abc
import random
import numpy
import math
import scipy.linalg
import scipy.stats
import sys

sys.path.append("./")
from codonZ import AA_std_table, reverse_AA_table

_DT = 0.01

AMINO_ACIDS = numpy.fromstring('ARNDCQEGHILKMFPSTWYV', dtype='S1')
NUCLEOTIDES = numpy.fromstring('AGCT', dtype='S1')


class BaseEvolutionSimulator(metaclass=abc.ABCMeta):
    """An abstract sequence simulator"""

    def __init__(self, rates, mat, ancestral_seq, freq, substitution_model):
        """Create a BaseSequenceSimualtor

        Parameters
        ----------
        site_count_range: tuple[int]
            The upper and lower bounds of subsititution sites
        subsititution_model: str
            A subsititution model
        """
        self.site_count = len(rates)
        self.rates = numpy.array(rates)
        self.substitution_model = substitution_model
        self.mat = mat
        self.ancestral_seq = ancestral_seq
        self.freq = freq


    def generate_sequences(self, tree, leaves=None):
        """Generate sequences and assign them to the tree nodes.

        Parameters
        -----------
        tree : ete3.PhyloTree
            A phylogenetic tree.
        leaves: List[ete3.PhyloNode] | NoneType
            Limiting simulation to the specified list of leaf nodes,
            ignoring all other leaves. If it is None, simulate all
            nodes. Default is None.
        """
        parameters = self.prepare_parameters()
        # 将所有指定叶节点及其父节点也加入，一直回溯到祖先
        # 将所有这些节点进行序列的模拟
        if leaves is not None:
            simulating_nodes = set()
            for leaf in leaves:
                node = leaf
                while node is not None:
                    simulating_nodes.add(node)
                    node = node.up
        else:
            simulating_nodes = None
        # 对每一棵树进行的每个节点模拟序列
        for node in tree.traverse('preorder'):
            if simulating_nodes is not None and node not in simulating_nodes:
                continue
            if node == tree:
                # 如果为根节点则初始化一条序列
                self._initialize_root(node, tree, self.ancestral_seq)
            else:
                # 如果是子代节点则生成序列
                self._simulate_branch(node, **parameters)
        self.decode_sequences(tree)
        # self.clean_up(tree)
        return tree


    def prepare_parameters(self):
        """Prepare parameters for evolution simulation.

        Returns
        -------
        dict
            A dictionary of simulation parameters
        """
        # 生成前景枝与背景枝所用矩阵
        # q = self._calculate_q_matrix(self.mat, self.freq)
        # p = self._calculate_p_matrix(q)
        parameters = {
            # 在site_count_range范围内任意选取数字代表序列长度
            'site_count': self.site_count,
            'r': self.mat,
            # 'q': q,
            # 'p': p,
            'equilibrium_frequencies': self.freq,
        }
        return parameters


    @abc.abstractmethod
    def load_substitution_model(self, model):
        """Load the specified substitution model

        Parameters
        ----------
        model : str
            The name of the substitution model

        Returns
        -------
        numpy.ndarray[float]
            The R matrix
        numpy.ndarray[float]
            The equilibrium frequency matrix
        """


    @classmethod
    def _calculate_q_matrix(cls, r, equilibrium_frequencies):
        """Generate Q matrix from R matrix

        Parameters
        ----------
        r:numpy.ndarray
            R matrix, symmetric exchange rates
        equilibrium_frequencies: numpy.ndarray
            specific equilibrium frequency vector

        Returns
        -------
        numpy.ndarray[float]
            Q matrix in sequence evolution
        """
        q = r * equilibrium_frequencies
        numpy.fill_diagonal(q, 0)
        numpy.fill_diagonal(q, -q.sum(axis=1))
        q /= numpy.abs((equilibrium_frequencies * numpy.diag(q)).sum())
        return q


    @classmethod
    def _calculate_p_matrix(cls, q):
        """Generate P matrix

        Parameters
        ----------
        q: numpy.ndarray
            instant rate matrix Q

        Returns
        -------
        numpy.ndarray
            P (transition) matrix in sequence evolution
        """
        return scipy.linalg.expm(q * _DT)


    @classmethod
    @abc.abstractmethod
    def _initialize_root(cls, node, tree, ancestral_seq):
        """Generate a sequence for the tree root.

        Parameters
        ----------
        tree : ete3.PhyloTree
            A phylogenetic tree.
        """


    @classmethod
    @abc.abstractmethod
    def _simulate_branch(self, node, **__):
        """Simulate a sequence for a non-root tree node.

        Parameters
        ----------
        node : ete3.PhyloTreeNode
            A phylogenetic tree
        """


    @classmethod
    @abc.abstractmethod
    def decode_sequences(cls, tree):
        """Encode all sequences in a tree to bytes.

        Parameters
        ----------
        tree : ete3.PhyloTree
            A phylogenetic tree with all sequences simulated.
        """


    @classmethod
    def clean_up(cls, tree):
        """Remove extra attributes"""


class BranchEvolutionSimulator(BaseEvolutionSimulator,
                                            metaclass=abc.ABCMeta):
    """A sequence simulator for heterogeneous branches"""

    def __init__(self, *args, profile, profile_resampler=None,
                 heterogeneous_branch_ratio="random", rate_swap_ratio='random',
                 profile_swap_model=None, **kwargs):
        """Create a BaseSequenceSimulator.

                Parameters
                ----------
                site_count_range : tuple[int]
                    The upper and lower bounds of substitution sites.
                substitution_model : str
                    A substitution model.
                profile : pathlib.Path | str
                    The frequency profile file or name.
                profile_resampler : tuple[str, int] | NoneType
                    If it is None, no profile resampling will be performed.
                    Otherwise a new profile matrix is generated by resampling.
                heterogeneous_branch_ratio : float | str
                    The probability to regenerate profile matrices and
                    site-specific evolution rates on a branch. If it is
                    "random", a random number between 0 and 1 is used instead.
                    Default is "random".
                rate_swap_ratio : float | str
                    The fraction of site-specific rates to shuffle on a branch.
                    If it is "random", a random number between 0 and 1 is used
                    instead. Default is "random".
                profile_swap_model : scipy.stats.rv_discrete | None
                    A statistical model to generate the number of swapping
                    site-specific profile elements.
        """
        super().__init__(*args, **kwargs)
        self.profile = profile
        self.profile_resampler = profile_resampler
        self.heterogeneous_branch_ratio = heterogeneous_branch_ratio
        self.rate_swap_ratio = rate_swap_ratio
        self.profile_swap_model = profile_swap_model


    def prepare_parameters(self):
        """Prepare parameters for evolution simulation.

        Returns
        -------
        dict
            A dictionary of simulation parameters.
        """
        parameters = super().prepare_parameters()
        shape = 20
        site_count = parameters['site_count']
        equilibrium_frequency = parameters['equilibrium_frequencies']
        # 对profile文件中加入噪声以及调换平衡频率排列顺序
        profile_content = self._generate_profile(
            site_count, shape, equilibrium_frequencies=equilibrium_frequency,
            profile=self.profile, resampler=self.profile_resampler,
        )
        # evolutionary rate
        rates = self.rates
        heterogeneous_branch_ratio = self.heterogeneous_branch_ratio
        if heterogeneous_branch_ratio == 'random':
            heterogeneous_branch_ratio = random.random()
        parameters.update({
            'rates': rates,
            'profile': profile_content,
            'heterogeneous_branch_ratio': heterogeneous_branch_ratio,
            'rate_swap_ratio': self.rate_swap_ratio,
            'profile_swap_model': self.profile_swap_model
        })
        return parameters


    """每次提取profile时生成不一样的文件(含噪声以及顺序不同)"""
    @classmethod
    def _generate_profile(cls, site_count, shape, *, equilibrium_frequencies,
                          profile, resampler):
        """Create a BaseSequenceSimulator.

        Parameters
        ----------
        site_count : int
            The number of substitution sites.
        equilibrium_frequencies : numpy.ndarray[float]
            A equilibrium frequency matrix.
        profile : pathlib.Path | str
            The frequency profile file or name.
        resampler : tuple[str, int] | NoneType
            If it is None, no profile resampling will be performed.
            Otherwise a new profile matrix is generated by resampling.
        """
        if profile is None:
            profile = numpy.array(equilibrium_frequencies).reshape(-1, shape)
            # 出错bug原因：repeat函数中重复次数只能保留整数部分,但在此例必须向上取整达到sequence长度与rates长度保持一致
            return numpy.repeat(profile, math.ceil(site_count / profile.shape[0]), axis=0)[0:site_count]
        # 获取profile路径
        path = f"../heterogeneity_profiles/{profile}"
        count_array = numpy.loadtxt(path)
        # 判断是否只有一条平衡频率
        if len(count_array.shape) == 1:
            count_array = count_array.reshape(1, -1)
        # 判断是否进行了归一化操作
        if count_array.max() >= 1.0001:
            # 对其中的0进行非0化处理
            count_array.clip(numpy.finfo('f4').eps)
        # 归一化处理
        profile = count_array / count_array.sum(axis=1, keepdims=True)
        # 随机取出多个平衡频率并赋值给profile
        sample = numpy.random.choice(profile.shape[0], site_count)
        profile = profile[sample]
        # 添加噪音
        if resampler is None:
            return profile
        elif not isinstance(resampler, tuple) or len(resampler) != 2:
            raise ValueError(f"invalid resampler: {resampler}")
        elif resampler[0] == 'multinomial':
            for i in range(profile.shape[0]):
                sample = numpy.random.multinomial(
                    resampler[1], profile[i, :]
                ) + resampler[1] / 1e5
                sample /= sample.sum(axis=1)
                profile[i, :] = sample
            return profile
        elif resampler[0] == 'dirichlet':
            # 使得每一列均值保持不变，但是服从gamma分布
            sample = numpy.random.gamma(profile * resampler[1])
            sample = sample / sample.sum(axis=1, keepdims=True)
            sample = sample.clip(numpy.finfo("f4").eps)
            return sample / sample.sum(axis=1,keepdims=True)
        else:
            raise ValueError(f"invalid resampler: {resampler}")


    @classmethod
    def _initialize_root(cls, node, tree, ancestral_seq):
        """Generate a sequence for the tree root.

        tree: ete3.PhyloTree
            A phylogenetic tree.
        profile: numpy.ndarray[float]
            The frequency profile matrix
        """
        if ancestral_seq is None:
            tree.add_feature('sequence', _sample_sequences(parameters['profile']))
        node.add_feature('sequence', numpy.array([AA_std_table[i] for i in list(ancestral_seq)]))

    @classmethod
    def _simulate_branch(cls, node, **parameters):
        """Simulate a sequence for a non-root tree node.

        Parameters
        ----------
        node : ete3.PhyloNode
            A phylogenetic tree.
        rates : numpy.ndarray[float]
            The site-specific evolution rates.
        profile : numpy.ndarray[float]
            The frequency profile matrix.
        heterogeneous_branch_ratio : float
            The probability to regenerate profile matrices and
            site-specific evolution rates on a branch.
        rate_swap_ratio : float
            The fraction of site_specific rates to shuffle on a branch.
        profile_swap_model : scipy.stats.rv_discrete | None
            A statistical model to generate the number of swapping
            site-specific profile elements. If it is None, no swap is performed.
        """
        # 获取父节点的序列
        sequence = node.up.sequence.copy()
        heterogeneous_branch_ratio = parameters['heterogeneous_branch_ratio']
        # 使得每个位点上生成的枝长正比于演化速率
        t = node.dist * cls._simulate_branch_rates(
            parameters['rates'], heterogeneous_branch_ratio,
            parameters['rate_swap_ratio'],
        )
        mask = t > 0
        mask_count = mask.sum()
        if mask_count == 0:
            node.add_feature('sequence', sequence)
            return
        # 氨基酸使用频率变换
        profile = cls._generate_branch_profile(
            parameters['profile'], heterogeneous_branch_ratio,
            parameters['profile_swap_model']
        )
        # 跳跃式演化策略
        indices = numpy.arange(profile.shape[1])
        q = parameters['r'] * profile[:, None, :]
        q += numpy.finfo('float32').eps
        q_diag = q[:, indices, indices] - q.sum(axis=2)
        q[:, indices, indices] = q_diag
        scale_f = - (q_diag * profile).sum(axis=-1)
        q = q / scale_f[:, None, None]
        q_diag = q.sum(axis=2) - q[:, indices, indices]
        # q_diag.shape  (4986, 20, )
        jump = q / q_diag[..., None]
        jump[:, indices, indices] = 0
        # 减少的时间为 1/qi
        t[mask] -= numpy.random.exponential(numpy.reciprocal(q_diag[mask, sequence[mask]]))
        mask = t > 0
        mask_count = mask.sum()
        while mask_count > 0:
            # 筛取mask>0位点的1 * 20的矩阵,做sample_sequence运算
            sequence[mask] = _sample_sequences(jump[mask, sequence[mask]])
            t[mask] -= numpy.random.exponential(numpy.reciprocal(q_diag[mask, sequence[mask]]))
            mask = t > 0
            mask_count = mask.sum()
        node.add_feature('sequence', sequence)


    @classmethod
    def _simulate_branch_rates(cls, rates, heterogeneous_branch_ratio,
                               swap_ratio):
        """Create a BaseSequenceSimulator. rates distributed as gamma distribution.

            Parameters
            ----------
            rates : numpy.ndarray[float]
                The site-specific evolution rates.
            heterogeneous_branch_ratio : float
                The probability to regenerate profile matrices and
                site-specific evolution rates on a branch.
            swap_ratio : float
                The fraction of site-specific rates to shuffle on a branch.
        """
        # 加入位点异质性
        rates = rates.copy()
        if random.random() < heterogeneous_branch_ratio:
            print("heterogeneous branch")
            if swap_ratio == 'random':
                swap_ratio = random.random()
                print("swap branch")
            length = rates.shape[0]
            # 抽取length * swap_ratio个数的位点的演化速率，无放回
            mask = numpy.random.choice(length, round(length * swap_ratio),
                                       replace=False)
            # 将速率值随机排列后放回去
            rates[mask] = numpy.random.permutation(rates[mask])
        return rates


    # 氨基酸使用频率差异
    @classmethod
    def _generate_branch_profile(cls, profile, heterogeneous_branch_ratio,
                                 swap_model):
        """Generate a profile matrix for a branch.

        Parameters
        ----------
            profile : numpy.ndarray[float]
                The frequency profile matrix.
            heterogeneous_branch_ratio : float
                The probability to regenerate profile matrices and
                site-specific evolution rates on a branch.
            swap_model : scipy.stats.rv_discrete | None
                A statistical model to generate the number of swapping
                site-specific profile elements. If it is None, no swap is
                performed.
        """
        profile = profile.copy()
        # print(profile, profile.shape)      (608, 20)
        if (swap_model is not None
                and random.random() < heterogeneous_branch_ratio):
            # 生成每个平衡频率的氨基酸频率所要调换的次数
            swap_count = swap_model.rvs(profile.shape[0])
            # 交换利用cython
            # _common.repeated_swap(profile, swap_count)
        return profile

class ProteinSimulatorMixin:
    """A protein simulator mix-in."""

    _SUBSTITUTION_MODELS = ('cpREV64', 'dayhoff', 'jones', 'lg', 'wag',
                            'mtArt', 'mtmam', 'mtREV24', 'MtZoa')

    def load_substitution_model(self, name):
        """Read a PAML amino acid substitution model.

        Parameters
        ----------
        name: str
            A PAML substitution model name.

        Returns
        -------
        numpy.ndarray[float]
            The R matrix
        numpy.ndarray[float]
            The equilibrium frequency matrix
        """
        if name == 'random':
            name = random.choice(self._SUBSTITUTION_MODELS)
        codon_count = AMINO_ACIDS.size
        path = f'../substitution_models/{name}.dat'
        string_list = open(path).read().split()
        i = 0
        r = numpy.zeros((codon_count, codon_count))
        for j in range(codon_count):
            for k in range(j):
                r[j, k] = float(string_list[i])
                i += 1
        r = r + r.T
        equilibrium_frequencies = numpy.zeros(codon_count)
        for j in range(codon_count):
            try:
                equilibrium_frequencies[j] = float(string_list[i])
            except ValueError:
                equilibrium_frequencies = None
                break
            else:
                i += 1
        return r, equilibrium_frequencies


    @classmethod
    def decode_sequences(cls, tree):
        """Encode all sequences in a tree to bytes.

        Parameters
        ----------
        tree : ete3.PhyloTree
            A phylogenetic tree with all sequences simulated.
        """
        AA_table = reverse_AA_table()
        for node in tree.traverse('preorder'):
            if hasattr(node, 'sequence'):
                decoded = "".join([AA_table[i] for i in list(node.sequence)])
                node.add_feature('sequence', decoded)


    # 去除每个node中的profile与rates
    @classmethod
    def clean_up(cls, tree):
        for node in tree.traverse():
            node.del_feature('profile')
            node.del_feature('rates')


def _sample_sequences(p, n=None):
    if p.ndim == 2:
        n = p.shape[0]
    sequence = (numpy.random.random((n, 1)) > p.cumsum(axis=1)).sum(axis=1)
    return sequence.clip(max=(p.shape[-1] - 1)).astype('i1')


class HeterogeneousProteinSequenceSimulator(
    ProteinSimulatorMixin,
    BranchEvolutionSimulator,
):
    """A protein sequence simulator with heterogeneity"""
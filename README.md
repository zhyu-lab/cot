# CoT

Copy number Transformer for scDNA-seq data.

## Requirements

* Python 3.9+.

# Installation

## Clone repository

First, download CoT from github and change to the directory:

```bash
git clone https://github.com/zhyu-lab/cot
cd cot
```

## Create conda environment (optional)

Create a new environment named "cot":

```bash
conda create --name cot python=3.9
```

Then activate it:

```bash
conda activate cot
```

## Install requirements

```bash
python -m pip install -r ./tr_ae/requirements.txt
cd prep
cmake .
make
cd ..
chmod +x run_cot.sh ./hmm/run_CloneHMM.sh ./hmm/CloneHMM
```

# Usage

## Step 1: prepare input data

The “./prep/bin/prepInput” command is used to obtain read counts, GC-content and mappability data. 

To successfully run the command you will need to obtain/create these items:

* A merged BAM file (10X Genomics) containing sequencing data of all cells or seperate BAM files of single cells
* A FASTA file defining reference sequence
* A BIGWIG file for calculating mappability scores 
* A barcode file listing the barcodes of all cells or names of all BAM files to be analyzed

Reference sequence file formatted as .fasta can be downloaded from [UCSC genome browser](http://hgdownload.soe.ucsc.edu/downloads.html).

Mappability files formatted as .bw for human genomes are available from UCSC genome browser (https://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/ for hg18, https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/ for hg19, and https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/ for hg38). 

Users can also generate their own mappability files using [gem-library](https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/) and [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) utility.

Here is an example for creating mappability file from reference sequence 
(suppose “$gemlibrary” and “$bigwig” are the pathes where gem-library and wigToBigWig are installed, respectively).

```
chmod +x $gemlibrary/bin/gem* $bigwig/wigToBigWig
export PATH=$PATH:$gemlibrary/bin:$bigwig
gem-indexer -T 10 -c dna -i ./testData/refs/example.fa -o ./testData/refs/example_index
gem-mappability -T 10 -I ./testData/refs/example_index.gem -l 36 -o ./testData/refs/example_36
gem-2-wig -I ./testData/refs/example_index.gem -i ./testData/refs/example_36.mappability -o ./testData/refs/example_36
wigToBigWig ./testData/refs/example_36.wig ./testData/refs/example.sizes ./testData/refs/example_36.bw
```

The all arguments of the “./prep/bin/prepInput” command are as follows:

| Parameter     | Description                                                                                      | Possible values                                                                         |
| ------------- | ------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------- |
| -b, --bam     | a merged BAM file (10X Genomics) or a directory containing BAM files of the cells to be analyzed | Ex: /path/to/sample.bam                                                                 |
| -r, --ref     | genome reference file (.fasta)                                                                   | Ex: /path/to/hg19.fa                                                                    |
| -m, --map     | mappability file (.bw)                                                                           | Ex: /path/to/hg19.bw                                                                    |
| -B, --barcode | a file listing the barcodes of all cells or names of all BAM files to be analyzed                | Ex: /path/to/barcodes.txt                                                               |
| -c, --chrlist | the entries chromosomes to be analyzed (should be separated by commas)                           | Ex: chrlist=1,2,3,4,5  default:1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 |
| -s, --binsize | set the size of bin to count reads                                                               | Ex: binsize=200000  default:500000                                                      |
| -q, --mapQ    | threshold value for mapping quality                                                              | Ex: mapQ=10  default:20                                                                 |
| -o, --output  | output file to save results                                                                      | Ex: /path/to/example.txt                                                                |

Example:

```
./prep/bin/prepInput -b /path/to/sample.bam -r /path/to/hg19.fasta -m /path/to/hg19.bw -B /path/to/barcodes.txt -o example.txt
```

## Step 2: train the CoT model

The “./tr_ae/train.py” Python script is used to learn latent representations of cells and cluster cells into distinct subpopulations.

The arguments to run “./tr_ae/train.py” are as follows:

| Parameter    | Description                                            | Possible values                    |
| ------------ | ------------------------------------------------------ | ---------------------------------- |
| --input      | input file generated by “./prep/bin/prepInput” command | Ex: /path/to/example.txt           |
| --output     | a directory to save results                            | Ex: /path/to/results               |
| --epochs     | number of epoches to train the CoT                     | Ex: epochs=200  default:30         |
| --batch_size | batch size                                             | Ex: batch_size=32  default:32      |
| --lr         | learning rate                                          | Ex: lr=0.0005  default:0.0001      |
| --max_k      | maximum number of clusters to consider                 | Ex: max_k=30  default:30           |
| --latent_dim | latent dimension                                       | Ex: latent_dim=4  default:10       |
| --d_seg      | the dimension of model or the size of the data segment | Ex: d_seg=256  default:1024        |
| --cov_type   | the covariance type of each mixture component of GMM   | Ex: cov_type='tied' default:'full' |

Example:

```
python ./tr_ae/train.py --input ./data/example.txt --epochs 200 --batch_size 32 --lr 0.0001 --latent_dim 10 --d_seg 256 --seed 0 --max_k -1 --output data
```

```
python ./tr_ae/train.py --input ./data/10X_Genomics_dataset_C.txt --epochs 50 --batch_size 32 --lr 0.0001 --latent_dim 10 --d_seg 1024 --seed 0 --output data --cov_type='tied'
```

## Step 3: detect single-cell CNAs

The “./hmm/CloneHMM.m” MATLAB script is used to call single-cell CNAs. 

The arguments run “./hmm/CloneHMM.m” are as follows:

| Parameter | Description                                             | Possible values        |
| --------- | ------------------------------------------------------- | ---------------------- |
| inputFile | “lrc.txt” file generated by “./tr_ae/train.py” script   | Ex: /path/to/lrc.txt   |
| labelFile | “label.txt” file generated by “./tr_ae/train.py” script | Ex: /path/to/label.txt |
| outputDir | a directory to save results                             | Ex: /path/to/results   |
| maxCN     | maximum copy number to consider                         | Ex: 6  default:10      |

Example:

```
CloneHMM('../data/lrc.txt','../data/label.txt','../data',10)
plot_cn_results('../data')
```

**We also provide a script “run_cot.sh” to integrate all three steps to run CoT.**
This script requires that MATLAB Compiler Runtime (MCR) v91 (R2016b) is installed in user's machine. 
The MCR can be downloaded from [MathWorks Web site](https://www.mathworks.com/products/compiler/matlab-runtime.html). 

Example:

```
./run_cot.sh /path/to/bam /path/to/ref.fa /path/to/ref.bw /path/to/barcodes.txt /path/to/results
```

Type ./run_cot.sh to learn details about how to use this script.

## Citation

Please cite CoT in your publications if it helps your research:

``` bibtex
@article{liu2024cot,
  title={CoT: a transformer-based method for inferring tumor clonal copy number substructure from scDNA-seq data},
  author={Furui Liu, Fangyuan Shi, Fang Du, Xiangmei Cao and Zhenhua Yu},
  journal={Briefings in Bioinformatics, revision submitted},
  year={2024}
}
```

# Contact

If you have any questions, please contact lfr_nxu@163.com.
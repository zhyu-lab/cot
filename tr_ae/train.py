import argparse
import os
import datetime
import torch
import numpy as np
from matplotlib import pyplot as plt
from genomedata import GenomeData
from datasets import *
from statistics import mean
from transformer import *
from torch.nn import MSELoss
from clustering import GMM
from sklearn.decomposition import PCA


os.environ["CUDA_VISIBLE_DEVICES"] = "0"


def main(args):
    start_t = datetime.datetime.now()

    gd = GenomeData(args.input)
    gd.load_data()
    gd.preprocess_data()

    data = gd.data_rc_all.copy()
    data = np.transpose(data)
    data_bk = data.copy()

    def setup_seed(seed):
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        np.random.seed(seed)

    setup_seed(args.seed)
    train_set = get_seg_data(data_bk, args.d_seg)
    loader = Data.DataLoader(CellDataSet(train_set), args.batch_size, True)
    # define and create Transformer architecture
    set_tr_parameter(args.d_seg, train_set, args.latent_dim)
    model = Transformer().cuda()
    criterion = MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, betas=(0.9, 0.98), eps=1e-08, weight_decay=0,
                                 amsgrad=False)
    # Start training the model
    model.train()
    train_loss = []
    epochList = []
    for epoch in range(args.epochs):
        losses = []
        epochList.append(epoch)
        for enc_inputs in loader:
            ori_inputs = enc_inputs.cuda()
            enc_inputs = ori_inputs.clone()
            outputs = model(enc_inputs)

            # The loss is calculated by intercepting the unfilled data segment
            unfilled_outputs = ((outputs.clone()).reshape(args.batch_size, -1))[..., 0:(data_bk.shape[1])]
            unfilled_ori_inputs = ((ori_inputs.clone()).reshape(args.batch_size, -1))[..., 0:(data_bk.shape[1])]
            loss = criterion(unfilled_outputs, unfilled_ori_inputs)

            losses.append(loss.item())
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        mean_loss = mean(losses)
        train_loss.append(np.array(np.array(mean_loss / args.batch_size)))
        print("epo: " + str(epoch) + " loss:" + str(mean_loss / args.batch_size))

    output_dir = args.output
    if not os.path.isdir(output_dir):
        os.mkdir(args.output)
    ll_file = output_dir + '/loss.txt'
    if os.path.isfile(ll_file):
        os.remove(ll_file)
    file_o = open(ll_file, 'w')
    np.savetxt(file_o, np.c_[np.reshape(train_loss, (1, len(train_loss)))], fmt='%f', delimiter=',')
    file_o.close()

    # get latent representation of single cells after CoT training is completed
    model.eval()
    loader2 = Data.DataLoader(CellDataSet(train_set), args.batch_size, False)
    feature = None
    feature_List = []
    torch.cuda.empty_cache()
    for enc_inputs_data in loader2:
        ori_inputs = enc_inputs_data.cuda()
        enc_inputs = ori_inputs.clone()
        encodings = model.Encoder(enc_inputs)
        if feature is None:
            feature = encodings
        else:
            feature = torch.cat([feature, encodings], dim=0)

    feature_List = feature.cpu().detach().numpy()

    # use PCA to reduce the data dimension: from latent dimension to 3D
    pca = PCA(n_components=3)
    pca.fit(feature_List)
    feature_List = pca.transform(feature_List)

    # use Gaussian mixture model to cluster the single cells
    print('clustering the cells...')
    if args.max_k <= 0:
        max_k = np.max([1, feature_List.shape[0] // 5])
    else:
        max_k = np.min([args.max_k, feature_List.shape[0] // 5])
    label_p, number = GMM(feature_List, max_k).cluster(output_dir)
    print('inferred number of clusters: {}'.format(number))

    # save results
    lrc_file = output_dir + '/lrc.txt'
    if os.path.isfile(lrc_file):
        os.remove(lrc_file)
    file_o = open(lrc_file, 'a')
    np.savetxt(file_o, np.c_[gd.bin_size], fmt='%d', delimiter=',')
    np.savetxt(file_o, np.c_[np.reshape(gd.data_chr_all, (1, len(gd.data_chr_all)))], fmt='%d', delimiter=',')
    np.savetxt(file_o, np.c_[np.reshape(gd.data_bin_all, (1, len(gd.data_bin_all)))], fmt='%d', delimiter=',')
    np.savetxt(file_o, np.c_[np.transpose(gd.data_lrc_all)], fmt='%.3f', delimiter=',')
    file_o.close()

    label_file = output_dir + '/label.txt'
    file_o = open(label_file, 'w')
    np.savetxt(file_o, np.c_[np.reshape(gd.barcodes, (1, len(gd.barcodes)))], fmt='%s', delimiter=',')
    np.savetxt(file_o, np.c_[np.reshape(label_p, (1, len(label_p)))], fmt='%d', delimiter=',')
    file_o.close()

    latent_file = output_dir + '/latent.txt'
    file_o = open(latent_file, 'w')
    np.savetxt(file_o, np.c_[feature_List], fmt='%.3f', delimiter=',')
    file_o.close()

    end_t = datetime.datetime.now()
    print('elapsed time: ', (end_t - start_t).seconds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="CoT")
    parser.add_argument('--epochs', type=int, default=20, help='number of epoches to train the CoT.')
    parser.add_argument('--batch_size', type=int, default=32, help='batch size.')
    parser.add_argument('--lr', type=float, default=0.0001, help='learning rate.')
    parser.add_argument('--max_k', type=int, default=10, help='the maximum number of clusters to consider.')
    parser.add_argument('--latent_dim', type=int, default=10, help='the latent dimension.')
    parser.add_argument('--seed', type=int, default=0, help='random seed.')
    parser.add_argument('--input', type=str, default='', help='a file containing read counts, GC-content and mappability data.')
    parser.add_argument('--output', type=str, default='', help='a directory to save results.')
    parser.add_argument('--d_seg', type=int, default=512,
                        help='the dimension of model or the size of the data segment.')
    args = parser.parse_args()
    main(args)
